#' Load and Prepare PAF File for Gridplot Transformation
#'
#' This function loads and prepares the input PAF file. The preparations include filtering alignments by length, rounding start and endpoint to the nearest multiple of the compression factor, enforcing slope = 1 for all alignments, and transforming start/end in negative alignments. Finally, it adds qlow and qhigh columns to the PAF file.
#'
#' @param inpaf Character/link: Path to the input PAF file.
#' @param minlen Numeric: Minimum length of alignments to keep.
#' @param compression Numeric: Compression factor to round the start and endpoint to the nearest multiple of the compression factor.
#' @param quadrantsize Numeric: Size of the quadrant used to group alignments.
#' @param inparam_chunklen Numeric: Length of sequence chunks.
#'
#' @return The modified PAF file.
#'
#' @author Wolfram Hoeps
#' @export
load_and_prep_paf_for_gridplot_transform <-
  function(inpaf, minlen, compression, quadrantsize = 1e5, inparam_chunklen = 10000) {
    paf <- read.table(inpaf)
    colnames(paf) <- c(
      "qname",
      "qlen",
      "qstart",
      "qend",
      "strand",
      "tname",
      "tlen",
      "tstart",
      "tend",
      "nmatch",
      "alen",
      "mapq"
    )
    
    ### PAF prep ###
    
    
    # In paf, start is always smaller than end. For our slope etc calculation, it will be easier
    # to change the order of start and end, if the orientation of the alignment is negative.
    paf <- transform(
      paf,
      qend = ifelse(strand == "-", qstart, qend),
      qstart = ifelse(strand == "-", qend, qstart)
    )

    # Merge before compression.
    paf <- compress_paf_fnct(
      inpaf_df = paf,
      save_outpaf = F,
      second_run = T,
      inparam_compression = compression,
      inparam_chunklen = inparam_chunklen
    )
    paf <- paf[paf$alen > minlen, ]
    
    if (dim(paf)[1] == 0) {
      return(paf)
    }
    # (Error out if no paf survives initial filtering)
    stopifnot("Error: no alignment longer than minlen" = dim(paf)[1] > 0)
    
    # Round start/end by compression factor
    paf[, c("qstart", "qend", "tstart", "tend")] <- round(data.frame(paf[, c("qstart", "qend", "tstart", "tend")] / compression), 0) * compression
    
    # In paf, start is always smaller than end. For our slope etc calculation, it will be easier
    # to change the order of start and end, if the orientation of the alignment is negative.
    paf <- transform(
      paf,
      qend = ifelse(strand == "-", qstart, qend),
      qstart = ifelse(strand == "-", qend, qstart)
    )
    
    
    # If any entry is not slope 1 yet (which is unlikely after compression), then make it so.
    paf <- enforce_slope_one(paf)
    
    # In paf, start is always smaller than end. For our slope etc calculation, it will be easier
    # to change the order of start and end, if the orientation of the alignment is negative.
    paf <- transform(
      paf,
      qend = ifelse(strand == "-", qstart, qend),
      qstart = ifelse(strand == "-", qend, qstart)
    )

    # Now, in case any two ends have come close together, re-merge them.
    # paf <- compress_paf_fnct(
    #   inpaf_df = paf,
    #   save_outpaf = F,
    #   second_run = T,
    #   inparam_compression = compression,
    #   inparam_chunklen = inparam_chunklen
    # )
    
    
    # Add slope information to the paf.
    paf <- cbind(paf, add_slope_intercept_info(paf))
    
    # The Transform from earlier shoots back at us here. We now have to
    # make a new columns for qlow and qhigh - basically what was start and
    # end originally.
    paf$qlow <- -(apply(-paf[, c("qstart", "qend")], 1, function(x) {
      max(x)
    }))
    paf$qhigh <- apply(paf[, c("qstart", "qend")], 1, function(x) {
      max(x)
    })
    
    ### PAF prep over ###
    return(paf)
  }

#' @export
update_params <- function(minlen, compression) {
  next_power_of_10 <- 10^(floor(log10(minlen)) + 1)
  
  # If minlen is already at a power of 10 or has surpassed it by doubling, double the power of 10 value
  if (minlen >= next_power_of_10) {
    new_val <- next_power_of_10 * 2
  } else {
    # Otherwise, simply double the current value
    new_val <- minlen * 2
  }
  
  # Ensure that new_val does not exceed the next power of 10, if it does, set it to the next power of 10
  if (new_val >= next_power_of_10) {
    new_val <- next_power_of_10
  }
  
  # Return updated values
  list(minlen = new_val, compression = new_val)
}


#' Generate a bitlocus plot from a PAF file.
#'
#' This function takes a PAF file, and generates a bitlocus plot from it.
#' The generated bitlocus plot is based on a grid of x and y coordinates that are
#' calculated from the PAF file. The x and y coordinates of each alignment in the PAF
#' file are then calculated, and the corresponding positions in the grid are marked
#' with a "1". The resulting grid is then plotted as a bitlocus plot.
#'
#' @param inpaf A PAF file to be plotted.
#' @param params A list of parameters used to control the compression of the bitlocus plot.
#' @param gridplot_save If set to TRUE, the resulting bitlocus plot is saved to a file.
#' @param pregridplot_save If set to TRUE, the resulting pre-bitlocus plot is saved to a file.
#' @param pregridplot_nolines If set to TRUE, the resulting pre-bitlocus plot is displayed without grid lines.
#' @param saveplot If set to TRUE, the resulting bitlocus plot is saved to a file.
#'
#' @return A list containing the x and y coordinates of the grid, and the coordinates of the alignments.
#'
#' @author Wolfram Hoeps
#' @export
wrapper_paf_to_bitlocus <-
  function(inpaf,
           params,
           gridplot_save = F,
           pregridplot_save = F,
           pregridplot_nolines = F,
           saveplot = F) {
 
  minlen <- params$minlen
  compression <- params$compression

  # Initialize PAF matrix
  paf <- matrix(NA, nrow = 1e3, ncol = 1e3)
  gridlines.x <- rep(NA, 1e4)
  gridlines.y <- rep(NA, 1e4)

  while (TRUE) {
    paf <- load_and_prep_paf_for_gridplot_transform(inpaf, minlen, compression, inparam_chunklen=params$chunklen)
    
    # Check for max number of alignments
    if (dim(paf)[1] > params$max_n_alns) {
      minlen_compression <- update_params(minlen, compression)
      minlen = minlen_compression$minlen
      compression = minlen_compression$compression
      next()
    }

    # Check if PAF is cluttered
    if (is_cluttered_paf(paf)) {
      cat("This PAF file seems very cluttered and/or repetitive. Abort mission.\n")
      return()
    }

    # Create the x and y grid
    gxy <- make_xy_grid_fast(paf)
    gridlines.x <- gxy[[1]]
    gridlines.y <- gxy[[2]]

    #cat(paste0("Gridline dimensions: ", length(gridlines.x), " and ", length(gridlines.y), "\n"))
    
    # Check for max grid size
    if ((length(gridlines.x) + length(gridlines.y)) > params$max_size_col_plus_rows) {
      minlen_compression <- update_params(minlen, compression)
      minlen = minlen_compression$minlen
      compression = minlen_compression$compression
      next()
    }
    
    # Optional: Update log, if log_collection exists
    if (exists("log_collection")) {
      log_collection$compression <- compression
    }
        
    # Run bressi to fill the grid
    # 'bressi' is the bresenham algorithm. It can also handle vectors
    # with slope != 1. This used to be common in a previous version of ntk.
    # Right now, bressi is a bit overkill but we keep it in anyway.
    
    grid_list = bressiwrap(paf, gxy)
    grid_list <- grid_list[order(grid_list$x), ]
    
    ### MAKE AND SAVE PLOTS
    pregridplot_paf <- plot_paf(paf, gridlines.x, gridlines.y, linewidth = 0.1)
    
    n_pairs = nrow(find_sv_opportunities(gridlist_to_gridmatrix(list(gridlines.x, gridlines.y, grid_list))))
    if (n_pairs > 2000){
      print(paste0('Too many pairs: ', n_pairs))
      minlen_compression <- update_params(minlen, compression)
      minlen = minlen_compression$minlen
      compression = minlen_compression$compression
      next()
    }

    if (pregridplot_save == F) {
      print(pregridplot_paf)
    } else {
      ggplot2::ggsave(
        pregridplot_paf,
        file = pregridplot_save,
        height = 10,
        width = 10,
        units = "cm",
        device = "pdf"
      )
    }

    return(list(gridlines.x, gridlines.y, grid_list))
  }
    # If neither condition is met, break the loop
    break
  }

#' Ensure that the slope of each alignment is equal to one.
#'
#' This function takes a data frame of alignments, and ensures that the slope of each
#' alignment is equal to one. This is done by first calculating the length of the
#' alignment along the query and target sequences, and then setting the length along
#' the query and target sequences to be equal to the negative of the maximum of these
#' two lengths. The resulting alignment is then returned, with the length of the alignment
#' along the query and target sequences, and the squaresize, removed.
#'
#' @param df A data frame containing alignments.
#'
#' @return A modified data frame containing alignments, with the slope of each alignment equal to one.
#'
#' @author Wolfram Hoeps
#' @export
enforce_slope_one <- function(df) {
  df$qlen_aln <- abs(df$qend - df$qstart)
  df$tlen_aln <- abs(df$tend - df$tstart)
  
  df$squaresize <- -(apply(-df[, c("qlen_aln", "tlen_aln")], 1, function(x) {
    max(x)
  }))
  
  df$qend <- df$qstart + df$squaresize
  df$tend <- df$tstart + df$squaresize
  
  df[, c("qlen_aln", "tlen_aln", "squaresize")] <- NULL
  
  df <- df[(df$qstart != df$qend) & (df$tstart != df$tend), ]
  
  return(df)
}

#' Generate a grid of x and y coordinates for a bitlocus plot.
#'
#' This function takes a PAF file, and generates a grid of x and y coordinates for a
#' bitlocus plot. The grid is generated by first bouncing the start and end points of
#' each alignment in the PAF file, and then bouncing the resulting points a few more
#' times. The bounced points are then collected, and the unique points are returned
#' as the x and y coordinates of the grid.
#'
#' @param paf a dataframe of alignments which is in paf-like format.
#' @param generations. The max number of generations to run this thing. Default infinity.
#' @return A list containing the x and y coordinates of the grid, and a logical indicating whether the grid was successfully generated or not.
#'
#' @author Wolfram Hoeps
#' @export
make_xy_grid_fast <- function(paf_df, generations=100000) {
  
  if (exists("log_collection")) {
    log_collection$grid_inconsistency <<- T
  }
  
  all_x <- sort(unique(c(paf_df$tstart, paf_df$tend)))
  all_y <- sort(unique(c(paf_df$qstart, paf_df$qend)))
  
  slopes = (paf_df$qend - paf_df$qstart) / (paf_df$tend - paf_df$tstart)
  
  paf_tstarts = paf_df$tstart
  paf_tends = paf_df$tend
  paf_qstarts = paf_df$qstart
  paf_qends = paf_df$qend
  
  for (i in 1:generations) {
    #print(paste0('Generation ', i, '. n_xgridlines: ', length(all_x)))
    # Only process the most recent gridlines
    
    if (i == 1){
      process_x <- all_x
      process_y <- all_y
    } else {
      process_x <- new_x
      process_y <- new_y
    }
    
    # Pre-allocate
    max_size <- nrow(paf_df) * max(length(process_x), length(process_y))
    new_x <- numeric(max_size)
    new_y <- numeric(max_size)
    
    # Index counters
    idx_x <- 1
    idx_y <- 1
    
    for (row in 1:nrow(paf_df)) {
      
      slope = slopes[row] 
      paf_tstarts_row = paf_tstarts[row]
      paf_tends_row = paf_tends[row]
      paf_qstarts_row = paf_qstarts[row]
      paf_qends_row = paf_qends[row]
      
      current_xs = process_x[(process_x > paf_tstarts_row) & (process_x < paf_tends_row)]
      current_ys = process_y[(process_y > min(c(paf_qstarts_row, paf_qends_row))) & (process_y < max(c(paf_qstarts_row, paf_qends_row)))]
      
      # Identify overlaps between current_xs and pairwise alignments
      delta_x = current_xs - paf_tstarts_row
      new_y_vals = paf_qstarts_row + slope * delta_x
      n = length(new_y_vals)
      if (n > 0){
        new_y[idx_y:(idx_y + n - 1)] <- new_y_vals
        idx_y <- idx_y + n
      }
      
      # Identify overlaps between current_ys and pairwise alignments
      delta_y = current_ys - paf_qstarts_row
      new_x_vals = paf_tstarts_row + delta_y / slope
      n = length(new_x_vals)
      if (n > 0){
        new_x[idx_x:(idx_x + n - 1)] <- new_x_vals
        idx_x <- idx_x + n
      }
    }
    
    # Trim the vectors to the actual sizes
    new_x <- new_x[1:(idx_x-1)]
    new_y <- new_y[1:(idx_y-1)]
    
    # New_x and new_y should be subset to keep only those values which are not part of all_x and all_y
    new_x <- unique(new_x[!(new_x %in% all_x)])
    new_y <- unique(new_y[!(new_y %in% all_y)])
    
    prev_x <- all_x
    prev_y <- all_y
    
    all_x <- c(all_x, new_x)
    all_y <- c(all_y, new_y)
    
    # If no new gridlines are identified, break
    if (length(all_x) == length(prev_x) && length(all_y) == length(prev_y)) {
      #print(paste0('Grid has converged after ', i-1, ' iterations with ', length(all_x), ' lines'))
      
      if (exists("log_collection")) {
        log_collection$grid_inconsistency <<- F
      }
      
      break
    }
  }
  
  return(list(sort(unique(all_x)), sort(unique(all_y))))
}


#' Add slope and intercept information to a data frame of vectors
#'
#' This function calculates the slope, inverse slope, y-intercept, and x-intercept of a vector given its start and end points.
#' The resulting slope and intercept values are added to the input data frame as additional columns.
#'
#' @param vector_df A data frame containing the start and end points of each vector.
#' @param xstart The name of the column in `vector_df` containing the x-coordinate of the start point.
#' @param xend The name of the column in `vector_df` containing the x-coordinate of the end point.
#' @param ystart The name of the column in `vector_df` containing the y-coordinate of the start point.
#' @param yend The name of the column in `vector_df` containing the y-coordinate of the end point.
#'
#' @return A data frame containing the slope, inverse slope, y-intercept, and x-intercept of each vector.
#'
#' @author Wolfram Hoeps
#' @export
add_slope_intercept_info <- function(vector_df, xstart = "tstart", xend = "tend", ystart = "qstart", yend = "qend") {
  dx <- vector_df[xend] - vector_df[xstart]
  dy <- vector_df[yend] - vector_df[ystart]
  
  # Slope
  m <- dy / dx
  
  # Inverse slope
  m_inv <- 1 / m
  
  # y-intercept
  y_intercept <- vector_df[ystart] - (m * vector_df[xstart])
  
  # x-intercept
  x_intercept <- y_intercept / m
  
  slope_df <- data.frame(m, m_inv, y_intercept, x_intercept)
  colnames(slope_df) <- c("slope", "slope_inv", "y_intercept", "x_intercept")
  
  return(slope_df)
}

#   bounce_point(paf, c(points_current_gen[i,'t'], points_current_gen[i,'q']))