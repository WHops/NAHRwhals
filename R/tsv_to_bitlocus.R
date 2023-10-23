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

    # Filter alignments by length
    if (nrow(paf) > 10000){
      minlen = minlen * 2
    }
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
    paf <- compress_paf_fnct(
      inpaf_df = paf,
      save_outpaf = F,
      second_run = T,
      inparam_compression = compression,
      inparam_chunklen = inparam_chunklen
    )



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


#' Calculate the "bounced" points of a vector given a point
#'
#' This function calculates the "bounced" points of a vector given a point, where the bounced points are the points where the vector intersects the x- and y-axes.
#' The function returns a data frame containing the bounced points.
#'
#' @param vectors A data frame containing the start and end points of each vector, as well as slope and intercept information.
#' @param point A numeric vector containing the x- and y-coordinates of the point.
#'
#' @return A data frame containing the bounced points of the vector.
#'
#' @author Wolfram Hoeps
#' @export
bounce_point <- function(vectors, point) {
  newpoints <- data.frame(matrix(ncol = 2, nrow = 0))
  newpoints[1, ] <- point
  colnames(newpoints) <- c("x", "y")

  if (!(point[1] %in% x_orig_visited)) {
    # Check for vertical overlaps
    y_overlap_vecs <- vectors[((vectors$tstart < as.numeric(point[1])) &
      (vectors$tend > as.numeric(point[1]))), ]
    # If there is overlap, calculate overlap of point with every overlapper
    if (dim(y_overlap_vecs)[1] > 0) {
      for (i in 1:dim(y_overlap_vecs)[1]) {
        newpoints <- rbind(
          newpoints,
          data.frame(x = point[1], y = ((y_overlap_vecs$slope) * as.numeric(point[1])) + y_overlap_vecs$y_intercept)
        )
      }
    }

    x_orig_visited <<- c(x_orig_visited, point[1])
  }

  if (!(point[2] %in% y_orig_visited)) {
    # Check for horizontal overlaps
    x_overlap_vecs <- vectors[((vectors$qlow < as.numeric(point[2])) &
      (vectors$qhigh > as.numeric(point[2]))), ]
    # If there is overlap, calculate overlap of point with every overlapper
    if (dim(x_overlap_vecs)[1] > 0) {
      for (i in 1:dim(x_overlap_vecs)[1]) {
        newpoints <- rbind(
          newpoints,
          data.frame(x = (point[2] - x_overlap_vecs$y_intercept) / x_overlap_vecs$slope, y = point[2])
        )
      }
    }

    y_orig_visited <<- c(y_orig_visited, point[2])
  }

  return(newpoints)
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
    # Load paf. If compression is too low change it automatically.
    # Compression is good if
    # 1) there are less than 50 alignments, and
    # 2) There are no incongruencies in the grid.

    # if (compression_params$auto_find == T) {
    #  wiggled_params = find_minlen_compression_params_wiggle(
    #    inpaf, compression_params
    #  )
    #  minlen = wiggled_params[[1]]
    #  compression = wiggled_params[[2]]
    # } else {
    print("Minlen/Compression manually chosen. Testing viability")
    minlen <- params$minlen
    compression <- params$compression
    # }

    print("Making the final grid with:")
    print(paste0("Minlen: ", minlen))
    print(paste0("Compression: ", compression))

    paf <- matrix(NA, nrow = 1e3, ncol = 1e3)
    gridlines.x <- rep(NA, 1e4)
    gridlines.y <- rep(NA, 1e4)
    # SUPER weird construct here. Should be re-written.
    # I want to re-run stuff until:
    # 1) The max number of alignments is kept.
    # 2) The max size of the grid is kept.
    while (((length(gridlines.x) + length(gridlines.y)) > params$max_size_col_plus_rows) | ((dim(paf)[1]) > params$max_n_alns)) {
      # Stop here and find out why with comp/minlen 200, the big piece gets lost.
      paf <- load_and_prep_paf_for_gridplot_transform(inpaf, minlen, compression)

      if ((dim(paf)[1]) > params$max_n_alns) {
        print("Too many alignments. Increasing minlen conditions")
        minlen <- minlen * 2
        compression <- compression * 2


        next()
      }

      if (is_cluttered_paf(paf)) {
        print("This paf file seems very cluttered and/or repetitive. Abort mission. ")
        return()
      }
      print(paste0("Leading to a paf of dimensions: ", dim(paf)[1]))
      gxy <- make_xy_grid(paf, n_additional_bounces = 50)

      # Make an entry to the output logfile #
      if (exists("log_collection")) {
        log_collection$compression <<- compression
      }

      gridlines.x <- gxy[[1]]
      gridlines.y <- gxy[[2]]

      print(paste0("Gridline dimensions: ", length(gridlines.x), " and ", length(gridlines.y)))

      if ((length(gridlines.x) + length(gridlines.y)) > params$max_size_col_plus_rows) {
        print("Too large grids. Repeating.")
        minlen <- minlen * 2
        compression <- compression * 2
        next()
      }
    }

    # Run bressi to fill the grid
    # 'bressi' is the bresenham algorithm. It can also handle vectors
    # with slope != 1. This used to be common in a previous version of ntk.
    # Right now, bressi is a bit overkill but we keep it in anyway.
    grid_list <- data.frame()

    for (i in 1:dim(paf)[1]) {
      grid_list <- rbind(grid_list, as.data.frame(
        bresenham(
          x = as.numeric(paf[i, c("tstart", "tend")]),
          y = as.numeric(paf[i, c("qstart", "qend")]),
          gridlines.x,
          gridlines.y,
          debug = F
        )
      ))
    }


    grid_list <- remove_duplicates_triple(grid_list)

    # grid_list_bu = grid_list
    # Sweep clean (?)
    # grid_list = clean_sweep_matrix(grid_list)

    # Sort by x
    grid_list <- grid_list[order(grid_list$x), ]
    # grid_list_bu = grid_list_bu[order(grid_list_bu$x),]

    grid_list <- remove_duplicates_triple(grid_list)
    # grid_list_bu = remove_duplicates_triple(grid_list_bu)

    ### MAKE AND SAVE PLOTS
    pregridplot_paf <- plot_paf(paf, gridlines.x, gridlines.y, linewidth = 0.1)
    # pp  + ggplot2::coord_cartesian(ylim=c(1e5,2e5), xlim = c(2.5e5, 3.5e5)) +
    #   ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())
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
    # browser()
    # p = plot_matrix_ggplot_named(grid_list, gridlines.x, gridlines.y)
    # print(p)
    # if (gridplot_save == F) {
    #   print(p)
    # } else {
    #   ggplot2::ggsave(
    #     p,
    #     file = gridplot_save,
    #     height = 15,
    #     width = 15,
    #     units = 'cm',
    #     device = 'pdf'
    #   )
    # }


    return(list(gridlines.x, gridlines.y, grid_list))
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
#' @param n_additional_bounces The number of additional times the points are bounced.
#'
#' @return A list containing the x and y coordinates of the grid, and a logical indicating whether the grid was successfully generated or not.
#'
#' @author Wolfram Hoeps
#' @export
make_xy_grid <- function(paf, n_additional_bounces = 10) {
  grid_successful <- T
  # Points is 'gridpoints'
  points_overall <- data.frame()

  # Bounce points a couple of times. For the bouncing algorithm,
  # see separate description.

  x_orig_visited <<- c()
  y_orig_visited <<- c()

  for (i in 1:dim(paf)[1]) {
    # Bounce start once
    # p1_bounced = bounce_point(paf, c(paf[i,]$tstart, paf[i,]$qlow))
    # # Bounce end once
    # p2_bounced = bounce_point(paf, c(paf[i,]$tend, paf[i,]$qhigh))

    # Bounce start once
    p1_bounced <- bounce_point(paf, c(paf[i, ]$tstart, paf[i, ]$qstart))
    # Bounce end once
    p2_bounced <- bounce_point(paf, c(paf[i, ]$tend, paf[i, ]$qend))
    points_overall <- rbind(rbind(points_overall, p1_bounced), p2_bounced)

    # Bounce bounced #1
    for (j in 1:dim(p1_bounced)[1]) {
      points_overall <- rbind(points_overall, bounce_point(paf, as.numeric(p1_bounced[j, ])))
    }
    # Bounce bounced #2
    for (k in 1:dim(p2_bounced)[1]) {
      points_overall <- rbind(points_overall, bounce_point(paf, as.numeric(p2_bounced[k, ])))
    }

    points_overall <- unique(points_overall[order(points_overall$x), ])
  }

  # Unique and sort
  points_overall <- unique(points_overall[order(points_overall$x), ])
  plot_helper_debug(paf, points_overall)
  # browser()

  points_last_iteration <- points_overall
  # # bounce some more
  for (j in (1:n_additional_bounces)) {
    print(paste0("Additional bounce: ", j, " out of ", n_additional_bounces))

    # If we have converged after 1st iteration, exit here.
    if ((j == 2) & (all(points_overall %in% points_last_iteration))) {
      print("Grid has converged. All fine.")
      # Make an entry to the output logfile #
      if (exists("log_collection")) {
        log_collection$grid_inconsistency <<- F
      }

      gridlines.y <- unique(round(sort(c(0, points_overall$y))))
      gridlines.x <- unique(round(sort(c(0, points_overall$x))))

      return(list(gridlines.x, gridlines.y, T))
    }

    if ((j > 2) & (all(points_overall %in% points_last_iteration))) {
      print(paste0("Grid has converged after ", j - 1, " additional bounces."))
      grid_successful <- F
      # Make an entry to the output logfile #
      if (exists("log_collection")) {
        log_collection$grid_inconsistency <<- T
      }

      plot_helper_debug(paf, points_overall)

      gridlines.y <- unique(round(sort(c(0, points_overall$y))))
      gridlines.x <- unique(round(sort(c(0, points_overall$x))))

      return(list(gridlines.x, gridlines.y, F))
    }

    # Run one more bounding point.
    if (j == 1) {
      points_for_consideration <- points_overall
    } else {
      points_for_consideration <- dplyr::anti_join(points_overall, points_last_iteration)
    }
    points_last_iteration <- points_overall
    # browser()

    # Bounce every point once
    # [note: we could def constrict this to *new* points (last generation)]
    for (i in 1:nrow(points_for_consideration)) {
      points_overall <- unique(rbind(points_overall, bounce_point(paf, as.numeric(points_for_consideration[i, ]))))
    }
  }

  print(
    "WARNING. Grid has not converged. Small alignment incongruencies are likely, results may not be reliable. Consider re-running with larger conversion-factor (typically >100)"
  )

  if (exists("log_collection")) {
    log_collection$grid_inconsistency <<- F
  }

  gridlines.y <- unique(round(sort(c(0, points_overall$y))))
  gridlines.x <- unique(round(sort(c(0, points_overall$x))))

  return(list(gridlines.x, gridlines.y, F))
}

#' Debugging function for plotting points and PAF data
#'
#' This function takes a data frame of alignments in .paf-like format, and a data frame
#' of points, and generates a plot to help with debugging.
#'
#' @param paf A data frame of alignments in .paf-like format.
#' @param points_overall A data frame of points.
#'
#' @export
plot_helper_debug <- function(paf, points_overall) {
  pp <- plot_paf(
    paf,
    unique(round(sort(c(0, points_overall$x)))),
    unique(round(sort(c(0, points_overall$y)))),
    linewidth = 0.5
  )
  print(pp)
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
