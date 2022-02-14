#' Find SV opportunities
#'
#' @description Building block: Take in a 'bitlocus',
#' and return a list of NAHR-driven rearrangements that
#' are possible on the reference sequence.
#'
#' @param sample  an nxm matrix (n=gridpoints_x, m=gridpoints_y)
#' @return A dataframe with columns 'p1', 'p2', 'sv'
#'
#' @author Wolfram Höps
#' @rdname find_sv_opportunities
#' @export
find_sv_opportunities <- function(sample) {
  double_SD_rows = which(rowSums(sample != 0) >= 2)
  all_opportunities = data.frame(matrix(ncol = 3, nrow = 0))
  
  
  # Rather disgusing block of if, else, ....
  # Probably can be implemented way more efficient.
  
  # For every row of the dot matrix, check first if this row
  # has more than one non-zero value in the matrix.
  for (nrow in 1:dim(sample)[1]) {
    if (nrow %in% double_SD_rows) {
      # ... If if does, now find the colnames of all those.
      identical_columns = as.numeric(which(sample[nrow,] != 0))
      
      # make a matrix containing every combination of identical columns.
      # Imagine we have one segment repeated in columns 4, 5 and 9.
      # Then we get NAHR possibilities in 4-5, 4-9 and 5-9.
      # If this is only one pair, then the one pair stays the same.
      combs = combn(identical_columns, 2)
      
      for (ncomb in (1:dim(combs)[2])) {
        identical_column_pair = combs[, ncomb]
        
        # Determine relative orientation
        orientation = sign(min(sample[nrow, identical_column_pair[1]], sample[nrow, identical_column_pair[2]]))
        
        result_pair = c(identical_column_pair, orientation)
        all_opportunities = rbind(all_opportunities, result_pair)
        
      }
    }
  }
  
  
  colnames(all_opportunities) = c('p1', 'p2', 'sv')
  all_opportunities = unique(all_opportunities)
  
  
  ## INV or DEL/DUP pair? ##
  # If at least one of the pair is negative, it's an 'inv'.
  if (dim(all_opportunities[all_opportunities$sv == -1,])[1] > 0) {
    all_opportunities[all_opportunities$sv == -1 ,]$sv = 'inv'
  }
  # If at least one of the pair is positive, it's a 'del+dup' (???)
  if (dim(all_opportunities[all_opportunities$sv == 1,])[1] > 0) {
    all_opportunities[all_opportunities$sv == 1 ,]$sv = 'del'
    # Duplicate plus pairs for dup and del.
    pluspairs = all_opportunities[all_opportunities$sv == 'del',]
    pluspairs$sv = 'dup'
    all_opportunities = rbind(all_opportunities, pluspairs)
  }
  
  
  return(all_opportunities)
  
}


# ----- Define a function for plotting a matrix ----- #
#' myImagePlot
#'
#' @description Found on stackexchange
#'
#' @param x a matrix
#' @return none, but image is plotted
#'
#' @author unknown (http://www.phaget4.org/R/image_matrix.html)
#' @rdname myImagePlot
#' @export
myImagePlot <- function(x, ...) {
  min <- min(x)
  max <- max(x)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <- c()
  # check for additional function arguments
  if (length(list(...))) {
    Lst <- list(...)
    if (!is.null(Lst$zlim)) {
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if (!is.null(Lst$yLabels)) {
      yLabels <- c(Lst$yLabels)
    }
    if (!is.null(Lst$xLabels)) {
      xLabels <- c(Lst$xLabels)
    }
    if (!is.null(Lst$title)) {
      title <- Lst$title
    }
  }
  # check for null values
  if (is.null(xLabels)) {
    xLabels <- c(1:ncol(x))
  }
  if (is.null(yLabels)) {
    yLabels <- c(1:nrow(x))
  }
  
  layout(matrix(
    data = c(1, 2),
    nrow = 1,
    ncol = 2
  ),
  widths = c(4, 1),
  heights = c(1, 1))
  
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- rgb(seq(0, 1, length = 256),
                   # Red
                   seq(0, 1, length = 256),
                   # Green
                   seq(1, 0, length = 256))  # Blue
  ColorLevels <- seq(min, max, length = length(ColorRamp))
  
  # # Reverse Y axis
  # reverse <- nrow(x) : 1
  # yLabels <- yLabels[reverse]
  # x <- x[reverse,]
  
  # Data Map
  par(mar = c(3, 5, 2.5, 2))
  image(
    1:length(xLabels),
    1:length(yLabels),
    t(x),
    col = ColorRamp,
    xlab = "",
    ylab = "",
    axes = FALSE,
    zlim = c(min, max)
  )
  if (!is.null(title)) {
    title(main = title)
  }
  axis(
    BELOW <- 1,
    at = 1:length(xLabels),
    labels = xLabels,
    cex.axis = 0.7
  )
  axis(
    LEFT <- 2,
    at = 1:length(yLabels),
    labels = yLabels,
    las = HORIZONTAL <- 1,
    cex.axis = 0.7
  )
  
  # Color Scale
  par(mar = c(3, 2.5, 2.5, 2))
  image(
    1,
    ColorLevels,
    matrix(
      data = ColorLevels,
      ncol = length(ColorLevels),
      nrow = 1
    ),
    col = ColorRamp,
    xlab = "",
    ylab = "",
    xaxt = "n"
  )
  
  layout(1)
}
# ----- END plot function ----- #


#' Find SV carry_out_compressed_sv
#'
#' @description TBD
#'
#' @param bitl  an nxm matrix (n=gridpoints_x, m=gridpoints_y)
#' @param input_ins a vector with instruction: p1, p2, sv.
#'
#' @author Wolfram Höps
#' @rdname find_sv_opportunities
#' @export
carry_out_compressed_sv <- function(bitl, input_ins) {
  pair = input_ins[1:2]
  action = input_ins[3]
  
  # Modify the matrix accordingly
  if (action == 'del') {
    bitl_mut = bitl[,-c(as.numeric(pair[1]):(as.numeric(pair[2]) - 1))]
  } else if (action == 'dup') {
    bitl_mut = cbind(bitl[, 1:as.numeric(pair[2])],
                     bitl[, as.numeric(pair[1] + 1):dim(bitl)[2]])
  } else if (action == 'inv') {
    bitl_mut = cbind(cbind(bitl[, 1:as.numeric(pair[1])], -bitl[, (as.numeric(pair[2]) -
                                                                     1):(as.numeric(pair[1]) + 1)]),
                     bitl[, as.numeric(pair[2]):dim(bitl)[2]])
  }
  
  
  
}



#' From Stackoverflow
#' @export
colMax <- function(data)
  sapply(data, max, na.rm = TRUE)


#' eval_mutated_seq
#'
#' @description Evaluation function: how similar are the input and output?
#' This needs work still...
#'
#' @param bitlocus  an nxm matrix (n=gridpoints_x, m=gridpoints_y)
#' @return numeric [0-1] with evaluation score.
#'
#' @author Wolfram Höps
#' @rdname eval_mutated_seq
#' @export
#' @export
eval_mutated_seq <- function(bitlocus) {
  # Old stuff here:
  #coverage_x =  sum(colMax(as.data.frame(bitlocus)) > 0) / dim(bitlocus)[2]
  #coverage_y =  sum(colMax(as.data.frame(t(bitlocus))) > 0) / dim(bitlocus)[1]
  
  
  # A matrix is perfect if it has a non-zero, non-negative diagonal.
  # We also want the matrix to be symmetrical.
  
  # Symmetry needs some work. It's not reciprocal (2-1 != 1-2).
  symmetry = -(abs((dim(bitlocus)[2] / dim(bitlocus)[1]) - 1)) + 1
  diag_filled = sum(diag(bitlocus) > 0) / min(dim(bitlocus))
  
  return(round(symmetry * diag_filled, 3))
}


#' gimme_sample_matrix
#'
#' @description Think of this as a function used for development and debugging.
#' @param mode diff/same: should x and y be the same sequence?
#' @return matrix (bitlocus)
#'
#' @author Wolfram Höps
#' @rdname gimme_sample_matrix
#' @export
gimme_sample_matrix <- function(mode = 'diff') {
  #samplefasta_link = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/vignettes/simulated_seq_10kb_4SDs.fa'
  #samplemutfasta_link = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/vignettes/simulated_seq_10kb_del_trim.fa'
  samplefasta_link = system.file('extdata', 'simulated_seq_10kb_4SDs.fa', package =
                                   'nahrtoolkit')
  samplemutfasta_link = system.file('extdata', 'simulated_seq_10kb_dup.fa', package =
                                      'nahrtoolkit')
  
  #samplemutfasta_link = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/vignettes/simulated_seq_10kb_del.fa'
  
  if (mode == 'same') {
    samplemutfasta_link = samplefasta_link
  }
  
  samplepaf_link = paste0('blub33.paf')
  make_chunked_minimap_alnment(
    samplefasta_link,
    samplemutfasta_link,
    samplepaf_link,
    chunklen = 500,
    minsdlen = 0,
    saveplot = F,
    hllink = F,
    hltype = F,
    quadrantsize = 10000
  )
  grid = wrapper_paf_to_bitlocus(samplepaf_link, compression = 100, minlen = 0)[[3]]
  
  
  x_missing = which(min(grid$x):max(grid$x) %in% grid$x == F)
  y_missing = which(min(grid$y):max(grid$y) %in% grid$y == F)
  
  for (xm in x_missing) {
    grid = rbind(grid, c(xm, xm, 0))
  }
  for (ym in y_missing) {
    grid = rbind(grid, c(ym, ym, 0))
  }
  
  
  sample = reshape2::dcast(grid, y ~ x, fill = 0)
  sample$x = NULL
  sample$y = NULL
  sample = as.matrix(sample)
  
  colnames(sample) = 1:dim(sample)[2]
  rownames(sample) = 1:dim(sample)[1]
  
  return(sample)
}

#' ### UNDER ACTIVE DEVELOPMENT ###
#' explore_mutation_space
#'
#' @description Main workhorse, tying together the pieces.
#' @param bitlocus matrix, nxm
#' @param depth How many SVs in sequence should be simulated?
#' @return evaluation matrix
#'
#' @author Wolfram Höps
#' @rdname gimme_sample_matrix
#' @export
explore_mutation_space <- function(bitlocus, depth) {
  #sample = gimme_sample_matrix()
  sample = bitlocus
  pairs = find_sv_opportunities(sample)
  
  count = 1
  depth = 3
  
  res = data.frame(matrix(ncol = depth + 1, nrow = 0))
  res = rbind(res, unlist(c(
    eval_mutated_seq(bitlocus), 'ref', rep('NA', depth - 1)
  )))
  colnames(res) = c('eval', paste0('mut', 1:depth))
  
  
  
  for (npair_level1 in 1:dim(pairs)[1]) {
    print(npair_level1)
    pair_level1 = pairs[npair_level1,]
    bitl_mut = carry_out_compressed_sv(sample, pair_level1)
    
    # add
    res_add = unlist(c(
      eval_mutated_seq(bitl_mut),
      paste(pair_level1$p1, pair_level1$p2, pair_level1$sv, sep = '_')
    ))
    names(res_add) = c('eval', 'mut1')
    res = dplyr::bind_rows(res, res_add)
    if (depth > 1) {
      newpairs = find_sv_opportunities(bitl_mut)
      for (npair_level2 in 1:dim(newpairs)[1]) {
        pair_level2 = newpairs[npair_level2,]
        bitl_mut2 = carry_out_compressed_sv(bitl_mut, pair_level2)
        
        # add
        res_add = unlist(c(
          eval_mutated_seq(bitl_mut2),
          paste(pair_level1$p1, pair_level1$p2, pair_level1$sv, sep = '_'),
          paste(pair_level2$p1, pair_level2$p2, pair_level2$sv, sep =
                  '_')
        ))
        names(res_add) = c('eval', 'mut1', 'mut2')
        res = dplyr::bind_rows(res, res_add)
        
      }
      
      if (depth > 2) {
        newpairs3 = find_sv_opportunities(bitl_mut2)
        for (npair_level3 in 1:dim(newpairs3)[1]) {
          pair_level3 = newpairs3[npair_level3,]
          bitl_mut3 = carry_out_compressed_sv(bitl_mut2, pair_level3)
          
          # add
          res_add = unlist(c(
            eval_mutated_seq(bitl_mut3),
            paste(
              pair_level1$p1,
              pair_level1$p2,
              pair_level1$sv,
              sep = '_'
            ),
            paste(
              pair_level2$p1,
              pair_level2$p2,
              pair_level2$sv,
              sep = '_'
            ),
            paste(
              pair_level3$p1,
              pair_level3$p2,
              pair_level3$sv,
              sep = '_'
            )
          ))
          names(res_add) = c('eval', 'mut1', 'mut2', 'mut3')
          res = dplyr::bind_rows(res, res_add)
        }
      }
    }
  }
  return(res)
}


# Sample case
#
#C = gimme_sample_matrix(mode = 'diff')
#a = explore_mutation_space(C)
# eval_mutated_seq(C)
#
# samplei = ginv(sample)
#
# sample %*% samplei
#
# myImagePlot(sample)
# myImagePlot(bitl_mut)
#
# myImagePlot(sample %*% t(sample))
# library(MASS)
#
# s2 = sample
# s3 = cbind(cbind(s2[, 1:4], s2[, 4]), s2[, 5:14])
# colnames(s3) = 1:dim(s3)[2]
# myImagePlot(s3)
#
# find_sv_opportunities(s2)
#
#
#
#
#
# A = rbind(c(2, 3), c(4, 5))
# inv(A)
# solve(A)
#
# Ai = MASS::ginv(A)
#
# A %*% Ai
