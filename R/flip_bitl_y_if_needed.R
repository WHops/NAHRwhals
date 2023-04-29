#' flip_bitl_y_if_needed
#'
#' Determines whether a bitlocus is 'upside down', which happens when the
#' alignment is negative. To determine this, look only at values on the outer
#' 25% edge of the plot (assuming that SVs are mainly confined to the middle),
#' and determine if this outer part has more positive or negative alignments.
#' If there are more negative, then we flip the y-axis.
#' @param bitl a gridmatrix / bitlocus.
#' @return the same gridmatrix, flipped or unflipped.
#' @author  Wolfram Hoeps
#' @export
flip_bitl_y_if_needed <- function(bitl) {
  # Create a copy of the input dataframe.
  # We will mess with that copy, and trash it.
  bitl_f <- bitl



  # If the matrix is large enough, define a frame.
  if ((dim(bitl)[1] >= 4) | (dim(bitl)[2] >= 4)) {
    # But every dimension must be larger than 1.
    if ((dim(bitl)[1] >= 2) & (dim(bitl)[2] >= 2)) {
      # Define frame corners: Which values do we want to consider?
      # Minimum framesize is now 2. This is a bugfix.
      y_min <- max(round(dim(bitl_f)[1] * (1 / 4)), 2)
      y_max <- min(round(dim(bitl_f)[1] * (3 / 4)), dim(bitl_f)[1] - 2)

      x_min <- max(round(dim(bitl_f)[2] * (1 / 4)), 2)
      x_max <- min(round(dim(bitl_f)[2] * (3 / 4)), dim(bitl_f)[2] - 2)


      # W, 14th March 2022. Changing that part a bit.
      # Alternatively: set the non-corners of the frame to 0.
      bitl_f[y_min:y_max, ] <- 0
      bitl_f[, x_min:x_max] <- 0
    }
  }
  # # Set all inner values to zero. We don't care about them.
  # bitl_f[y_min:y_max, x_min:x_max] = 0



  # Count positive and negative alignments in the frame.
  pos_aln_sum <- sum(bitl_f[bitl_f > 0])
  neg_aln_sum <- abs(sum(bitl_f[bitl_f < 0]))

  # Emergency outclause. If all is zero, do not flip.
  if (pos_aln_sum + neg_aln_sum == 0) {
    pos_aln_sum <- 1
  }

  # If one is not at least double the other, we are unsure about the outcome.
  if ((max(pos_aln_sum, neg_aln_sum) / min(pos_aln_sum, neg_aln_sum)) < 2) {
    print("Warning. Frame-flipping determination potentially unsure. Consider using a wider interval. ")
    # Make an entry to the output logfile #
    if (exists("log_collection")) {
      log_collection$flip_unsure <<- T
    }
    # Log file entry done #
  }

  # Flip if needed (but on the bitl, not on the messed bitl_f.)
  if (neg_aln_sum > pos_aln_sum) {
    print("Inverse y axis detected. Flipping ... ")

    # If the matrix has only one row, we need a special case because apply
    # return a vector in this case, not a matrix.
    if (dim(bitl_f)[1] == 1) {
      bitl <- -bitl
    } else {
      bitl <- -apply(bitl, 2, rev)
    }
  }

  return(bitl)
}
