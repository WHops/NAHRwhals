#' Merge neighboring alignments in a PAF file for stitching together chunked alignments
#'
#' This function identifies alignments in a PAF file that are potential neighbors by checking if their distance in any direction (+-x, +-y) is less than 5% of their alignment length. Alignments are considered as neighbors if they have the same strand and the end of one is the start of the other. The function returns a data.frame of row pairs that are considered neighbors based on the above criteria. This function is called in the process of stitching together chunked alignments.
#'
#' @param inpaf A PAF file containing alignments
#' @param second_run A logical indicating if this is the second run of the stitching process
#' @param inparam_chunklen An optional parameter to specify the length of the input chunks
#' @param inparam_compression An optional parameter to specify the compression factor to be used in the merging process
#'
#' @return A data.frame of row pairs that are considered neighbors
#'
#' @author Wolfram Hoeps
#' @export
merge_paf_entries_intraloop <- function(inpaf, second_run = F, inparam_chunklen = NULL, inparam_compression = NULL) {
  
  
  # For safety: sort entries by qstart. Reset row names so they start at 1.
  #inpaf <- inpaf[order(inpaf$tstart), ]
  inpaf_rownames <- row.names(inpaf)

  rownames(inpaf) <- NULL

  # We consider alignments as 'potential neighbours' if their distance in any direction
  # (+-x, +-y) is less than 5% of their alignment length.
  # tolerance_bp = 10 #
  if (second_run) {
    # Not setting to 1 because i'm scared of edge cases.. But should be one.
    # tolerance_bp = 2

    # Fixing a bug here which caused problems with wrong stuff being merged,
    # If their start AND end points are both in reach of one point.
    chunklen <- min(inpaf$qlen)
    tolerance_bp_unlimited <- 0.5 * inparam_compression
    tolerance_bp <- min(tolerance_bp_unlimited, chunklen * 0.9)
  } else {
    # On second thought, let's have this cap on the length of the input chunklen.
    tolerance_bp <- 0.05 * inparam_chunklen
  }
  # Identify alignments that border each other: same strand, and end of of is the start
  # of the other. With some tolerance
  # Take only one half of the minus matrix so pairs dont appear twice.
  inpaf_t <- inpaf
  # inpaf_t =  transform(
  #   inpaf_t2,
  #   qend = ifelse(strand == '-', qstart, qend),
  #   qstart = ifelse(strand == '-', qend, qstart)
  # )








  rowpairs <- data.frame(which(
    
    (
      (abs(outer(inpaf_t$tend, inpaf_t$tstart, "-")) < tolerance_bp) &
        (abs(outer(inpaf_t$qend, inpaf_t$qstart, "-")) < tolerance_bp) & 
        (outer(inpaf_t$strand, inpaf_t$strand, "=="))
    ),
    arr.ind = T
  ))
  
  # Remove pairs with very small alignments.
  if (tolerance_bp > 100) {
    rowpairs <- rowpairs[inpaf_t[rowpairs$row, ]$alen > tolerance_bp, ]
    rowpairs <- rowpairs[inpaf_t[rowpairs$col, ]$alen > tolerance_bp, ]
  }

  rowpairs$row <- as.numeric(inpaf_rownames[rowpairs$row])
  rowpairs$col <- as.numeric(inpaf_rownames[rowpairs$col])




  return(rowpairs)
}
