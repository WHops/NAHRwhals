
#' flip_bitl_y_if_needed
#' 
#' Determines whether a bitlocus is 'upside down', which happens when the
#' alignment is negative. To determine this, look only at values on the outer
#' 25% edge of the plot (assuming that SVs are mainly confined to the middle),
#' and determine if this outer part has more positive or negative alignments. 
#' If there are more negative, then we flip the y-axis. 
#' @param bitl_f a gridmatrix / bitlocus. 
#' @return the same gridmatrix, flipped or unflipped. 
#' @author  Wolfram Hoeps
#' @export
flip_bitl_y_if_needed <- function(bitl){
  
  # Create a copy of the input dataframe. 
  # We will mess with that copy, and trash it.
  bitl_f = bitl
  
  # Define frame corners: Which values do we want to consider?
  y_min = round(dim(bitl_f)[1] * (1/4))
  y_max = round(dim(bitl_f)[1] * (3/4))
  
  x_min = round(dim(bitl_f)[2] * (1/4))
  x_max = round(dim(bitl_f)[2] * (3/4))
  
  # Set all other values to zero. We don't care about them. 
  bitl_f[y_min:y_max, x_min:x_max] = 0
  
  # Count positive and negative alignments in the frame. 
  pos_aln_sum = sum(bitl_f[bitl_f > 0] )
  neg_aln_sum = abs(sum(bitl_f[bitl_f < 0] ))
  
  # If one is not at least double the other, we are unsure about the outcome. 
  if ((max(pos_aln_sum, neg_aln_sum) / min(pos_aln_sum, neg_aln_sum)) < 2){
    print("Warning. Frame-flipping determination potentially unsure. Consider using a wider interval. ")
  }
  
  # Flip if needed (but on the bitl, not on the messed bitl_f.)
  if (neg_aln_sum > pos_aln_sum){
    print('Inverse y axis detected. Flipping ... ')
    bitl = -bitl[nrow(bitl):1,]
  }
  
  return(bitl)
}