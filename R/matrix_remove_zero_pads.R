#" matrix_remove_zero_pads
#' Quick and dirty way of getting rid of padding zero columns or rows. 
#' @author Wolfram Häps
#' @export
matrix_remove_zero_pads <- function(bitl_f){
   

  # W, 16th March 2022. After some troubles, we now do not remove zero pads
  # For very small matrices (<5). This is unnecessary there, and saves us trouble. 
  if ((dim(bitl_f)[1] < 5) | (dim(bitl_f)[2] < 5)){
    return(bitl_f)
  }
  
  # Make sure everything is good. 
  stopifnot("Error: Removing of zero pads only on matrix class." = class(bitl_f) == c('matrix', 'array'))
  stopifnot("Error: Removing of zero pads only on matrices with dim > 1" = all(dim(bitl_f)[1] > 1))
  
  # Attempt to remove first and last cols and rows, if they are 0-only. 
  # stop that process once the matrix does not change anymore, or once
  # it has <5 rows or columns.
  converged = F
  while (converged==F){
    
    bitl_pre = bitl_f
    
    if (all(bitl_f[1,] == 0)){
      bitl_f = bitl_f[2:nrow(bitl_f),]
    }
    if (all(bitl_f[nrow(bitl_f),] == 0)){
      bitl_f = bitl_f[1:nrow(bitl_f)-1,]
    }
    if (all(bitl_f[,1] == 0)){
      bitl_f = bitl_f[,2:ncol(bitl_f)]
    }
    if (all(bitl_f[,ncol(bitl_f)] == 0)){
      bitl_f = bitl_f[,1:ncol(bitl_f)-1]
    }
    
    # If matrix does not change anymore, or if it is becoming very small, then stop.
    if (identical(bitl_pre,bitl_f) | (dim(bitl_f)[1] < 5) | (dim(bitl_f)[2] < 5)){
      converged = T
    }
    
  }

  return(bitl_f)
  
}
