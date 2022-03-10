#" matrix_remove_zero_pads
#' Quick and dirty way of getting rid of padding zero columns or rows. 
#' @author Wolfram Häps
#' @export
matrix_remove_zero_pads <- function(bitl_f, iterations = 100){
  
  # Make sure everything is good. 
  stopifnot("Error: Removing of zero pads only on matrix class." = class(bitl_f) == c('matrix', 'array'))
  stopifnot("Error: Removing of zero pads only on matrices with dim > 1" = all(dim(bitl_f)[1] > 1))
  
  # Attempt to remove first and last cols and rows, if they are 0-only. 
  # stop that process once the matrix does not change anymore. 
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
    
    if (identical(bitl_pre,bitl_f)){
      converged = T
    }
    
  }
  return(bitl_f)
  
}