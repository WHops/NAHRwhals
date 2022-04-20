#" matrix_remove_zero_pads
#' Quick and dirty way of getting rid of padding zero columns or rows. 
#' @author Wolfram Häps
#' @export
matrix_remove_zero_pads <- function(bitl_f){
  bitl_f = as.matrix(bitl_f)

  empty_rows_bool = as.numeric(apply(as.matrix(bitl_f)==0, 1, all))
  empty_cols_bool = as.numeric(apply(as.matrix(bitl_f)==0, 2, all))
  
  #Filter down to only intervals touching the border. 
  leading_rows_to_rm = which(cumprod(empty_rows_bool) > 0)
  trailing_rows_to_rm = (length(empty_rows_bool) + 1) - which(cumprod(rev(empty_rows_bool)) > 0)
  
  leading_cols_to_rm = which(cumprod(empty_cols_bool) > 0)
  trailing_cols_to_rm = (length(empty_cols_bool) + 1) - which(cumprod(rev(empty_cols_bool)) > 0)
  
  rows_to_keep = (1:dim(bitl_f)[1])[!((1:dim(bitl_f)[1]) %in% c(leading_rows_to_rm, trailing_rows_to_rm))]
  cols_to_keep = (1:dim(bitl_f)[2])[!((1:dim(bitl_f)[2]) %in% c(leading_cols_to_rm, trailing_cols_to_rm))]
  
  # Filter bitl
  bitl_f_filter = as.matrix(bitl_f[rows_to_keep, cols_to_keep])
  
  if (all(dim(bitl_f_filter) == c(0,0))){
    return(NULL)
  }
    
  return(bitl_f_filter)
}
# 
#   if ((dim(bitl_f)[1] == 1) & (dim(bitl_f)[2] > 1)){
#     
#   }
#   # W, 16th March 2022. After some troubles, we now do not remove zero pads
#   # For very small matrices (<5). This is unnecessary there, and saves us trouble. 
#   if ((dim(bitl_f)[1] < 5) | (dim(bitl_f)[2] < 5)){
#     return(bitl_f)
#   }
#   
#   # Make sure everything is good. 
#   stopifnot("Error: Removing of zero pads only on matrix class." = class(bitl_f) == c('matrix', 'array'))
#   stopifnot("Error: Removing of zero pads only on matrices with dim > 1" = all(dim(bitl_f)[1] > 1))
#   
#   # Attempt to remove first and last cols and rows, if they are 0-only. 
#   # stop that process once the matrix does not change anymore, or once
#   # it has <5 rows or columns.
#   converged = F
#   while (converged==F){
#     
#     bitl_pre = bitl_f
#     
#     if (all(bitl_f[1,] == 0)){
#       bitl_f = bitl_f[2:nrow(bitl_f),]
#     }
#     if (all(bitl_f[nrow(bitl_f),] == 0)){
#       bitl_f = bitl_f[1:nrow(bitl_f)-1,]
#     }
#     if (all(bitl_f[,1] == 0)){
#       bitl_f = bitl_f[,2:ncol(bitl_f)]
#     }
#     if (all(bitl_f[,ncol(bitl_f)] == 0)){
#       bitl_f = bitl_f[,1:ncol(bitl_f)-1]
#     }
#     
#     # If matrix does not change anymore, or if it is becoming very small, then stop.
#     if (identical(bitl_pre,bitl_f) | (dim(bitl_f)[1] < 5) | (dim(bitl_f)[2] < 5)){
#       converged = T
#     }
#     
#   }
# 
#   return(bitl_f)
#   
# }
