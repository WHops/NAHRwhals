#" matrix_remove_zero_pads
#' Quick and dirty way of getting rid of padding zero columns or rows. 
#' @author Wolfram Häps
#' @export
matrix_remove_zero_pads <- function(bitl_f){
  
  # If all outside cols and rows have entries, we dont have to cut anything...
  if ((any(bitl_f[1,] != 0)) & 
      (any(bitl_f[nrow(bitl_f),] != 0)) & 
      (any(bitl_f[,1] != 0)) &
      (any(bitl_f[,ncol(bitl_f)] != 0)) 
  ){
    return(bitl_f)
  }
  
  bitl_f = as.matrix(bitl_f)

  empty_rows_bool = as.numeric(apply(bitl_f==0, 1, all))
  empty_cols_bool = as.numeric(apply(bitl_f==0, 2, all))
  
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
