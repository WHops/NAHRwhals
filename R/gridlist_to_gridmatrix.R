#' gridlist_to_gridmatrix
#' 
#' Convert gridlist to grimatrix
#' @param grid_list a dataframe, columns (x,y,z).
#' @return A matrix of the bitplot
#' @author Wolfram Höps
#' @export
gridlist_to_gridmatrix <- function(grid_list){
  
  # Sort by x
  
  grid_list = grid_list[order(grid_list$x),]
  
  # EXPERIMENTAL PART #
  # Remove duplicates TRIPLE. First, remove obvious duplicates (all 3 same)
  # Then, identify those with same x/y but different z. The get 0.
  # Finally, the second step has produced new duplicates, so be mean again. 
  
  # Remove duplicates. This happens if SDs overlap (which can be a result
  # of the compression step)
  grid_list = grid_list[!duplicated(grid_list),]
  
  # If a value appears both positive and negative, we have to set it to 0, 
  # rather than keeping randomly one of the two.
  grid_list[duplicated(grid_list[,c('x','y')]),'z'] = 0
  
  # Remove duplicates again. 
  grid_list = grid_list[!duplicated(grid_list[,c('x','y')], fromLast=T),]


  
  x_missing = which(1:max(grid_list$x) %in% grid_list$x == F)
  y_missing = which(1:max(grid_list$y) %in% grid_list$y == F)
  
  for (xm in x_missing) {
    grid_list = rbind(grid_list, c(xm, xm, 0))
  }
  for (ym in y_missing) {
    grid_list = rbind(grid_list, c(ym, ym, 0))
  }
  
  
  gridmatrix = reshape2::dcast(grid_list, y ~ x, fill = 0)
  gridmatrix$x = NULL
  gridmatrix$y = NULL
  gridmatrix = as.matrix(gridmatrix)
  
  colnames(gridmatrix) = 1:dim(gridmatrix)[2]
  rownames(gridmatrix) = 1:dim(gridmatrix)[1]
  
  return(gridmatrix)
}
