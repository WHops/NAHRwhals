
#' Remove duplicates TRIPLE. First, remove obvious duplicates (all 3 same)
#' Then, identify those with same x/y but different z. The get 0.
#' Finally, the second step has produced new duplicates, so be mean again.
#' @author Wolfram Hoeps
#' @export
remove_duplicates_triple <- function(grid_list_f){
  
  stopifnot("Error in remove_duplicates_triple: Input dataframe does not contain
            columns x,y,z" = colnames(grid_list_f) == c('x','y','z'))
  
  # Remove duplicates. This happens if SDs overlap (which can be a result
  # of the compression step)
  grid_list_f = grid_list_f[!duplicated(grid_list_f),]
  
  # If a value appears both positive and negative, we have to set it to 0,
  # rather than keeping randomly one of the two.
  grid_list_f[duplicated(grid_list_f[, c('x', 'y')]), 'z'] = 0
  
  # Remove duplicates again.
  grid_list_f = grid_list_f[!duplicated(grid_list_f[, c('x', 'y')], fromLast =
                                      T),]
  return(grid_list_f)
  
}