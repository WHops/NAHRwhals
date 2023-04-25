
#' Remove duplicate values from a triple data frame
#'
#' This function takes a triple data frame containing columns x, y, and z, and removes duplicate values. First, obvious duplicates (where all three values are the same) are removed. Then, values with the same x/y but different z are identified and set to 0. Finally, the second step may produce new duplicates, so the function removes those again. This function is typically used to simplify a grid_list.
#'
#' @param grid_list_f A data frame containing columns x, y, and z
#'
#' @return A data frame with duplicate values removed
#'
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