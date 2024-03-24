#' gridlist_to_gridmatrix
#'
#' Convert gridlist to grimatrix
#' @param grid a dataframe, columns (x,y,z).
#' @return A matrix of the bitplot
#' @author Wolfram Hoeps
#' @export
gridlist_to_gridmatrix <- function(grid) {
  # Sort by x
  grid_list <- grid[[3]]
  dim_x <- length(grid[[1]]) - 1
  dim_y <- length(grid[[2]]) - 1

  grid_list <- grid_list[order(grid_list$x), ]

  # Remove duplicates if there are any (there shouldn't be...)
  grid_list <- remove_duplicates_triple(grid_list)

  x_missing <- which(1:dim_x %in% grid_list$x == F)
  y_missing <- which(1:dim_y %in% grid_list$y == F)

  # W, 7th March 2022.
  # Fixed a bug here, where previously points might be appended
  # to the grid_list that were inappropriate if the grid was not
  # symmetrical.
  for (xm in x_missing) {
    grid_list <- rbind(grid_list, c(xm, 1, 0))
  }
  for (ym in y_missing) {
    grid_list <- rbind(grid_list, c(1, ym, 0))
  }

  # Remove duplicates AGAIN.
  # afaik, the only possible duplicates is if '1,1' was appended
  # twice.
  grid_list <- remove_duplicates_triple(grid_list)

  # list to matrix
  suppressMessages(suppressWarnings({
    gridmatrix <- reshape2::dcast(grid_list, y ~ x, fill = 0)
  }))
  gridmatrix$x <- NULL
  gridmatrix$y <- NULL
  gridmatrix <- as.matrix(gridmatrix)

  # Columns and rownames show the size of the grid.
  colnames(gridmatrix) <- diff(grid[[1]])
  row.names(gridmatrix) <- diff(grid[[2]])

  return(gridmatrix)
}
