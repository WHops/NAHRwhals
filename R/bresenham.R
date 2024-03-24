#' @export
bressiwrap <- function(paf_df, gxy) {
  
  griddiffs_x = diff(gxy[[1]])
  griddiffs_y = diff(gxy[[2]])
  # Create an empty list to store dataframes
  df_list <- vector("list", 1)#length(griddiffs_x)*length(griddiffs_y)*0.1)
  x_all = list()
  y_all = list()
  z_all = list()
  
  x_idxs = match(paf_df$tstart, gxy[[1]])
  y_idxs = match(paf_df$qstart, gxy[[2]])

  pafdf_tstarts = as.numeric(paf_df$tstart)
  pafdf_tends = as.numeric(paf_df$tend)
  pafdf_qstarts = as.numeric(paf_df$qstart)
  pafdf_qends = as.numeric(paf_df$qend)

  for (i in 1:dim(paf_df)[1]) {
    xyz <- 
      bresenham(
        x_in = c(pafdf_tstarts[i], pafdf_tends[i]),
        y_in = c(pafdf_qstarts[i], pafdf_qends[i]),
        gxy[[1]],
        gxy[[2]],
        griddiffs_x,
        griddiffs_y,
        x_idxs[i],
        y_idxs[i]
      )
    
    x_all[[i]] = xyz[[1]]
    y_all[[i]] = xyz[[2]]
    z_all[[i]] = xyz[[3]]
  }
  
  # Combine all x, y, and z values together
  x_combined = unlist(x_all)
  y_combined = unlist(y_all)
  z_combined = unlist(z_all)
  
  # Convert to a data frame
  grid_list <- unique(data.frame(x = x_combined, y = y_combined, z = z_combined))
  
  return(grid_list)
}


#' Bresenham's Integer Line Drawing Algorithm
#'
#' Generate integer x,y points between vertices by Bresenham's algorithm.
#'
#' @param x_in numeric vector of the x coordinates of the points to be joined.
#' @param y_in numeric vector of the y coordinates of the points to be joined.
#' @param gridpoints_x numeric vector of the x coordinates of the grid points.
#' @param gridpoints_y numeric vector of the y coordinates of the grid points.
#' @param griddiffs_x numeric vector of the differences between the x coordinates
#' @param griddiffs_y numeric vector of the differences between the y coordinates
#'
#' @return A list of length 2 with \code{integer} \code{x,y} coordinates connecting
#'   the vertices
#'
#' @export
#'
bresenham <-
  function(x_in,
           y_in,
           gridpoints_x,
           gridpoints_y,
           griddiffs_x,
           griddiffs_y,
           x_grididx,
           y_grididx) {


    # process all vertices in pairs
    # Find points for the next x value.
    x <- x_in[1] # coordinates updated in x, y
    y <- y_in[1]
    x.end <- x_in[2]
    y.end <- y_in[2]

    ans_list <- list()
    x_all = numeric()
    y_all = numeric()
    z_all = numeric()
    counter <- 1
    
    pos_sign = (y < y.end)
    # What about stopping as soon as one of the contestants has reached their goal?
    while ((x != x.end) & (y != y.end)){
      # We will have to initiate the stepcount at the right position in the future.
      # Or +1 and 0.
      
      sx <- griddiffs_x[x_grididx]
      
      if (pos_sign){
        sy <- griddiffs_y[y_grididx]
      } else {
        sy <- -griddiffs_y[y_grididx - 1]
      }

      # increment x,y
      x <- x + sx
      y <- y + sy

      # Not down that we have incremented
      x_all[counter] = x_grididx
      
      if (!pos_sign) {
        y_all[counter] <- y_grididx - 1
        z_all[counter] <- -sx
        
      } else if (pos_sign) {
        y_all[counter] <- y_grididx
        z_all[counter] <- sx
      }

      x_grididx <- x_grididx + 1
      y_grididx <- y_grididx + sign(sy)

      counter <- counter + 1
    }

    
    return(list(x_all, y_all, z_all))

  }
