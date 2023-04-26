#' Bresenham's Integer Line Drawing Algorithm
#'
#' Generate integer x,y points between vertices by Bresenham's algorithm.
#'
#' @param x numeric vector of the x coordinates of the points to be joined.
#' @param y numeric vector of the y coordinates of the points to be joined.
#' @param gridpoints_x numeric vector of the x coordinates of the grid points.
#' @param gridpoints_y numeric vector of the y coordinates of the grid points.
#' @param close logical value indicating if the points form a closed polygon.
#' @param debug logical value indicating whether to enter debugging mode or not.
#'
#' @return A list of length 2 with \code{integer} \code{x,y} coordinates connecting
#'   the vertices
#'
#' @export
#'
bresenham <-
  function(x,
           y = NULL,
           gridpoints_x,
           gridpoints_y,
           close = TRUE,
           debug = F)
  {
    if (debug) {
      browser()
    }
    griddiffs_x = diff(gridpoints_x)
    griddiffs_y = diff(gridpoints_y)
    # accept any coordinate structure
    v <- xy.coords(
      x = x,
      y = y,
      recycle = TRUE,
      setLab = FALSE
    )
    if (!all(is.finite(v$x), is.finite(v$y)))
      stop("finite coordinates required")
    
    # These are coordinates.
    # We adjust the ends to cling to the nearest corner.
    v$x = c(gridpoints_x[abs(gridpoints_x - x[1]) == min(abs(gridpoints_x - x[1]))],
            gridpoints_x[abs(gridpoints_x - x[2]) == min(abs(gridpoints_x - x[2]))])
    v$y = c(gridpoints_y[abs(gridpoints_y - y[1]) == min(abs(gridpoints_y - y[1]))],
            gridpoints_y[abs(gridpoints_y - y[2]) == min(abs(gridpoints_y - y[2]))])
    
    
    # process all vertices in pairs
    # Find points for the next x value.
    i = 1
    x <- v$x[i] # coordinates updated in x, y
    y <- v$y[i]
    x.end <- v$x[i + 1]
    y.end <- v$y[i + 1]
    
    dx <- abs(x.end - x)
    dy <- -abs(y.end - y)
    
    err = dx + dy
    x_grididx = which(gridpoints_x == v$x[1])
    y_grididx = which(gridpoints_y == v$y[1])
    
    # collect result in 'ans, staring with 1st point
    
    # v$x <- c(v$x, v$x[1])
    # v$y <- c(v$y, v$y[1])
    
    # Hacky version to print the pixel below, for negative ranges.
    if (y > y.end) {
      ans <-
        data.frame(x = x_grididx,
                   y = y_grididx - 1,
                   z = -(griddiffs_x[x_grididx]))
    } else {
      ans <-
        data.frame(x = x_grididx,
                   y = y_grididx,
                   z = (griddiffs_x[x_grididx]))
    }
    # lapply(v[1:2], "[", 1)
    # process one
    # We do this until we are EXACTLY at x.end and y.end
    #while(!(isTRUE(all.equal(x, x.end)) && isTRUE(all.equal(y, y.end)))) {
    
    # What about stopping as soon as one of the contestants has reached their goal?
    while (!(isTRUE(all.equal(x, x.end)) |
             isTRUE(all.equal(y, y.end)))) {

      # Here is where the omnidirectionality is encoded.
      # s perhaps for 'step'?
      
      # We will have to initiate the stepcount at the right position in the future.
      # Or +1 and 0.
      sx <-
        ifelse(x < x.end, griddiffs_x[x_grididx],-griddiffs_x[x_grididx - 1])
      sy <-
        ifelse(y < y.end, griddiffs_y[y_grididx],-griddiffs_y[y_grididx - 1])

      e2 <- 2 * err
      if (e2 >= dy | e2 <= dx) {
        # increment x
        
        # increment x
        x <- x + sx
        x_grididx <- x_grididx + sign(sx)
        
        
        # recalc error
        dx <- abs(x.end - x)
        dy <- -abs(y.end - y)
        err <- dx + dy
        
        #}
        #if (e2 <= dx) { # increment y
        
        y <- y + sy
        y_grididx <- y_grididx + sign(sy)
        
        dx <- abs(x.end - x)
        dy <- -abs(y.end - y)
        err <- dx + dy
        
      } else {
        print('Algorithm is likely stuck')
        
      }
      if (y > y.end) {
        ans = rbind(ans, c(x_grididx, y_grididx - 1,-(griddiffs_x[x_grididx])))
        
      } else {
        ans = rbind(ans, c(x_grididx, y_grididx, (griddiffs_x[x_grididx])))
      }
      
    }
    
    # remove duplicated points (typically 1st and last)
    dups <- duplicated(do.call(cbind, ans), MARGIN = 1)
    
    to_return_nodup = as.data.frame(lapply(ans, "[",!dups))
    to_return <- to_return_nodup[-nrow(to_return_nodup), ]
    #to_return = to_return_nodup
    return(to_return)
  }
