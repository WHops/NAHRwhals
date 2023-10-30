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
           griddiffs_x,
           griddiffs_y,
           close = TRUE,
           debug = F) {
    if (debug) {
      browser()
    }

    # accept any coordinate structure
    v <- xy.coords(
      x = x,
      y = y,
      recycle = TRUE,
      setLab = FALSE
    )
    if (!all(is.finite(v$x), is.finite(v$y))) {
      stop("finite coordinates required")
    }

    # process all vertices in pairs
    # Find points for the next x value.
    i <- 1
    x <- v$x[i] # coordinates updated in x, y
    y <- v$y[i]
    x.end <- v$x[i + 1]
    y.end <- v$y[i + 1]

    dx <- abs(x.end - x)
    dy <- -abs(y.end - y)

    err <- dx + dy
    x_grididx <- which(gridpoints_x == v$x[1])
    y_grididx <- which(gridpoints_y == v$y[1])

    # Hacky version to print the pixel below, for negative ranges.
    # if (y > y.end) {
    #   ans <-
    #     data.frame(
    #       x = x_grididx,
    #       y = y_grididx - 1,
    #       z = -(griddiffs_x[x_grididx])
    #     )
    # } else {
    #   ans <-
    #     data.frame(
    #       x = x_grididx,
    #       y = y_grididx,
    #       z = (griddiffs_x[x_grididx])
    #     )
    # }

    ans_list <- list()
    counter <- 1
    
    # What about stopping as soon as one of the contestants has reached their goal?
    while ((x != x.end) & (y != y.end)) {
      # We will have to initiate the stepcount at the right position in the future.
      # Or +1 and 0.
      sx <-
        ifelse(x < x.end, griddiffs_x[x_grididx], -griddiffs_x[x_grididx - 1])
      sy <-
        ifelse(y < y.end, griddiffs_y[y_grididx], -griddiffs_y[y_grididx - 1])

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

        # }
        # if (e2 <= dx) { # increment y

        y <- y + sy
        y_grididx <- y_grididx + sign(sy)

        dx <- abs(x.end - x)
        dy <- -abs(y.end - y)
        err <- dx + dy
      } else {
        print("Algorithm is likely stuck")
      }
      if (y > y.end) {
        ans_list[[counter]] <- c(x_grididx, y_grididx - 1, -(griddiffs_x[x_grididx]))
      } else {
        ans_list[[counter]] <- c(x_grididx, y_grididx, (griddiffs_x[x_grididx]))
      }
      
      counter <- counter + 1
    }

    ans <- unique(do.call(rbind, ans_list))
    colnames(ans) = c('x','y','z')
    return(ans)

  }
