# bresenham_orig <- function(x, y = NULL, close = TRUE)
# {
#   # accept any coordinate structure
#   v <- xy.coords(
#     x = x,
#     y = y,
#     recycle = TRUE,
#     setLab = FALSE
#   )
#   if (!all(is.finite(v$x), is.finite(v$y)))
#     stop("finite coordinates required")
#   
#   v[1:2] <-
#     lapply(v[1:2], round) # Bresenham's algorithm IS for integers
#   nx <- length(v$x)
#   if (nx == 1)
#     return(list(x = v$x, y = v$y)) # just one point
#   if (nx > 2 &&
#       close == TRUE) {
#     # close polygon by replicating 1st point
#     v$x <- c(v$x, v$x[1])
#     v$y <- c(v$y, v$y[1])
#     nx <- nx + 1
#   }
#   # collect result in 'ans, staring with 1st point
#   ans <- lapply(v[1:2], "[", 1)
#   
#   # process all vertices in pairs
#   for (i in seq.int(nx - 1)) {
#     x <- v$x[i] # coordinates updated in x, y
#     y <- v$y[i]
#     x.end <- v$x[i + 1]
#     y.end <- v$y[i + 1]
#     
#     dx <- abs(x.end - x)
#     dy <- -abs(y.end - y)
#     sx <- ifelse(x < x.end, 1,-1)
#     sy <- ifelse(y < y.end, 1,-1)
#     err <- dx + dy
#     
#     # process one segment
#     while (!(isTRUE(all.equal(x, x.end)) &&
#              isTRUE(all.equal(y, y.end)))) {
#       e2 <- 2 * err
#       if (e2 >= dy) {
#         # increment x
#         err <- err + dy
#         x <- x + sx
#       }
#       if (e2 <= dx) {
#         # increment y
#         err <- err + dx
#         y <- y + sy
#       }
#       ans$x <- c(ans$x, x)
#       ans$y <- c(ans$y, y)
#     }
#   }
#   # remove duplicated points (typically 1st and last)
#   dups <- duplicated(do.call(cbind, ans), MARGIN = 1)
#   return(lapply(ans, "[",!dups))
# }



#' Bresenham's Integer Line Drawing Algorithm
#'
#' Generate integer x,y points between vertices by Bresenham's algorithm.
#'
#' @param x,y the x and y coordinates a points to be joined where y can
#'  can be missing and the argument x is processed by \code{\link{xy.coords}}.
#' @param close \code{logical} value indicating if the points form a closed
#'  polygon (without duplicating the first and last points)
#'
#' @return
#'
#' A list of length 2 with \code{integer} \code{x,y} coordinates connecting
#'   the vertices
#' @examples
#' # simple line
#'   bresenham(x = c(1, 4), y = c(1, 12))
#' # closed and open polygon
#'   verts <- list(x = c(1, 9, 6, 4, 2), y = c(9, 9, 2, 3, 1))
#'   plot(verts, type = "l", ylim = c(0, 10), xlim = c(0, 10))
#'   points(bresenham(verts)) # closed
#'   points(bresenham(verts, close = FALSE), col = 2, pch = 16)
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
                   z = -log10(griddiffs_x[x_grididx]))
    } else {
      ans <-
        data.frame(x = x_grididx,
                   y = y_grididx,
                   z = log10(griddiffs_x[x_grididx]))
    }
    # lapply(v[1:2], "[", 1)
    # process one
    # We do this until we are EXACTLY at x.end and y.end
    #while(!(isTRUE(all.equal(x, x.end)) && isTRUE(all.equal(y, y.end)))) {
    
    # What about stopping as soon as one of the contestants has reached their goal?
    while (!(isTRUE(all.equal(x, x.end)) |
             isTRUE(all.equal(y, y.end)))) {
      # print(x)
      # print(x.end)
      # Here is where the omnidirectionality is encoded.
      # s perhaps for 'step'?
      
      # We will have to initiate the stepcount at the right position in the future.
      # Or +1 and 0.
      sx <-
        ifelse(x < x.end, griddiffs_x[x_grididx],-griddiffs_x[x_grididx - 1])
      sy <-
        ifelse(y < y.end, griddiffs_y[y_grididx],-griddiffs_y[y_grididx - 1])
      # print(sx)
      # print(dx)
      # print(sy)
      # print(dx)
      # print(dy)
      # browser()
      e2 <- 2 * err
      #print(paste0('e2 >= dy: ', e2 >= dy))
      if (e2 >= dy | e2 <= dx) {
        # increment x
        
        # print('incrementing x')
        # increment x
        x <- x + sx
        x_grididx <- x_grididx + sign(sx)
        
        
        # recalc error
        dx <- abs(x.end - x)
        dy <- -abs(y.end - y)
        err <- dx + dy
        
        #}
        #if (e2 <= dx) { # increment y
        
        # print('incrementing y')
        y <- y + sy
        y_grididx <- y_grididx + sign(sy)
        
        dx <- abs(x.end - x)
        dy <- -abs(y.end - y)
        err <- dx + dy
        
      } else {
        print('Algorithm is likely stuck')
        
      }
      if (y > y.end) {
        ans = rbind(ans, c(x_grididx, y_grididx - 1,-log10(griddiffs_x[x_grididx])))
        
      } else {
        ans = rbind(ans, c(x_grididx, y_grididx, log10(griddiffs_x[x_grididx])))
      }
      
    }
    
    # remove duplicated points (typically 1st and last)
    dups <- duplicated(do.call(cbind, ans), MARGIN = 1)
    
    to_return_nodup = as.data.frame(lapply(ans, "[",!dups))
    to_return <- to_return_nodup[-nrow(to_return_nodup), ]
    #to_return = to_return_nodup
    return(to_return)
  }
#
# gridlines_x = c(0,999,2000,3999,5000,6003,6997,8002,8996,10000)
# gridlines_y = gridlines_x
#
#
# # library(ggplot2)



#
# 159921,
# 40802 ,
# 39989 ,
#
gridlines_x = c( 40004 , 40855 , 82718, 116112, 160000,160804  )
gridlines_y = gridlines_x #c( 0,  21714  ,39989  ,40802,  82718, 116112 ,159921 ,160804 ,179018)
x = c(116112, 160804)
y = c(82718, 40004)

plot = ggplot2::ggplot() +
  ggplot2::geom_hline(ggplot2::aes(yintercept=gridlines_y), color='grey') +
  ggplot2::geom_vline(ggplot2::aes(xintercept=gridlines_x), color='grey') +
  ggplot2::coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  ggplot2::theme_bw() +
  geom_segment(aes(x=x[1], xend=x[2], y=y[1], yend=y[2])) +
  xlim(c(110000,170000)) + ylim(c(35000, 85000))
plot
df = as.data.frame(bresenham(x = x, y = y, gridlines_x, gridlines_y, debug=F))
#
ggplot(df) + geom_tile(aes(x=x, y=y)) +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  theme_bw()
#
#
#
#
#
#
#
#
#
#
#
# gridlines_x = c( 0,  21714 , 40004 , 40855 , 82718 ,116112 ,160000 ,160819, 179018 ,193610)
# gridlines_y = c( 0,  21714  ,39989  ,40802,  82718, 116112 ,159921 ,160804 ,179018 ,193610)
# x = c(116112, 160804)
# y = c(82718, 40004)
#
# plot = ggplot2::ggplot() +
#   ggplot2::geom_hline(ggplot2::aes(yintercept=gridlines_y), color='grey') +
#   ggplot2::geom_vline(ggplot2::aes(xintercept=gridlines_x), color='grey') +
#   ggplot2::coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
#   ggplot2::theme_bw() +
#   geom_segment(aes(x=x[1], xend=x[2], y=y[1], yend=y[2])) +
#   xlim(c(110000,170000)) + ylim(c(35000, 85000))
# plot
# df = as.data.frame(bresenham(x = x, y = y, gridlines_x, gridlines_y, debug=F))
# #
# ggplot(df) + geom_tile(aes(x=x, y=y)) +
#   coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
#   theme_bw()
#
#
#
#
#
#
#
#
#
# gridlines_x = c(570 , 579 , 603 , 697,  802 , 896 , 909 , 930 ,1000)
# gridlines_y =  gridlines_x
#
# x = c(579, 600)
# y = c(930, 909)
#
# plot = ggplot2::ggplot() +
#   ggplot2::geom_hline(ggplot2::aes(yintercept=gridlines_y), color='grey') +
#   ggplot2::geom_vline(ggplot2::aes(xintercept=gridlines_x), color='grey') +
#   ggplot2::coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
#   ggplot2::theme_bw() +
#   geom_segment(aes(x=x[1], xend=x[2], y=y[1], yend=y[2]))
# plot
# df = as.data.frame(bresenham(x = x, y = y, gridlines_x, gridlines_y, debug=F))
# #
# ggplot(df) + geom_tile(aes(x=x, y=y)) +
#   coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
#   theme_bw()
# #
# # library(ggplot2)
# gridpoints_x = c(0,2,4)
# gridpoints_y = c(0,2,4)
# df = as.data.frame(bresenham(x = c(0,2), y = c(0,2), gridpoints_x, gridpoints_y, debug=F))
# #
# ggplot(df) + geom_tile(aes(x=x, y=y)) +
#   #coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
#   theme_bw()
#
# gridpoints_x = c(0,1000,1500,2000,2500, 3100)
# gridpoints_y = c(0,1000,1500,2000,2500, 3100)
# df1 = as.data.frame(bresenham(x = c(0,1500), y = c(0,1500), gridpoints_x, gridpoints_y, debug=F))
# df2 = as.data.frame(bresenham(x = c(2000,3000), y = c(2000,3000), gridpoints_x, gridpoints_y))
# df3 = as.data.frame(bresenham(x = c(1000,2500), y = c(2500,1000), gridpoints_x, gridpoints_y))
# df = rbind(df1, rbind(df2,df3))
#
# ggplot(df) + geom_tile(aes(x=x, y=y, fill=z)) +
#   coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
#   theme_bw()
#
#
#
df1 = as.data.frame(bresenham_orig(x = c(0,1500), y = c(1500,0)))
df2 = as.data.frame(bresenham_orig(x = c(2000,3000), y = c(2000,3000)))
df3 = as.data.frame(bresenham_orig(x = c(1000,1500), y = c(2500,2000)))
df = rbind(df1, rbind(df2,df3))
ggplot(df1) + geom_tile(aes(x=x, y=y)) +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  theme_bw()

# Found example:
