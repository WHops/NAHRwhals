#Gridlocus improve
library(ggplot2)
devtools::load_all()
paf = read.table('~/PhD/projects/nahrcall/revisions/WGS_mode/bigpaf.paf', header=T)
paf = read.table('res/chr1-144923346-150047635/diff/paf/T2T_NA12_xy.paf', header=T)
paf = read.table('~/Desktop/test.tsv', header=T)
# colnames_paf <- c(
#   "qname",
#   "qlen",
#   "qstart",
#   "qend",
#   "strand",
#   "tname",
#   "tlen",
#   "tstart",
#   "tend",
#   "nmatch",
#   "alen",
#   "mapq"
# )
# colnames(paf)[1:length(paf)] <- colnames_paf
# 
# paf <- transform(
#   paf,
#   qend = ifelse(strand == "-", qstart, qend),
#   qstart = ifelse(strand == "-", qend, qstart)
# )

paf_df = load_and_prep_paf_for_gridplot_transform('res/chr1-144923346-150047635/diff/paf/T2T_NA12_xy.paf', minlen=1000, compression=1000)



gxy = make_xy_grid_fast(paf_df)

xdf = data.frame(x=gxy[[1]])
ydf = data.frame(y=gxy[[2]])

ggplot(paf_df) + geom_segment(aes(x=tstart, xend=tend, y=qstart, yend=qend)) +  geom_hline(data = ydf, aes(yintercept=y)) + geom_vline(data= xdf, aes(xintercept=x))#+  geom_hline(data = yspdf, aes(yintercept=y), color='red') + geom_vline(data= xspdf, aes(xintercept=x), color='red')


microbenchmark(grid_list = fill_matrix_bressiexchange(paf_df, gxy),times = 2)
grid_list = fill_matrix_bressiexchange(paf_df, gxy)

grid_list = data.frame()
griddiffs_x = diff(gxy[[1]])
griddiffs_y = diff(gxy[[2]])
#microbenchmark(

prof = profvis(bressiwrap())
microbenchmark(gl = bressiwrap(), times = 2)
bressiwrap <- function() {
  
  # Create an empty list to store dataframes
  df_list <- vector("list", dim(paf_df)[1])
  
  for (i in 1:dim(paf_df)[1]) {
    df_list[[i]] <- as.data.frame(
      bresenham(
        x = as.numeric(paf_df[i, c("tstart", "tend")]),
        y = as.numeric(paf_df[i, c("qstart", "qend")]),
        gxy[[1]],
        gxy[[2]],
        griddiffs_x,
        griddiffs_y,
        debug = F
      )
    )
  }
  
  # Bind all dataframes together
  grid_list <- do.call(rbind, df_list)
  
  return(grid_list)
}
#)
gl = bressiwrap()
plot_matrix_ggplot(gl)

#' @export
fill_matrix_bressiexchange <- function(paf_df, gridlines_xy){
  
  #
  gridlines_x = gridlines_xy[[1]]
  gridlines_y = gridlines_xy[[2]]
  
  # Compute midpoints for all quadratic intervals
  mid_x <- (head(gridlines_x, -1) + tail(gridlines_x, -1)) / 2
  mid_y <- (head(gridlines_y, -1) + tail(gridlines_y, -1)) / 2
  
  midpoints = expand.grid(x = mid_x, y = mid_y)
  
  traversal_count <- data.frame(x = midpoints[,1],
                                y = midpoints[,2])
  
  return_grid = expand.grid(x=(1:(length(gridlines_x)-1)),
                            y=(1:(length(gridlines_y)-1)), 
                            z=0)
  return_grid$len = rep(diff(gridlines_x), length(gridlines_y)-1)
  
  # Check overlap for each paf entry
  
  paf_df_slopes = paf_df$slope
  paf_df_tstarts = paf_df$tstart
  paf_df_tends = paf_df$tend
  paf_df_qstarts = paf_df$qstart
  paf_df_qends = paf_df$qend

  midpoints_x = midpoints$x
  midpoints_y = midpoints$y
  
  for (row in 1:nrow(paf_df)) {
      print(row)
      c <- paf_df_qstarts[row] - (paf_df_slopes[row] * paf_df_tstarts[row])
      overlapped <- which(
          midpoints_x > paf_df_tstarts[row] & 
          midpoints_x < paf_df_tends[row] & 
          midpoints_y == ((midpoints_x * paf_df_slopes[row]) + c)
      )
      return_grid[overlapped, "z"] <- return_grid[overlapped, 'len'] * paf_df_slopes[row]
    
  }
  
  return_grid$len = NULL
  return_grid = return_grid[return_grid$z != 0,]
  
  return(return_grid)
}



# 
# points <- c(3, 6, 10)
# intervals <- matrix(c(1,2, 3,7, 2,4, 7,8), ncol=2)
# 
# points <- c(3, 6, 10)
# intervals <- matrix(c(1,2, 3,7, 2,4, 7,8), ncol=2)
# 
# # Sort points and intervals
# points <- sort(points)
# order_indices <- order(intervals[,1])
# sorted_intervals <- intervals[order_indices, ]
# 
# find_overlapping_intervals <- function(points, sorted_intervals, original_order) {
#   result <- integer(0)
#   p_ptr <- 1
#   i_ptr <- 1
#   
#   while (p_ptr <= length(points) && i_ptr <= nrow(sorted_intervals)) {
#     if (points[p_ptr] >= sorted_intervals[i_ptr, 1] && points[p_ptr] <= sorted_intervals[i_ptr, 2]) {
#       result <- c(result, which(original_order == i_ptr)) # Get the original index
#       p_ptr <- p_ptr + 1
#     } else if (points[p_ptr] < sorted_intervals[i_ptr, 1]) {
#       p_ptr <- p_ptr + 1
#     } else {
#       i_ptr <- i_ptr + 1
#     }
#   }
#   
#   return(result)
# }
# 
# overlapping_indices <- find_overlapping_intervals(points, sorted_intervals, order_indices)
# print(overlapping_indices)
# 
# intervals = matrix(c(paf$tstart, paf$tend), ncol = 2)
# order_indices <- order(intervals[,1])
# sorted_intervals <- intervals[order_indices, ]
# 
# 
# t0 = unique(paf$tstart)
# q0 = unique(paf$qstart)
# 
# i = 30
# t1 = find_overlapping_intervals(t0[i], sorted_intervals, order_indices)
# print(t1)
# 
# 
# 




















