#Gridlocus improve
library(ggplot2)
library(microbenchmark)
library(profvis)
options(scipen=99999)
mb = microbenchmark
devtools::load_all()
paf = read.table('~/PhD/projects/nahrcall/revisions/WGS_mode/bigpaf.paf', header=T)
#paf = read.table('res/chr1-144923346-150047635/diff/paf/T2T_NA12_xy.paf', header=T)
#paf = read.table('~/Desktop/test.tsv', header=T)

paf = read.table('~/Desktop/paffast.paf', header=T)

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
paf <- transform(
  paf,
  qend = ifelse(strand == "-", qstart, qend),
  qstart = ifelse(strand == "-", qend, qstart)
)

#paf_df = load_and_prep_paf_for_gridplot_transform('res/chr1-144923346-150047635/diff/paf/T2T_NA12_xy.paf', minlen=10000, compression=10000)

#paf_df = read.table('~/Desktop/daemon.tsv', header=T)

#paf_df = paf
gxy = make_xy_grid_fast(paf_df)
microbenchmark(make_xy_grid_fast(paf_df), times=2)
xdf = data.frame(x=gxy[[1]])
ydf = data.frame(y=gxy[[2]])

plot_x_y = ggplot(paf_df) + geom_segment(aes(x=tstart, xend=tend, y=qstart, yend=qend))# + xlim(c(37000000, 37500000))# +  geom_hline(data = ydf, aes(yintercept=y)) + geom_vline(data= xdf, aes(xintercept=x))#+  geom_hline(data = yspdf, aes(yintercept=y), color='red') + geom_vline(data= xspdf, aes(xintercept=x), color='red')
plot_x_y

#microbenchmark(grid_list = fill_matrix_bressiexchange(paf_df, gxy),times = 2)
#grid_list = fill_matrix_bressiexchange(paf_df, gxy)

#grid_list = data.frame()

#microbenchmark(

#prof = profvis({bressiwrap(paf_df, gxy)}, interval=0.001)
microbenchmark(gl = bressiwrap(paf_df, gxy), times = 2)
#res = solve_mutation(mat, maxdepth = maxdepth, solve_th = 101, compression = 50000, is_cluttered_already_paf = F) # , discovery_exact = params$discovery_exact)

#)
gl = bressiwrap(paf_df, gxy)
plot_matrix_ggplot(as.data.frame(gl))
mat = gridlist_to_gridmatrix(list(gxy[[1]], gxy[[2]], gl))


microbenchmark(make_segmented_pairwise_plot_home(gxy, plot_x_y), times=2)


make_segmented_pairwise_plot <- function(grid_xy, plot_x_y) {

  x_vals = grid_xy[[1]]
  length_x_vals = length(x_vals)
  
  y_vals = grid_xy[[2]]
  length_y_vals = length(y_vals)
  
  xstart <- x_vals[1:length_x_vals - 1]
  xend <-   x_vals[2:length_x_vals]
  
  ystart <- y_vals[1:length_y_vals - 1]
  yend <-  y_vals[2:length_y_vals]
  
  xmax <- max(x_vals)
  ymax <- max(y_vals)
  
  datx <- data.frame(
    xstart = xstart,
    xend = xend,
    xmax = xmax
  )
  daty <- data.frame(
    yend = yend,
    ymax = ymax,
    ystart = ystart
  )
  
  likely_stepsize <- min(c(diff(x_vals), diff(y_vals)))
  
  # Introducing special case for nrow(datx) = 1
  if (likely_stepsize == Inf) {
    likely_stepsize <- 0
  }
  
  if (length(xstart) > 433) {
    print("Too many segments to make colored plot")
    return()
  }
  
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == "qual", ]
  col_vector <- rep(unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))), 10)
  
  max_x <- max(ggplot2::layer_scales(plot_x_y)$x$range$range, max(xmax) + likely_stepsize)
  max_y <- max(ggplot2::layer_scales(plot_x_y)$y$range$range, max(ymax) + likely_stepsize)
  
  min_x <- min(ggplot2::layer_scales(plot_x_y)$x$range$range, min(xstart) - likely_stepsize)
  min_y <- min(ggplot2::layer_scales(plot_x_y)$y$range$range, min(ystart) - likely_stepsize)
  
  plot_x_y_segmented <- plot_x_y +
    ggplot2::geom_rect(
      data = datx,
      ggplot2::aes(xmin = xstart, xmax = xend, ymin = 0, ymax = ymax, fill = col_vector[1:length(xstart)]),
      alpha = 0.5
    ) +
    ggplot2::guides(fill = FALSE) +
    ggplot2::coord_cartesian(xlim = c(min_x, max_x), ylim = c(min_y, max_y)) +
    ggplot2::scale_x_continuous(labels = scales::comma) +
    ggplot2::scale_y_continuous(labels = scales::comma)
  # ggplot2::geom_segment(data=daty,
  #             ggplot2::aes(x=0, xend=xmax, y=ystart, yend=ystart), color='grey')
  print(plot_x_y_segmented)
  

}
