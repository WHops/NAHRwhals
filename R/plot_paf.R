#' plot_paf
#' Make a plot of a paf data frame. Should contain tstart, tend, qstart, qend. 
#' Optionally, gridlines can be plotted, too. 
#' @author Wolfram Hoeps
#' @export
plot_paf <- function(paf_df,
                     gridlines.x = NULL,
                     gridlines.y = NULL,
                     wiggle = 0,
                     plot_gridlines = T, 
                     linewidth = 1) {
  plot = ggplot2::ggplot() +
    ggplot2::geom_segment(data = paf_df,
                          ggplot2::aes(
                            x = tstart,
                            xend = tend,
                            y = qstart,
                            yend = qend
                          )) +
    ggplot2::coord_fixed(
      ratio = 1,
      xlim = NULL,
      ylim = NULL,
      expand = TRUE,
      clip = "on"
    ) +
    ggplot2::theme_bw() + 
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())
  
  if (plot_gridlines) {
    plot = plot +
      ggplot2::geom_hline(ggplot2::aes(yintercept = gridlines.y + runif(
        length(gridlines.y),-wiggle, wiggle
      ))
      , color =
        'grey', size = linewidth) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = gridlines.x + runif(
        length(gridlines.x), wiggle, wiggle
      ))
      , color =
        'grey', size = linewidth) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())
  }
  return(plot)
}
