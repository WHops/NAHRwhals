#' Plot_matrix_ggplot
#'
#' Make a plot of a dataframe matrix.
#' @param data_frame_xyz a dataframe with three coordinate columns (x,y,z)
#' @author Wolfram Hoeps
#' @export
plot_matrix_ggplot <- function(data_frame_xyz) {
  p = ggplot2::ggplot(data_frame_xyz) + ggplot2::geom_tile(ggplot2::aes(
    x = x,
    y = y,
    fill = sign(z) * log10(abs(z))
  )) +
    ggplot2::scale_fill_gradient2(
      low = 'red',
      mid = 'black',
      high = 'blue',
      # limits = c(
      #            min(sign(data_frame_xyz$z) * log10(abs(data_frame_xyz$z))),
      #            max(sign(data_frame_xyz$z) * log10(abs(data_frame_xyz$z)))
      #            ),
      # midpoint = 0
    ) +
    ggplot2::coord_fixed(
      ratio = 1,
      xlim = NULL,
      ylim = NULL,
      expand = TRUE,
      clip = "on"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  p
  
  return(p)
}
