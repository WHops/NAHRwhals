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
  
  return(p)
  
  
}


#' plot_matrix_ggplot_named
#'
#' Make a plot of a dataframe matrix.
#' @param data_frame_xyz a dataframe with three coordinate columns (x,y,z)
#' @author Wolfram Hoeps
#' @export
plot_matrix_ggplot_named <- function(data_frame_xyz, colnames_f, rownames_f) {

  # Gridpoint lenghts are calculated here
  diff_rownames = paste0(as.character(diff(rownames_f)), ' (g',1:length(diff(rownames_f)), ')')
  diff_colnames = paste0(as.character(diff(colnames_f)), ' (g',1:length(diff(colnames_f)), ')')

  # Now they are put in place in x and y 
  data_frame_xyz$x = as.factor(as.character(diff_colnames[data_frame_xyz$x]))
  data_frame_xyz$y = as.factor(as.character(diff_rownames[data_frame_xyz$y]))
  
  # Make the plot
  p = ggplot2::ggplot(data_frame_xyz) + ggplot2::geom_tile(ggplot2::aes(
    x = x,
    y = y,
    fill = sign(z) * log10(abs(z))
  )) +
    ggplot2::scale_fill_gradient2(
      low = 'red',
      mid = 'black',
      high = 'blue',
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
       panel.grid.minor = ggplot2::element_blank(),
       axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)
     )
  
  # ... and sort the x and y axes
  p = p + ggplot2::scale_x_discrete(limits=as.character(diff_colnames)) +
      ggplot2::scale_y_discrete(limits=as.character(diff_rownames)) +
      ggplot2::labs(x='target sequence', y='query sequence')
  

  
  return(p)
  
}

