#' little stuff
#' @export
pmsqrt <- function(x){
  return (sign(x) * sqrt(abs(x)))
}

#' little stuff
#' @export
pmsqrt_rev <- function(x){
  return(sign(x) * (x**2))
}


#' Plot_matrix_ggplot
#'
#' Make a plot of a dataframe matrix.
#' @param data_frame_xyz a dataframe with three coordinate columns (x,y,z)
#' @author Wolfram Hoeps
#' @export
plot_matrix_ggplot <- function(data_frame_xyz) {
  
  browser()
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

  trans_pmsq = scales::trans_new(
    'trans_pmsq',
    pmsqrt,
    pmsqrt_rev
  )
  
  # Gridpoint lenghts are calculated here
  diff_rownames = paste0(as.character(diff(rownames_f)), ' (g',1:length(diff(rownames_f)), ')')
  diff_colnames = paste0(as.character(diff(colnames_f)), ' (g',1:length(diff(colnames_f)), ')')

  # Now they are put in place in x and y 
  data_frame_xyz$x = as.factor(as.character(diff_colnames[data_frame_xyz$x]))
  data_frame_xyz$y = as.factor(as.character(diff_rownames[data_frame_xyz$y]))
  
  
  #limit <- max(sqrt(abs(data_frame_xyz$z))) * c(-1, 1)
  limits = unique(sort(c(sort(unique(sqrt(sort(abs(data_frame_xyz$z))))),
                  -sort(unique(sqrt(sort(abs(data_frame_xyz$z))))))))
  #limits = limits[abs(limits) > 1000]n
  # Make the plot
  p = ggplot2::ggplot(data_frame_xyz) + ggplot2::geom_tile(ggplot2::aes(
    x = x,
    y = y,
    fill = sign(z) * sqrt(abs(z))
  )) +
    # ggplot2::scale_fill_steps2(
    #   low = 'red',
    #   high = 'blue',
    #   #midpoint = 0,
    #   limits=c(min(limits), max(limits)),
    #   breaks=limits,
    #   space='Lab'
    # ) +
    ggplot2::scale_fill_stepsn(
      #colors=c('#b2182b','#ef8a62','#fddbc7',
      #         '#f7f7f7','#d1e5f0','#67a9cf','#2166ac'),
      colors=c('#b2182b','#ef8a62','#fddbc7',
               '#d1e5f0','#67a9cf','#2166ac'),
      #limits = c(min(limits),max(limits)),
      n.breaks = 8,
      limits = c(min(limits), max(limits)),
      show.limits=T,
      #labels=round(seq(limits[1], limits[length(limits)],length.out=8)),
      trans=trans_pmsq
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
       axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1, size=3),
       axis.text.y = ggplot2::element_text(size=3)
     )
  
  # ... and sort the x and y axes
  p = p + ggplot2::scale_x_discrete(limits=as.character(diff_colnames)) +
      ggplot2::scale_y_discrete(limits=as.character(diff_rownames)) +
      ggplot2::labs(x='target sequence', y='query sequence')
  

  p
  
  
  return(p)
  
}

