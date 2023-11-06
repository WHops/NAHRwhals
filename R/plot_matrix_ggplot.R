#' Compute the signed square root
#'
#' @param x A numeric value or vector
#' @return The signed square root of the input value(s)
#' @export
pmsqrt <- function(x) {
  return(sign(x) * sqrt(abs(x)))
}

#' Compute the reverse signed square root
#'
#' @param x A numeric value or vector
#' @return The reverse signed square root of the input value(s)
#' @export
pmsqrt_rev <- function(x) {
  return(sign(x) * (x**2))
}


#' Plot_matrix_ggplot
#'
#' Make a plot of a dataframe matrix.
#' @param data_frame_xyz a dataframe with three coordinate columns (x,y,z)
#' @author Wolfram Hoeps
#' @export
plot_matrix_ggplot <- function(data_frame_xyz) {
  p <- ggplot2::ggplot(data_frame_xyz) +
    ggplot2::geom_tile(ggplot2::aes(
      x = x,
      y = y,
      fill = sign(z) * log10(abs(z))
    )) +
    ggplot2::scale_fill_gradient2(
      low = "red",
      mid = "black",
      high = "blue",
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


#' Plot a matrix of data
#'
#' @param data_frame_xyz A dataframe with three coordinate columns (x, y, z)
#' @param colnames_f A vector of column names
#' @param rownames_f A vector of row names
#' @return A ggplot2 object representing the matrix plot
#' @author Wolfram Hoeps
#' @export
plot_matrix_ggplot_named <- function(data_frame_xyz, colnames_f, rownames_f) {
  trans_pmsq <- scales::trans_new(
    "trans_pmsq",
    pmsqrt,
    pmsqrt_rev
  )
  # Gridpoint lenghts are calculated here
  diff_rownames <- paste0(as.character(diff(rownames_f)), " (g", 1:length(diff(rownames_f)), ")")
  diff_colnames <- paste0(as.character(diff(colnames_f)), " (g", 1:length(diff(colnames_f)), ")")

  missing_y = which(!(1:max(data_frame_xyz$y) %in% data_frame_xyz$y))
  missing_x = which(!(1:max(data_frame_xyz$x) %in% data_frame_xyz$x))

  for (i in c(missing_y)){
    data_frame_xyz = rbind(data_frame_xyz, c(i, i, NA))
  }

  for (i in c(missing_x)){
    data_frame_xyz = rbind(data_frame_xyz, c(i, i, NA))
  }

  # Now they are put in place in x and y
  data_frame_xyz$x <- as.factor(as.character(diff_colnames[data_frame_xyz$x]))
  data_frame_xyz$y <- as.factor(as.character(diff_rownames[data_frame_xyz$y]))

  data_frame_xyz = data_frame_xyz[order(data_frame_xyz$x, data_frame_xyz$y),]

  # limit <- max(sqrt(abs(data_frame_xyz$z))) * c(-1, 1)
  # limits = unique(sort(c(sort(unique(sqrt(sort(abs(data_frame_xyz$z))))),
  #                -sort(unique(sqrt(sort(abs(data_frame_xyz$z))))))))
  data_frame_xyz$z <- data_frame_xyz$z / 1000
  limits <- unique(sort(c(
    sort(unique((sort(abs(data_frame_xyz$z))))),
    -sort(unique((sort(abs(data_frame_xyz$z)))))
  )))

  # limits = limits[abs(limits) > 1000]n
  # Make the plot
  p <- ggplot2::ggplot(data_frame_xyz) +
    ggplot2::geom_tile(ggplot2::aes(
      x = x,
      y = y,
      fill = z # sign(z) * sqrt(abs(z))
    )) +
    ggplot2::scale_fill_stepsn(
      "Blocksize [kbp]",
      # colors=c('#b2182b','#ef8a62','#fddbc7',
      #         '#f7f7f7','#d1e5f0','#67a9cf','#2166ac'),
      colors = c(
        "#b2182b", "#ef8a62", "#fddbc7",
        "#d1e5f0", "#67a9cf", "#2166ac"
      ),
      # limits = c(min(limits),max(limits)),
      n.breaks = 8,
      limits = c(min(limits), max(limits)),
      show.limits = T, # ,
      # labels=round(seq(limits[1], limits[length(limits)],length.out=8)),
      trans = trans_pmsq
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
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 3),
      axis.text.y = ggplot2::element_text(size = 3)
    )

  # ... and sort the x and y axes
  p <- p + ggplot2::scale_x_discrete(limits = as.character(diff_colnames)) +
    ggplot2::scale_y_discrete(limits = as.character(diff_rownames)) +
    ggplot2::labs(x = "target sequence", y = "query sequence")


  p


  return(p)
}


#' Create a segmented pairwise plot with grid lines
#'
#' This function takes a grid_xy and a plot, and adds segments to the plot based on the grid lines. The resulting plot is saved to a file. This function is considered quick-and-dirty and may need to be improved in the future.
#'
#' @param grid_xy A list containing two vectors of x and y values that represent the grid lines
#' @param plot_x_y A ggplot2 object to which the grid segments will be added
#' @param outlinks A list of output file paths
#'
#' @author Wolfram Hoeps
#' @export
make_segmented_pairwise_plot <- function(grid_xy, plot_x_y, outlinks) {

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
  
  ggplot2::ggsave(
    filename = paste0(outlinks$outfile_colored_segment),
    plot = plot_x_y_segmented,
    width = 15,
    height = 15,
    units = "cm",
    dpi = 300
  )
}




#' Create a modified grid plot
#'
#' @param res A data frame with results
#' @param gridmatrix A matrix representing the grid
#' @param outlinks A list containing output links
#' @export
make_modified_grid_plot <- function(res, gridmatrix, outlinks) {
  
  res$mut_max = NULL
  # Sort results by 'eval' in descending order
  res <- res[order(res$eval, decreasing = TRUE), ]
  # Modify the grid matrix using the first row of results
  grid_modified <- modify_gridmatrix(gridmatrix, res[1, ])

  # Calculate gridlines for x and y axis
  gridlines.x <- cumsum(c(0, as.numeric(colnames(grid_modified))))
  gridlines.y <- cumsum(c(0, as.numeric(row.names(grid_modified))))

  # Update column and row names of the modified grid
  colnames(grid_modified) <- seq(1, dim(grid_modified)[2])
  row.names(grid_modified) <- seq(1, dim(grid_modified)[1])

  # Melt the modified grid matrix
  gm2 <- reshape2::melt(grid_modified)
  colnames(gm2) <- c("y", "x", "z")

  # Create the modified grid plot using non-zero values
  grid_mut_plot <- plot_matrix_ggplot_named(gm2[gm2$z != 0, ], gridlines.x, gridlines.y)

  # Save the grid_mut_plot to a file
  ggplot2::ggsave(
    filename = paste0(outlinks$outfile_plot_grid_mut),
    plot = grid_mut_plot,
    width = 10,
    height = 10,
    units = "cm",
    dpi = 300
  )

  # Print the grid_mut_plot
  print(grid_mut_plot)
}
