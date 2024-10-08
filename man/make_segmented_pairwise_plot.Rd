% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_matrix_ggplot.R
\name{make_segmented_pairwise_plot}
\alias{make_segmented_pairwise_plot}
\title{Create a segmented pairwise plot with grid lines}
\usage{
make_segmented_pairwise_plot(grid_xy, plot_x_y, outlinks)
}
\arguments{
\item{grid_xy}{A list containing two vectors of x and y values that represent the grid lines}

\item{plot_x_y}{A ggplot2 object to which the grid segments will be added}

\item{outlinks}{A list of output file paths}
}
\description{
This function takes a grid_xy and a plot, and adds segments to the plot based on the grid lines. The resulting plot is saved to a file. This function is considered quick-and-dirty and may need to be improved in the future.
}
\author{
Wolfram Hoeps
}
