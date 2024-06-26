% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrapper_aln_and_analyse.R
\name{determine_plot_minlen}
\alias{determine_plot_minlen}
\title{Determine optimal minimum length for plotting a dot plot}
\usage{
determine_plot_minlen(start, end)
}
\arguments{
\item{start}{The starting coordinate of the input sequence.}

\item{end}{The ending coordinate of the input sequence.}
}
\value{
The optimal minimum length for plotting the dot plot.
}
\description{
This function determines the optimal minimum length for plotting a dot plot based on the length of the input sequence. It uses a set of pre-defined length cutoffs and minimum length values to choose the best minimum length for the input sequence.
}
