% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrapper_aln_and_analyse.R
\name{determine_xpad}
\alias{determine_xpad}
\title{Determine desired padding value for the x-axis of a dot plot}
\usage{
determine_xpad(start, end)
}
\arguments{
\item{start}{The starting coordinate of the input sequence.}

\item{end}{The ending coordinate of the input sequence.}
}
\value{
The desired padding value for the x-axis of the dot plot.
}
\description{
This function determines the optimal padding value for the x-axis of a dot plot based on the length of the input sequence. It uses a set of pre-defined length cutoffs and padding values to choose the best padding value for the input sequence.
}
