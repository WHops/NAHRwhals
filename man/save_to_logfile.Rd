% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/save_logfile.R
\name{save_to_logfile}
\alias{save_to_logfile}
\title{Save results to logfile.}
\usage{
save_to_logfile(log, res, logfile, params, alt_x = F)
}
\arguments{
\item{log}{A named list containing the metadata of the region analyzed, including the chromosome, start and end positions, the name of the sample, the width of the original region, the padding added to the x-axis, and several binary flags indicating whether certain tests were performed or not.}

\item{res}{A data frame containing the results of the analysis, including the evaluated score for each mutation and its position.}

\item{logfile}{The file path to which the results should be appended.}

\item{params}{A named list containing the parameters used for the analysis, including the name of the chromosome, the start and end positions, the name of the sample, the compression factor, and the mutations that were simulated and tested.}

\item{alt_x}{A logical indicating whether to include the alternative x-axis coordinates in the log file.}
}
\description{
This function saves the results and parameters to a log file, in a tab-separated format.
}
\author{
Wolfram Hoeps
}
