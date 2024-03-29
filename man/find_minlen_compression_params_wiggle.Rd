% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wiggle_find_parameters.R
\name{find_minlen_compression_params_wiggle}
\alias{find_minlen_compression_params_wiggle}
\title{Find parameters by wiggling.}
\usage{
find_minlen_compression_params_wiggle(inpaf, compression_params)
}
\arguments{
\item{inpaf}{a string representing the path to a PAF file.}

\item{compression_params}{a list containing parameters for the compression.
\itemize{
\item mode: character string, either "precise" or "compressed".
\item baseline_log_minsize_min: integer representing the minimum value of the log2
of the candidate minlen values.
\item baseline_log_minsize_max: integer representing the maximum value of the log2
of the candidate minlen values.
\item n_tests: integer representing the number of tests to run for each candidate
minlen and compression value.
\item max_n_alns: integer representing the maximum number of alignments that
should be included in a grid.
\item max_size_col_plus_rows: integer representing the maximum number of columns
plus the maximum number of rows that a grid should have.
}}
}
\value{
a list with the optimal minlen and compression values.
}
\description{
Finds optimal parameters for compression by performing n tests for each candidate
minlen and compression value in a range defined by the input parameters. The
optimal values are the ones that yield the lowest collapse penalty and are
successful (no incongruencies), and the size of the grid is not larger than the
value specified by the max_size_col_plus_rows parameter. If the function cannot
find at least 2 viable grids, it raises the criteria (max_n_alns and
max_size_col_plus_rows) and retries up to n rounds. If the function cannot find
any viable grids after n rounds, it returns the suboptimal values with the
smallest number of rows and columns.
}
