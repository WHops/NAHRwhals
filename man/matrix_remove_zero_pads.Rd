% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrix_remove_zero_pads.R
\name{matrix_remove_zero_pads}
\alias{matrix_remove_zero_pads}
\title{matrix_remove_zero_pads}
\usage{
matrix_remove_zero_pads(bitl_f)
}
\arguments{
\item{bitl_f}{a matrix containing the data to be processed.}
}
\value{
A matrix with zero-padding columns and rows removed. If the resulting matrix is empty, the function returns NULL.
}
\description{
This function removes padding zero columns or rows from a matrix. If all outside columns and rows have entries, the function returns the original matrix.
}
\author{
Wolfram Hoeps
}
