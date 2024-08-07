% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/paf_to_bed.R
\name{paf_write_bed}
\alias{paf_write_bed}
\title{Wrapperfunction to go from a paf file (fresh out of minimap2) to a
filtered paf file (only nondiagonal) and then
to a bedfile.}
\usage{
paf_write_bed(inpaf_link, outsdbed_link)
}
\arguments{
\item{inpaf_link}{A link to the inpaf to be converted to bed.}

\item{outsdbed_link}{A link to the output bedfile. Can be NULL, then
the output will be returned}
}
\value{
a bed dataframe if outsdbed_link is NULL.
}
\description{
The Bedfile will have the wrong colnames
}
\author{
Wolfram Hoeps
}
