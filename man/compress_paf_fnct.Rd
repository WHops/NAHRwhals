% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compress_paf.R
\name{compress_paf_fnct}
\alias{compress_paf_fnct}
\title{Melt Alignment Chunks Back Together}
\usage{
compress_paf_fnct(
  inpaf_link = NULL,
  outpaf_link = NULL,
  inpaf_df = NULL,
  save_outpaf = T,
  n_quadrants_per_axis = NULL,
  second_run = F,
  inparam_chunklen = NULL,
  inparam_compression = NULL,
  qname = NULL
)
}
\arguments{
\item{inpaf_link}{A file path or link to the input chunked PAF file.}

\item{outpaf_link}{A file path or link to the output "molten" PAF file.}

\item{inpaf_df}{An optional input PAF data.frame.}

\item{save_outpaf}{Logical value whether to save the output PAF file (default: TRUE).}

\item{n_quadrants_per_axis}{The number of quadrants per axis (determines chunking).}

\item{second_run}{Logical value whether this is the second run of the function (default: FALSE).}

\item{inparam_chunklen}{The input parameter for chunk length in bp.}

\item{inparam_compression}{The input parameter for compression in bp.}
}
\description{
This function is a crucial part of the process of identifying synteny blocks (SDs)
and is the last step of the chunk-minimap2-edit-melt workflow. It takes the chunked PAF file
produced earlier in the pipeline and "melts" the alignment pieces back together.
Note that improvements may be needed in the future, as the melting process is not
entirely straightforward.
}
\author{
Wolfram Hoeps
}
