% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/explore_mutation_space_v2_helpers.R
\name{sort_new_by_penalised}
\alias{sort_new_by_penalised}
\title{sort_new_by_penalised}
\usage{
sort_new_by_penalised(bitlocus_f, res_df_f, compression)
}
\arguments{
\item{bitlocus_f}{matrix, representing a bitlocus}

\item{res_df_f}{a data frame containing results from a tree search run.}
}
\description{
After computing results, we want to penalize results that needed more steps; e.g, if a depth 1 and a depth 3 result
gave the same result, we want the depth 1 result to appear first. This also helps avoiding over-correcting of the model,
if it e.g. uses 3 steps to correct a small alignment artifact. Each step therefore gets a penalty that is 0.9 times the size
of the smallest dot of a bitlocus.
}
