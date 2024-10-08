% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_whole_genome_list.R
\name{extract_test_list_from_paf}
\alias{extract_test_list_from_paf}
\title{extract_test_list_from_paf}
\usage{
extract_test_list_from_paf(
  all_vs_all_paf,
  out_dir,
  out_file,
  genome_path,
  bedtools_bin,
  merge_distance,
  indel_ignore_distance,
  exclusion_mask = "none"
)
}
\arguments{
\item{all_vs_all_paf}{Path to the paf file containing all-vs-all alignments}

\item{out_dir}{Path to the output directory}

\item{genome_path}{Path to the genome fasta}

\item{bedtools_bin}{Path to the bedtools binary}

\item{merge_distance}{Distance to merge overlapping intervals}
}
\value{
Nothing
}
\description{
extract_test_list_from_paf
}
