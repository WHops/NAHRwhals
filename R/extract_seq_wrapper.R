#' A wrapper for extracting sequence
#'
#' This function extracts a sequence from the reference genome and the assembly genome based on the input parameters.
#'
#' @param params A list containing parameters for the sequence extraction
#' @param outlinks A list containing links for the output files
#'
#' @return A list containing the start and end positions of the extracted sequence and the parameters for the sequence extraction
#'
#' @author Wolfram Hoeps
#' @export
extract_sequence_wrapper <- function(params, outlinks) {

  # If we have a pre-computed coarse alignment, then we can use this to find out
  # which region we are talking about.
  if (!is.null(params$conversionpaf_link)) {
    # Pad-sequence
    start_end_pad <- enlarge_interval_by_factor(params$start_x,
      params$end_x,
      params$xpad,
      seqname_f = params$seqname_x,
      conversionpaf_f = params$conversionpaf_link
    )
    start_x_pad <- start_end_pad[1]
    end_x_pad <- start_end_pad[2]

    chr_start_end_pad <- c(params$seqname_x, start_x_pad, end_x_pad)

    # First, write the asm y and hg38 x.

    # First, write the asm y and hg38 x.
    # if (params$samplename_y == "hg38") {
    #   # If we have hg38 as y, we write hg38 into the 'x' and 'y' fasta.
    #   coords = write_x_y_sequences(params$seqname_x,
    #     start_x_pad,
    #     end_x_pad,
    #     outlinks$genome_x_fa_subseq, # hg38 sequence gets directed here
    #     outlinks$genome_y_fa_subseq, # hg38seq gets directed here, TOO
    #     NULL,
    #     NULL,
    #     NULL,
    #     params,n
    #     extract_only = T
    #   )
    # } else {

    coords = write_x_y_sequences(
      params$seqname_x,
      start_x_pad,
      end_x_pad,
      outlinks$genome_x_fa_subseq, # hg38 seq goes here
      outlinks$genome_y_fa_subseq, # asm seq goes here
      params$genome_y_fa,
      params$conversionpaf_link,
      outlinks,
      params
    )
    #}
    
    if (is.null(coords)){
      print("Liftover unsuccessful/inconclusive. Try larger frame!")
      return(NULL)
    }

    if (exists("log_collection")) {
      log_collection$y_seqname <<- coords['new_seqname']
      log_collection$y_start <<- coords['new_x_start']
      log_collection$y_end <<- coords['new_x_end']
    }
    
  } else if (is.null(params$conversionpaf_link)) {
    # Huh?
    run_silent(paste0("cp ", genome_x_fa, " ", outlinks$genome_x_fa_subseq))
    run_silent(paste0("cp ", genome_y_fa, " ", outlinks$genome_y_fa_subseq))

    start_x_pad <- params$start_x
    end_x_pad <- params$end_x
    chr_start_end_pad <- c(params$seqname_x, start_x_pad, end_x_pad)

    run_silent(paste0("cp ", genome_x_fa, " ", outlinks$genome_x_fa_subseq))
    run_silent(paste0("cp ", genome_y_fa, " ", outlinks$genome_y_fa_subseq))
  }
  return(list(chr_start_end_pad, params))
}

