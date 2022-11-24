#' write_x_y_sequences
#' This was previously part of the wrapper_aln_and_analyse. 
#' However, over time, this grew into more and more complex code, 
#' So it is for sure worth its own function by now. 
#' @export
write_x_y_sequences <- function(seqname_x, 
                                start_x_pad, 
                                end_x_pad, 
                                genome_x_fa_subseq, 
                                genome_y_fa_subseq, 
                                genome_y_fa, 
                                conversionpaf_link, 
                                outlinks, 
                                params){

  if ((end_x_pad - start_x_pad) < params$maxlen_refine){
    refine_runnr = 1
  } else {
    refine_runnr = 0
  }
  # Get coordinates in y
  coords_liftover = liftover_coarse(seqname_x,
                                    start_x_pad,
                                    end_x_pad,
                                    conversionpaf_link,
                                    lenfactor = 1, # Unneeded parameter
                                    whole_chr = (params$start_x %in% c(0, 1)),
                                    refine_runnr = refine_runnr)
  
  
  
  # Get subseq-fastas in x and y
  extract_subseq_bedtools(params$genome_x_fa,
                          params$seqname_x,
                          start_x_pad,
                          end_x_pad,
                          genome_x_fa_subseq)
  
  extract_subseq_bedtools(genome_y_fa,
                          coords_liftover$lift_contig,
                          coords_liftover$lift_start,
                          coords_liftover$lift_end,
                          genome_y_fa_subseq)
  
  if (refine_runnr == 1){
    # Special stuff
    paf = make_chunked_minimap_alnment(
      genome_x_fa_subseq,
      genome_y_fa_subseq,
      outlinks$outpaf_link_x_y,
      chunklen = params$chunklen,
      minsdlen = params$plot_minlen,
      saveplot = F,
      hllink = F,
      hltype = F,
      hlstart = F,#start_x - start_x_pad,
      hlend = F,
      x_start = start_x_pad,
      x_end = end_x_pad,
      x_seqname = params$seqname_x,
      anntrack = params$anntrack,
      hltrack = params$hltrack,
      onlypafreturn = T
    )
    coords_liftover_2nd = liftover_coarse('None',
                                          'none',
                                          'nonepaflink', 
                                          paf,
                                          refine_runnr = 2)
    
    print(coords_liftover)
    print(coords_liftover_2nd)
    
  extract_subseq_bedtools(genome_y_fa,
                              coords_liftover_2nd$lift_contig,
                              coords_liftover_2nd$lift_start,
                              coords_liftover_2nd$lift_end,
                              genome_y_fa_subseq)
    
  return( list(new_seqname = coords_liftover_2nd$lift_contig,
               new_x_start = coords_liftover_2nd$lift_start,
               new_x_end = coords_liftover_2nd$lift_end))
  }
  
  return( list(new_seqname = coords_liftover$lift_contig,
               new_x_start = coords_liftover$lift_start,
               new_x_end = coords_liftover$lift_end))

}

