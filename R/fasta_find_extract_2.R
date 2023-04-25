#' Write x and y subsequences based on given coordinates.
#'
#' This function extracts subsequences of genome_x_fa and genome_y_fa based on the given coordinates and parameters. It also converts the coordinates from genome_x_fa to genome_y_fa using a liftover_coarse function. The function returns the new start and end coordinates of the regions on genome_x_fa and genome_y_fa as a list. If extract_only is set to true, it only extracts the subsequence of genome_x_fa and returns the new start and end coordinates of the subsequence.
#'
#' @param seqname_x character string specifying the sequence name in genome_x_fa.
#' @param start_x_pad numeric value indicating the start position on genome_x_fa.
#' @param end_x_pad numeric value indicating the end position on genome_x_fa.
#' @param genome_x_fa_subseq character string specifying the output file for the subsequence of genome_x_fa.
#' @param genome_y_fa_subseq character string specifying the output file for the subsequence of genome_y_fa.
#' @param genome_y_fa character string specifying the path to the genome_y_fa file.
#' @param conversionpaf_link character string specifying the path to the conversion PAF file.
#' @param outlinks a list containing output file paths.
#' @param params a list containing parameters for the function.
#' @param extract_only logical value indicating whether or not to only extract the subsequence of genome_x_fa. Default is FALSE.
#'
#' @return a list containing the new sequence name, new start position on genome_x_fa, and new end position on genome_x_fa.
#'
#' @export
#' @author Wolfram Höps
write_x_y_sequences <- function(seqname_x, 
                                start_x_pad, 
                                end_x_pad, 
                                genome_x_fa_subseq, 
                                genome_y_fa_subseq, 
                                genome_y_fa, 
                                conversionpaf_link, 
                                outlinks, 
                                params, 
                                extract_only = F){

  if ((end_x_pad - start_x_pad) < params$maxlen_refine){
    refine_runnr = 1
  } else {
    refine_runnr = 0
  }
  
  if (extract_only == F){
    # Get coordinates in y
    coords_liftover = liftover_coarse(seqname_x,
                                      start_x_pad,
                                      end_x_pad,
                                      conversionpaf_link,
                                      lenfactor = 1, # Unneeded parameter
                                      whole_chr = F,#(params$start_x %in% c(0, 1)),
                                      refine_runnr = refine_runnr)
    if (is.null(coords_liftover)){
      return(NULL)
    }
  }
  

  # Get subseq-fastas in x and y
  extract_subseq_bedtools(params$genome_x_fa,
                          params$seqname_x,
                          start_x_pad,
                          end_x_pad,
                          genome_x_fa_subseq,
                          params)
  
  if (extract_only == T){
    
    extract_subseq_bedtools(params$genome_x_fa,
                            params$seqname_x,
                            start_x_pad,
                            end_x_pad,
                            genome_y_fa_subseq,
                            params)
    
    return( list(new_seqname = seqname_x,
                 new_x_start = start_x_pad,
                 new_x_end = end_x_pad))
    }
  
  extract_subseq_bedtools(genome_y_fa,
                          coords_liftover$lift_contig,
                          coords_liftover$lift_start,
                          coords_liftover$lift_end,
                          genome_y_fa_subseq,
                          params)
  
  if (refine_runnr == 1){
    # Special stuff
    paf = make_chunked_minimap_alnment(
      params,
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
    if (is.null(coords_liftover_2nd)){
      return(NULL)
    }
    
    
    print(coords_liftover)
    print(coords_liftover_2nd)
    
    extract_subseq_bedtools(genome_y_fa,
                                coords_liftover_2nd$lift_contig,
                                coords_liftover_2nd$lift_start,
                                coords_liftover_2nd$lift_end,
                                genome_y_fa_subseq,
                                params)
      
    return( list(new_seqname = coords_liftover_2nd$lift_contig,
                 new_x_start = coords_liftover_2nd$lift_start,
                 new_x_end = coords_liftover_2nd$lift_end))
  }
  
  return( list(new_seqname = coords_liftover$lift_contig,
               new_x_start = coords_liftover$lift_start,
               new_x_end = coords_liftover$lift_end))

}

