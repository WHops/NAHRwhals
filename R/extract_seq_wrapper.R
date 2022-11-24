
#' A wrapper for extracting sequence
#' @author Wolfram Höps
#' @export
extract_sequence_wrapper <- function(params, outlinks, start_x_pad, start_y_pad){
  # If we have a pre-computed coarse alignment, then we can use this to find out 
  # which region we are talking about. 
  if (!is.null(params$conversionpaf_link)) {
    
    # Pad-sequence
    print('hi')
    start_end_pad = enlarge_interval_by_factor(params$start_x,
                                               params$end_x,
                                               params$xpad,
                                               seqname_f = params$seqname_x,
                                               conversionpaf_f = params$conversionpaf_link)
    start_x_pad = start_end_pad[1]
    end_x_pad = start_end_pad[2]
    
    # First, write the asm y and hg38 x. 
    write_x_y_sequences(params$seqname_x, 
                        start_x_pad, 
                        end_x_pad, 
                        outlinks$genome_x_fa_subseq, 
                        outlinks$genome_y_fa_subseq, 
                        params$genome_y_fa,
                        params$conversionpaf_link,
                        outlinks, 
                        params)   
    
    if (params$alt_ref_sample != F){
      
      # Run a second time, this tome overwriting the x sequence!
      
      print('Detected "alt_ref_sample" != F. Using an alternative sequence to plot on x-axis)')
      
      # Second, search for the T2T region. 
      new_coords = write_x_y_sequences(params$seqname_x, start_x_pad, end_x_pad, 
                                       outlinks$genome_x_fa_altref_subseq, 
                                       outlinks$genome_x_fa_subseq, 
                                       params$genome_alt_ref_fa,
                                       params$conversionpaf_alt_ref_link,
                                       outlinks, 
                                       params)    
      
      params$seqname_x = new_coords[['new_seqname']]
      start_x_pad = new_coords[['new_x_start']]
      end_x_pad = new_coords[['new_x_end']]
      
      #params[c('anntrack', 'hltrack') ] = NULL
      
    }
    
    
  } else {
    
    # Huh? 
    system(paste0('cp ', genome_x_fa, ' ', outlinks$genome_x_fa_subseq))
    system(paste0('cp ', genome_y_fa, ' ', outlinks$genome_y_fa_subseq))
    
    start_x_pad = params$start_x
    end_x_pad = params$end_x
    
    system(paste0('cp ', genome_x_fa, ' ', outlinks$genome_x_fa_subseq))
    system(paste0('cp ', genome_y_fa, ' ', outlinks$genome_y_fa_subseq))
    
  }
}