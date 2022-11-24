
#' Produce up to three minimap2 based alignments. 
#' @author Wolfram Höps
#' @export
produce_pairwise_alignments_minimap2 <- function(params, outlinks, start_x_pad, end_x_pad){
  # Run alignments.
  # Run REF self alignment only if it hasn't been run before.
  if (params$self_plots) {
    if (is.na(file.size(outlinks$outfile_plot_self_x))) {
      plot_self_x = make_chunked_minimap_alnment(
        outlinks$genome_x_fa_subseq,
        outlinks$genome_x_fa_subseq,
        outlinks$outpaf_link_self_x,
        chunklen = params$chunklen,
        minsdlen = params$plot_minlen,
        saveplot = F,
        hllink =   F,
        hltype =   F,
        hlstart = params$start_x - start_x_pad,
        hlend =   params$end_x - start_x_pad,
        x_start = start_x_pad,
        x_end =   end_x_pad,
        x_seqname = params$seqname_x,
        anntrack = params$anntrack,
        hltrack = params$hltrack
      )
      print(plot_self_x)
      # Save alignment
      save_plot_custom(plot_self_x, outlinks$outfile_plot_self_x, 'pdf')
      save_plot_custom(plot_self_x,
                       outlinks$outfile_plot_self_x,
                       'png',
                       width = 20,
                       height = 20)
      
      
      
      # # # Run y self alignment
      plot_self_y = make_chunked_minimap_alnment(
        outlinks$genome_y_fa_subseq,
        outlinks$genome_y_fa_subseq,
        outlinks$outpaf_link_self_y,
        chunklen = params$chunklen,
        minsdlen = params$plot_minlen,
        saveplot = F,
        hllink = F,
        hltype = F,
        hlstart = F,#NULL,
        hlend = F#NULL
      )
      save_plot_custom(plot_self_y, outlinks$outfile_plot_self_y, 'pdf')
      save_plot_custom(plot_self_y,
                       outlinks$outfile_plot_self_y,
                       'png',
                       width = 20,
                       height = 20)
      print(plot_self_y)
    }
  }
  #Run xy alignment
  plot_x_y = make_chunked_minimap_alnment(
    outlinks$genome_x_fa_subseq,
    outlinks$genome_y_fa_subseq,
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
    hltrack = params$hltrack#end_x - start_x_pad
  )
  # Save alignments
  print(plot_x_y)
  save_plot_custom(plot_x_y, outlinks$outfile_plot_x_y, 'pdf')
  save_plot_custom(plot_x_y,
                   outlinks$outfile_plot_x_y,
                   'png',
                   width = 20,
                   height = 20)
  
  # We want to use plot_x_y later again (maybe)
  return(plot_x_y)
  
}