

#' wrapper_aln_and_analyse
#' Description ... 
#' @export
wrapper_aln_and_analyse <- function(seqname_x,
                                    start_x,
                                    end_x,
                                    genome_x_fa,
                                    genome_y_fa,
                                    conversionpaf_link,
                                    chunklen = 1000,
                                    aln_pad_factor = 1,
                                    sd_minlen = 1000,
                                    compression = 1000,
                                    depth = 2,
                                    runname = 'test',
                                    include_grid = T,
                                    xpad = 1){
  
  
  sequence_name_output = paste(paste0('res/',seqname_x), format(start_x, scientific = F),  format(end_x, scientific = F), sep='-')
  dir.create('res')
  dir.create(sequence_name_output)
  dir.create(paste0(sequence_name_output, '/self'))
  dir.create(paste0(sequence_name_output, '/self/pdf'))
  dir.create(paste0(sequence_name_output, '/self/paf'))
  dir.create(paste0(sequence_name_output, '/diff'))
  dir.create(paste0(sequence_name_output, '/diff/pdf'))
  dir.create(paste0(sequence_name_output, '/diff/pdf/grid'))
  dir.create(paste0(sequence_name_output, '/diff/paf'))
  dir.create(paste0(sequence_name_output, '/fasta'))
  
  # Define output files
  outpaf_link_self_x =  paste0(sequence_name_output, '/self/paf/', runname, '_x.paf')
  outpaf_link_self_y =  paste0(sequence_name_output, '/self/paf/', runname, '_y.paf')
  outpaf_link_x_y =     paste0(sequence_name_output, '/diff/paf/', runname, '_xy.paf')
  
  res_table_xy =        paste0(sequence_name_output, '/diff/', runname, '_res.tsv')
  
  outfile_plot_self_x = paste0(sequence_name_output, '/self/pdf/', runname, '_x.pdf')
  outfile_plot_self_y = paste0(sequence_name_output, '/self/pdf/', runname, '_y.pdf')
  outfile_plot_x_y =    paste0(sequence_name_output, '/diff/pdf/', runname, '_x_y.pdf')
  
  outfile_plot_pre_grid = paste0(sequence_name_output, '/diff/pdf/grid/', runname, '_x_y_grid_pre.pdf')
  outfile_plot_grid =     paste0(sequence_name_output, '/diff/pdf/grid/', runname, '_x_y_grid.pdf')
  outfile_plot_grid_mut = paste0(sequence_name_output, '/diff/pdf/grid/', runname, '_x_y_grid_mut.pdf')
  
  genome_x_fa_subseq = paste0(sequence_name_output, '/fasta/', runname, '_x.fa')
  genome_y_fa_subseq = paste0(sequence_name_output, '/fasta/', runname, '_y.fa')
  
  # Pad-sequence
  start_end_pad = enlarge_interval_by_factor(start_x, end_x, xpad, seqname_f = seqname_x, conversionpaf_f = conversionpaf_link) 
  start_x_pad = start_end_pad[1]
  end_x_pad = start_end_pad[2]
  browser()
  # Get coordinates in y
  coords_liftover = liftover_coarse(seqname_x, start_x_pad, end_x_pad, conversionpaf_link, lenfactor = aln_pad_factor)
  
  # Get subseq-fastas in x and y
  extract_subseq_bedtools(genome_x_fa, seqname_x, start_x_pad, end_x_pad, genome_x_fa_subseq)
  extract_subseq_bedtools(genome_y_fa, coords_liftover$lift_contig, coords_liftover$lift_start, coords_liftover$lift_end, genome_y_fa_subseq)
  
  # Run alignments. 
  plot_self_x = make_chunked_minimap_alnment(genome_x_fa_subseq, genome_x_fa_subseq, outpaf_link_self_x,
                                             chunklen = chunklen, minsdlen = 2000, saveplot=F,
                                             hllink = F, hltype = F, hlstart = start_x - start_x_pad, hlend = end_x - start_x_pad)
  
  plot_self_y = make_chunked_minimap_alnment(genome_y_fa_subseq, genome_y_fa_subseq, outpaf_link_self_y,
                                             chunklen = chunklen, minsdlen = 2000, saveplot=F,
                                             hllink = F, hltype = F, hlstart = start_x - start_x_pad, hlend = end_x - start_x_pad)
  
  plot_x_y = make_chunked_minimap_alnment(genome_x_fa_subseq, genome_y_fa_subseq, outpaf_link_x_y,
                                             chunklen = chunklen, minsdlen = 2000, saveplot=F,
                                             hllink = F, hltype = F, hlstart = start_x - start_x_pad, hlend = end_x - start_x_pad)
  
  # Save orig alignments
  ggplot2::ggsave(filename = outfile_plot_self_x,
                  plot = plot_self_x, 
                  width = 20, 
                  height = 20, 
                  units = 'cm',
                  dpi = 300)
  ggplot2::ggsave(filename = outfile_plot_self_y,
                  plot = plot_self_y, 
                  width = 20, 
                  height = 20, 
                  units = 'cm',
                  dpi = 300)
  ggplot2::ggsave(filename = outfile_plot_x_y,
                  plot = plot_x_y, 
                  width = 20, 
                  height = 20, 
                  units = 'cm',
                  dpi = 300)
  
  if (include_grid){
    # Make an xy grid
    grid_xy = wrapper_paf_to_bitlocus(outpaf_link_x_y, minlen = sd_minlen, compression = compression,
                                      gridplot_save = outfile_plot_grid, pregridplot_save = outfile_plot_pre_grid )
    gridmatrix = gridlist_to_gridmatrix(grid_xy)
    
    res = explore_mutation_space(gridmatrix, depth = depth)

    # Make a grid after applying the top res
    grid_modified = modify_gridmatrix(gridmatrix, res[1,])
    gm2 = reshape2::melt(grid_modified)
    colnames(gm2) = c('x','y','z')
    grid_mut_plot = plot_matrix_ggplot(gm2[gm2$z != 0,])
    ggplot2::ggsave(filename = outfile_plot_grid_mut,
                    plot = grid_mut_plot,
                    width = 10,
                    height = 10,
                    units = 'cm',
                    dpi = 300)
    
    # Save res table
    write.table(res, file = res_table_xy,
                col.names = T,
                row.names = F,
                quote = F,
                sep='\t'
    )
  }


  system(paste0('rm ', genome_y_fa_subseq))

  
  
}
