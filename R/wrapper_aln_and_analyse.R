

#' wrapper_aln_and_analyse
#' @description The main wrapper function for running nahrtoolkit. 
#' This is where all the information flow together. 
#' @param seqname_x region of interest on ref: seqname [character]
#' @param start_x region of interest on ref: start coord [numeric]
#' @param end_x region of interest on ref: end coord [numeric]
#' @param genome_x_fa link to genome assembly fasta file serving as 'reference' (Typically hg38) [character/link]
#' @param genome_y_fa link to fasta file serving as 'query'. Typically a genome assembly [character/link]
#' @param conversionpaf_link link to paf containing a genome-wide alignment between genome_x_fa and genome_y_fa [character/link]
#' @param logfile link to output 'log'-file [character/link]
#' @param chunklen Length of sequence chunks used in fine alignment. Typically between 100 and 50k. Increase chunklen for faster runtime, sacrificing fine alignment details. [numeric]
#' @param sd_minlen Auto SD-detection: shortest alignment to consider for SD detection. Typically same number as compression parameter. [numeric]
#' @param compression Auto SD-detection: Round SD locations by maximum this number of basepairs. Typically same as sd_minlen. Should be magnitudes smaller than event size. [numeric]
#' @param depth Maximum number of chained mutations to model. Current maximum is 3. Has a big impact on runtime. [numeric]
#' @param xpad enlarge input sequence by a factor. e.g. region chr1:1000-2000 will turn to 500 - 2500 with xpad=2 (region is twice as large afterwards.) [numeric]
#' @param samplename Name of the sample/run. Used in the output filenames. [character]
#' @param include_grid T/F: do you want to auto-detect SDs and model mutations? [bool; T/F]
#' @param debug debug mode used for developing (default: F) [bool; T/F]
#' @param use_paf_library T/F: If this is off, input region is ignored, and alignment is computed on entire x vs y sequence. Use this if your input fastas contain only the region of interest. (default: F) [bool; T/F]
#' @param aln_pad_factor Factor to enlarge the query sequence by a factor. [numeric]
#' @param clean_all Should fastas be cleared? (default: T). Advisable for very large sequences [T/F]
#' @examples 
#' 
#' Extract region chr1:500k-600k in hg38, find the respective region in the assembly, 
#' make an alignment and print the results to res/chr1-500000-6000000/
#' wrapper_aln_and_analyse('chr1', 500000, 600000, 'my/genomes/hg38.fa', 'my/assembly.fa', 
#'                         'my/conversionpaf.fa', samplename='testsample', 
#'                         chunklen = 1000, sd_minlen = 100, compression = 100, 
#'                         depth = 2, xpad = 1.2)
#' @export
wrapper_aln_and_analyse <- function(seqname_x,
                                    start_x,
                                    end_x,
                                    genome_x_fa,
                                    genome_y_fa,
                                    conversionpaf_link,
                                    logfile,
                                    chunklen = 1000,
                                    aln_pad_factor = 1.0,
                                    sd_minlen = 1000,
                                    compression = 1000,
                                    depth = 2,
                                    samplename = 'test',
                                    include_grid = T,
                                    xpad = 1,
                                    debug = F,
                                    use_paf_library = T,
                                    clean_all = T){
  
  if (debug){
    print('Debug mode!')
    browser()
  }
  log_collection <<- init_log_with_def_values()
  log_collection[c('chr', 'start', 'end', 'xpad', 'chunklen', 'samplename', 'depth')] <<-
    c(seqname_x, start_x, end_x, xpad, chunklen, samplename, depth)
  
  

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
  outpaf_link_self_x =  paste0(sequence_name_output, '/self/paf/aln_ref', samplename,'.paf')
  outpaf_link_self_y =  paste0(sequence_name_output, '/self/paf/', samplename, '_y.paf')
  outpaf_link_x_y =     paste0(sequence_name_output, '/diff/paf/', samplename, '_xy.paf')
  
  res_table_xy =        paste0(sequence_name_output, '/diff/', samplename, '_res.tsv')
  
  outfile_plot_self_x = paste0(sequence_name_output, '/self/pdf/aln_ref.pdf')
  outfile_plot_self_y = paste0(sequence_name_output, '/self/pdf/', samplename, '_y.pdf')
  outfile_plot_x_y =    paste0(sequence_name_output, '/diff/pdf/', samplename, '_x_y.pdf')
  
  outfile_plot_pre_grid = paste0(sequence_name_output, '/diff/pdf/grid/', samplename, '_x_y_grid_pre.pdf')
  outfile_plot_grid =     paste0(sequence_name_output, '/diff/pdf/grid/', samplename, '_x_y_grid.pdf')
  outfile_plot_grid_mut = paste0(sequence_name_output, '/diff/pdf/grid/', samplename, '_x_y_grid_mut.pdf')
  
  genome_x_fa_subseq = paste0(sequence_name_output, '/fasta/', samplename, '_x.fa')
  genome_y_fa_subseq = paste0(sequence_name_output, '/fasta/', samplename, '_y.fa')
  
  if (use_paf_library){
    
    # Pad-sequence
    start_end_pad = enlarge_interval_by_factor(start_x, end_x, xpad, seqname_f = seqname_x, conversionpaf_f = conversionpaf_link) 
    start_x_pad = start_end_pad[1]
    end_x_pad = start_end_pad[2]
    
    # Get coordinates in y
    coords_liftover = liftover_coarse(seqname_x, start_x_pad, end_x_pad, conversionpaf_link, lenfactor = aln_pad_factor)
  
    # Get subseq-fastas in x and y
    extract_subseq_bedtools(genome_x_fa, seqname_x, start_x_pad, end_x_pad, genome_x_fa_subseq)
    extract_subseq_bedtools(genome_y_fa, coords_liftover$lift_contig, coords_liftover$lift_start, coords_liftover$lift_end, genome_y_fa_subseq)
  } else {
    system(paste0('cp ', genome_x_fa, ' ', genome_x_fa_subseq))
    system(paste0('cp ', genome_y_fa, ' ', genome_y_fa_subseq))
    
    start_x_pad = 0
    end_x_pad = 1
    start_x = 0
    end_x = 1
    system(paste0('cp ', genome_x_fa, ' ', genome_x_fa_subseq))
    system(paste0('cp ', genome_y_fa, ' ', genome_y_fa_subseq))
    
  }
  # Run alignments. 
  # Run REF self alignment only if it hasn't been run before. 
  if (is.na(file.size(outfile_plot_self_x))){
    plot_self_x = make_chunked_minimap_alnment(genome_x_fa_subseq, genome_x_fa_subseq, outpaf_link_self_x,
                                               chunklen = chunklen, minsdlen = 2000, saveplot=F,
                                               hllink = F, hltype = F, hlstart = start_x - start_x_pad, hlend = end_x - start_x_pad)
    
    # Save alignment
    ggplot2::ggsave(filename = outfile_plot_self_x,
                    plot = plot_self_x, 
                    width = 20, 
                    height = 20, 
                    units = 'cm',
                    dpi = 300)
  }
  
  # Run y self alignment
  plot_self_y = make_chunked_minimap_alnment(genome_y_fa_subseq, genome_y_fa_subseq, outpaf_link_self_y,
                                             chunklen = chunklen, minsdlen = 2000, saveplot=F,
                                             hllink = F, hltype = F, hlstart = start_x - start_x_pad, hlend = end_x - start_x_pad)
  
  # Run xy alignment
  plot_x_y = make_chunked_minimap_alnment(genome_x_fa_subseq, genome_y_fa_subseq, outpaf_link_x_y,
                                             chunklen = chunklen, minsdlen = 2000, saveplot=F,
                                             hllink = F, hltype = F, hlstart = start_x - start_x_pad, hlend = end_x - start_x_pad)
  
  # Save alignments
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
                                      gridplot_save = outfile_plot_grid, pregridplot_save = outfile_plot_pre_grid,
                                      max_n_alns = 150)
    gridmatrix = gridlist_to_gridmatrix(grid_xy)
    browser()
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
  
  if (clean_all){
    if (!is.na(file.size(genome_x_fa_subseq))){
        system(paste0('rm ', genome_x_fa_subseq))
    }
    if (!is.na(file.size(genome_y_fa_subseq))){
        system(paste0('rm ', genome_y_fa_subseq))
    }
    if (!is.na(file.size(paste0(genome_x_fa_subseq,'.chunk.fa')))){
        system(paste0('rm ', genome_x_fa_subseq,'.chunk.fa'))
    }
    if (!is.na(file.size(paste0(genome_x_fa_subseq, '.chunk.fa')))){
        system(paste0('rm ', genome_y_fa_subseq, '.chunk.fa'))
    }
  }
  save_to_logfile(get('log_collection', envir=globalenv()), res, logfile)

  
  
}



