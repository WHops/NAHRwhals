



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
wrapper_aln_and_analyse <- function(params) {
  
  # Directly enter debug mode? 
  if (params$debug) {
    print('Debug mode!')
    browser()
  }
  
  # Work out some reference genome things...
  if (params$reference_genome == 'hg38'){
    params$alt_ref_sample = F
  } else if (params$reference_genome == 'T2T'){
    print('T2T chosen. Translating reference coodinates...')
    params$alt_ref_sample = 'T2T'
  } else {
    print('Error: reference_genome must be "hg38" or "T2T".')
    return()
  }

  # Start writing a log file. 
  log_collection <<- init_log_with_def_values()
  log_collection[c('chr',
                   'start',
                   'end',
                   'xpad',
                   'chunklen',
                   'samplename',
                   'depth')] <<-
    c(params$seqname_x, params$start_x, params$end_x, params$xpad, 
      params$chunklen, params$samplename, params$depth)
  
  # Determine 'main' output name for this run
  sequence_name_output = manufacture_output_res_name(
    params$seqname_x, params$start_x, params$end_x
  )
  # Create output folder tree
  make_output_folder_structure(sequence_name_output)
  # Define output files
  outlinks = define_output_files(sequence_name_output, params$samplename)
  if (params$use_paf_library){
      
    # Step 1: Get the sequences (write to disk)
    chr_start_end_pad_params = extract_sequence_wrapper(params, outlinks)
    chr_start_end_pad = chr_start_end_pad_params[[1]]
    params = chr_start_end_pad_params[[2]]
  } else {
    chr_start_end_pad = list(chr='seqname', start='1', end=as.numeric(nchar(read.table(params$genome_x_fa))))
    system(paste0('cp ', params$genome_x_fa, ' ', outlinks$genome_x_fa_subseq))
    system(paste0('cp ', params$genome_y_fa, ' ', outlinks$genome_y_fa_subseq))

  }
  # Step 2: Run the alignments
  plot_x_y = produce_pairwise_alignments_minimap2(params, outlinks, chr_start_end_pad)

  # Step 2.1: If plot_only, exit. No SV calls. 
  if (params$plot_only){
    print('Plots done. Not attemptying to produce SV calls (plot_only is set to TRUE).')
    if (params$clean_after_yourself){
      clean_after_yourself(outlinks)
    }
    return()
  }
  
  # Step 2.2: If alignment has a contig break, exit. No SV calls.
  if (log_collection$exceeds_y){
    print('Assembly contig is broken. Not attempting to produce SV calls')
    res_empty = data.frame(eval=0, mut1='ref')
    res_empty = annotate_res_out_with_positions_lens(res_empty, NULL)
    
    write_results(res_empty, outlinks, params)
    return()
  }
  
  # Step 3: Condense and make a condensed plot
  grid_xy = wrapper_condense_paf(params, outlinks)
  
  # Step 3.1: If the alignment is cluttered, exit. No SV calls.
  if (is.null(grid_xy)){
    print('Dotplot is cluttered. SV calculation is not proceeded.')
    res_empty = data.frame(eval=0, mut1='ref')
    res_empty = annotate_res_out_with_positions_lens(res_empty, NULL)
    write_results(res_empty, outlinks, params)
    
    return()
  } else {
    # What do we need here? 
  }
  
  
  make_segmented_pairwise_plot(grid_xy, plot_x_y, outlinks)
  gridmatrix = gridlist_to_gridmatrix(grid_xy)
  #saveRDS(gridmatrix, file='~/Desktop/sec_advanced')
  # Step 4: Solve and make a solved plot
  #res = solve_mutation_old(gridmatrix, depth = params$depth, discovery_exact = params$discovery_exact)
  res = solve_mutation(gridmatrix, maxdepth = params$depth)#, discovery_exact = params$discovery_exact)

  make_modified_grid_plot(res, gridmatrix, outlinks)
  # Step 5: Save
  write_results(res, outlinks, params)
  

        
}




#' TODO: describe
#' @export
determine_xpad <- function(start, end) {
  # Play with padding values
  if ((end - start) < 1000) {
    xpad = 3
  } else if ((end - start) < 100000) {
    xpad = 3
  } else if ((end - start) > 5000000){
    xpad = 1
  } else {
    xpad = 2
  }
  return(xpad)
}

#' TODO: describe
#' @export
determine_plot_minlen <- function(start, end){
  if ((end - start) > 100000){
    minlen = 500
  }
  if ((end - start) > 5000000){
    minlen = 1000
  } else if ((end - start) > 10000000){
    minlen = 2000
  } else if ((end-start) < 5000){
    minlen = 100
  }else {
    minlen = 200
  }
  return(minlen)
}


#' TODO: describe
#' @export
determine_chunklen <- function(start, end) {
  
  # This is a bit poorly encoded. 
  length_upper_cutoffs = c(10000, 50000, 1000000, 2000000, 5000000, 1e10)
  chunklens =            c(100, 1000,  10000,   10000,   10000,   50000)
  
  interval_len = end - start
  
  dists = length_upper_cutoffs - interval_len
  dists[dists<0] = Inf
  
  chunklen = chunklens[which.min(dists[dists>0])]
  
  return(chunklen)
}


#' TODO: describe
#' @export
determine_compression <- function(start, end) {
  
  
  # if ((end - start) > 5000 * 1000) {
  #   compression = 100000
  # }
  # if ((end - start) > 500 * 1000) {
  #   compression = 10000
  # } else if ((end-start) < 5000){
  #   compression = 500
  # } else if (((end - start)) < 50 * 1000) {
  #   compression = 1000
  # } else {
  #   compression = 1000
  # }
  
  
  # Compression is chosen to be 500 bp for short stretches (< 5000 bp),
  # 1 kbp for medium sized (10 kbp - 100 kbp), 10 kbp for long (100 kbp - 1 Mbp) 
  # and 20 kbp for very long (> 1)
  
  
  # This is a bit poorly encoded. 
  # length_upper_cutoffs = c(5000, 500000, 1000000, 2000000, 5000000, 1e10)
  # compressions =         c(500,  1000,   10000,   20000,   50000,   100000)
  # 
  length_upper_cutoffs = c(50000, 500000, 5000000, 1e10)
  compressions =         c(100,  1000,   10000,   20000)
  
  interval_len = end - start
  
  dists = length_upper_cutoffs - interval_len
  dists[dists<0] = Inf
  
  compression = compressions[which.min(dists[dists>0])]
  
  return(compression)
}


#' TODO: describe
#' @export
save_plot_custom <-
  function(inplot,
           filename,
           device,
           width = 20,
           height = 20,
           units = 'cm') {
    ggplot2::ggsave(
      filename = paste0(filename, '.', device),
      plot = inplot,
      width = width,
      height = height,
      units = 'cm',
      dpi = 300,
      device = device
    )
    print('plot saved.')
  }



#' manufacture_output_res_name
#' @export
manufacture_output_res_name <- function(seqname_x, start_x, end_x){
  
  # Manufacture the name
  sequence_name_output = paste(
    paste0('res/', seqname_x),
    format(start_x, scientific = F),
    format(end_x, scientific = F),
    sep = '-'
  )
  
  
  return(sequence_name_output)
}

#' manufacture_output_res_name
#' @export
make_output_folder_structure <- function(sequence_name_output){
  
  # Create a lot of subfolders
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
  
}

#' Define a hell lot of output files. 
#' There must be a more elegant way to do this btw 
#' @export
define_output_files <- function(sequence_name_output, samplename){
  
  outlinks = list()
  outlinks$outpaf_link_self_x =  paste0(sequence_name_output,
                                        '/self/paf/aln_ref',
                                        samplename,
                                        '.paf')
  outlinks$outpaf_link_self_y =  paste0(sequence_name_output, '/self/paf/', samplename, '_y.paf')
  outlinks$outpaf_link_x_y =     paste0(sequence_name_output, '/diff/paf/', samplename, '_xy.paf')
  
  outlinks$res_table_xy =        paste0(sequence_name_output, '/diff/', samplename, '_res.tsv')
  
  outlinks$outfile_plot_self_x = paste0(sequence_name_output, '/self/pdf/', samplename, '_x_self')
  outlinks$outfile_plot_self_y = paste0(sequence_name_output, '/self/pdf/', samplename, '_y')
  outlinks$outfile_plot_x_y =    paste0(sequence_name_output, '/diff/pdf/', samplename, '_x_y')
  
  outlinks$outfile_plot_pre_grid = paste0(sequence_name_output,
                                          '/diff/pdf/grid/',
                                          samplename,
                                          '_x_y_grid_pre.pdf')
  outlinks$outfile_plot_grid =     paste0(sequence_name_output,
                                          '/diff/pdf/grid/',
                                          samplename,
                                          '_x_y_grid.pdf')
  outlinks$outfile_colored_segment =     paste0(sequence_name_output,
                                                '/diff/pdf/grid/',
                                                samplename,
                                                '_x_y_colored.pdf')
  outlinks$outfile_plot_grid_mut = paste0(sequence_name_output,
                                          '/diff/pdf/grid/',
                                          samplename,
                                          '_x_y_grid_mut.pdf')
  
  outlinks$genome_x_fa_subseq = paste0(sequence_name_output, '/fasta/', samplename, '_x.fa')
  outlinks$genome_y_fa_subseq = paste0(sequence_name_output, '/fasta/', samplename, '_y.fa')
  
  outlinks$genome_x_fa_altref_subseq = paste0(sequence_name_output, '/fasta/', samplename, 'altref_x.fa')
  outlinks$genome_y_fa_altref_subseq = paste0(sequence_name_output, '/fasta/', samplename, 'altref_y.fa')
  
  return(outlinks)
}

#' @export
clean_after_yourself <- function(outlinks){
  if (!is.na(file.size(outlinks$genome_x_fa_subseq))) {
    system(paste0('rm ', outlinks$genome_x_fa_subseq))
  }
  if (!is.na(file.size(outlinks$genome_y_fa_subseq))) {
    system(paste0('rm ', outlinks$genome_y_fa_subseq))
  }
  if (!is.na(file.size(paste0(outlinks$genome_x_fa_subseq, '.chunk.fa')))) {
    system(paste0('rm ', outlinks$genome_x_fa_subseq, '.chunk.fa'))
  }
  if (!is.na(file.size(paste0(outlinks$genome_x_fa_subseq, '.chunk.fa')))) {
    system(paste0('rm ', outlinks$genome_y_fa_subseq, '.chunk.fa'))
  }
}

#' @export
write_results <- function(res, outlinks, params){
  # Save res table
  write.table(
    res,
    file = outlinks$res_table_xy,
    col.names = T,
    row.names = F,
    quote = F,
    sep = '\t'
  )
  
  # Save to logfile
  save_to_logfile(get('log_collection', envir = globalenv()), res, params$logfile, params, alt_x = (params$alt_ref_sample != F))
  
}

#' @export 
wrapper_condense_paf <- function(params, outlinks){
  
  # Make the condensation
  grid_xy = wrapper_paf_to_bitlocus(
    outlinks$outpaf_link_x_y,
    params,
    gridplot_save = outlinks$outfile_plot_grid,
    pregridplot_save = outlinks$outfile_plot_pre_grid
  )
  
  return(grid_xy)
}


make_segmented_pairwise_plot <- function(grid_xy, plot_x_y, outlinks){
  # Make plot_xy_segmented. 
  # Needs debug. 
  xstart = (grid_xy[[1]][1:length(grid_xy[[1]])-1])
  xend = (grid_xy[[1]][2:length(grid_xy[[1]])])
  ystart = (grid_xy[[2]][1:length(grid_xy[[2]])-1])
  yend = (grid_xy[[2]][2:length(grid_xy[[2]])])
  xmax = max(grid_xy[[1]])
  ymax = max(grid_xy[[2]])
  datx = data.frame(xstart = xstart, 
                    xend = xend,
                    xmax = xmax
  )
  daty = data.frame(yend = yend,
                    ymax = ymax,
                    ystart = ystart
  )
  
  likely_stepsize = min(c(diff(datx$xstart), diff(daty$ystart)))
  
  # Introducing special case for nrow(datx) = 1
  if (likely_stepsize == Inf){
    likely_stepsize = 0
  }
  
  if (length(xstart) > 433){
    print('Too many segments to make colored plot')
    return()
  }
  
  qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector = rep(unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))), 10)
  
  plot_x_y_segmented = plot_x_y + 
    ggplot2::geom_rect(data=datx,
                       ggplot2::aes(xmin=xstart, xmax=xend, ymin=0, ymax=ymax, fill=col_vector[1:length(xstart)]),
                       alpha=0.5) + 
    ggplot2::guides(fill = FALSE) +
    ggplot2::xlim(c(0,max(ggplot2::layer_scales(plot_x_y)$x$range$range, xmax+likely_stepsize))) + # Range is the max of previous plot and new additions. So that nothing gets cut off. 
    ggplot2::ylim(c(0,max(ggplot2::layer_scales(plot_x_y)$y$range$range, ymax+likely_stepsize))) +
    ggplot2::scale_x_continuous(labels = scales::comma) +
    ggplot2::scale_y_continuous(labels = scales::comma) 
    # ggplot2::geom_segment(data=daty,
    #             ggplot2::aes(x=0, xend=xmax, y=ystart, yend=ystart), color='grey')
  print(plot_x_y_segmented)
  
  ggplot2::ggsave(filename = paste0(outlinks$outfile_colored_segment),
                  plot = plot_x_y_segmented,
                  width = 15,
                  height = 15,
                  units = 'cm',
                  dpi = 300)
  
}

#' @export
make_modified_grid_plot <- function(res, gridmatrix, outlinks){
  
  res = res[order(res$eval, decreasing = T),]
  # res = filter_res(res, threshold = params$eval_th)
  
  grid_modified = modify_gridmatrix(gridmatrix, res[1,])
  
  
  gridlines.x = cumsum(c(0,as.numeric(colnames(grid_modified))))
  gridlines.y = cumsum(c(0,as.numeric(row.names(grid_modified))))
  
  colnames(grid_modified) = seq(1:dim(grid_modified)[2])
  row.names(grid_modified) = seq(1:dim(grid_modified)[1])
  
  gm2 = reshape2::melt(grid_modified)
  colnames(gm2) = c('y','x','z')
  grid_mut_plot = plot_matrix_ggplot_named(gm2[gm2$z != 0,], gridlines.x, gridlines.y)
  ggplot2::ggsave(filename = paste0(outlinks$outfile_plot_grid_mut),
                  plot = grid_mut_plot,
                  width = 10,
                  height = 10,
                  units = 'cm',
                  dpi = 300)
  print(grid_mut_plot)
  
}
