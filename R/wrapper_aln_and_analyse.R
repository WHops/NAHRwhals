



#' wrapper_aln_and_analyse
#' @description The main wrapper function for running NAHRwhals.
#' This is where all the information flow together.
#' @param params a single list containing all input information.
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
  if (params$use_paf_library == T){
      
    # Step 1: Get the sequences (write to disk)
    chr_start_end_pad_params = extract_sequence_wrapper(params, outlinks)
    chr_start_end_pad = chr_start_end_pad_params[[1]]
    params = chr_start_end_pad_params[[2]]
    
  } else if (params$use_paf_library == 'auto'){
    
    create_mmi_if_doesnt_exists(params)
    # params$conversionpaf_link ; this is what we need. 
    params = make_params_conversionpaf(params, outlinks)
    chr_start_end_pad_params = extract_sequence_wrapper(params, outlinks)
    chr_start_end_pad = chr_start_end_pad_params[[1]]
    params = chr_start_end_pad_params[[2]]
    system(paste0('rm ', params$conversionpaf_link))
    
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


#' determine_compression
#' 
#' Compression parameter for segmenting the dotplot. 
#' We choose this based on the length of the sequence. 
#' @export
determine_compression <- function(start, end) {
  
  length_upper_cutoffs = c(50000, 500000, 5000000, 1e10)
  compressions =         c(100,  1000,   10000,   20000)
  
  interval_len = end - start
  
  dists = length_upper_cutoffs - interval_len
  dists[dists<0] = Inf
  
  compression_parameter_bp = compressions[which.min(dists[dists>0])]
  
  return(compression_parameter_bp)
}


#' save_plot_custom
#' 
#' Quick helperfunction for saving a ggplot plot. 
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
#' 
#' Creates an output file name that will be used by other functions.
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

#' make_output_folder_structure
#' 
#' Create empty output folders that can be filled with data later. 
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

#' clean_after_yourself
#' 
#' This function, if invoked, deletes intermediate fasta sequences created by NAHRwhals. 
#' By default, this is not invoked. 
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

#' write_results
#' 
#' Write results of a tree search to an output file. Uses the global 'log' variable.
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

#' wrapper_condense_paf
#' 
#' Wrapper around a wrapper ... consider to simplify / remove. 
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