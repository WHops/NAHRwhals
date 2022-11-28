



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
    
  
  # Step 1: Get the sequences (write to disk)
  start_end_pad = extract_sequence_wrapper(params, outlinks)
  
  # Step 2: Run the alignments
  plot_x_y = produce_pairwise_alignments_minimap2(params, outlinks, start_end_pad[1], start_end_pad[2])

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
    res = data.frame(eval = 0, mut1 = 'ref', mut2 = NA, mut3 = NA)
    write_results(res, outlinks, params)
    return()
  }
  
  # Step 3: Condense and make a condensed plot
  grid_xy = wrapper_condense_paf(params, outlinks)
  make_segmented_pairwise_plot(grid_xy, plot_x_y, outlinks)
  gridmatrix = gridlist_to_gridmatrix(grid_xy)
  #browser()
  # Step 4: Solve and make a solved plot
  #res = solve_mutation_old(gridmatrix, depth = params$depth, discovery_exact = params$discovery_exact)
  res = solve_mutation(gridmatrix, maxdepth = params$depth)#, discovery_exact = params$discovery_exact)
  make_modified_grid_plot(res, gridmatrix, outlinks)
    
  # Step 5: Save
  write_results(res, outlinks, params)
  

        
}

