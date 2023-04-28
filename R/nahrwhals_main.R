
#' Nahrwhals: A function for analyzing genomic alignments and visualizing results
#'
#' This function performs sequence alignment between two genomic sequences (fasta files) and
#' generates visualizations of the alignment results.
#'

#' @param genome_x_fa Path to reference fasta file (default: "testdata/assemblies/hg38_partial.fa").
#' @param genome_y_fa Path to query fasta file (default: "testdata/assemblies/assembly_partial.fa").
#' @param seqname_x Sequence name for the reference fasta file (default: "chr1_partial").
#' @param start_x Start position for the reference sequence (default: 1700000).
#' @param end_x End position for the reference sequence (default: 3300000).
#' @param genome_y_fa_mmi Path to the .mmi index of the query fasta file or "default" (default: "default").
#' @param anntrack Logical, whether to use annotation track (default: FALSE).
#' @param logfile Path to output log file (default: "res/unittest.tsv").
#' @param samplename_x Sample name for the reference fasta file (default: "Fasta_x").
#' @param samplename_y Sample name for the query fasta file (default: "Fasta_y").
#' @param compare_full_fastas Logical, whether to compare full fasta files (default: FALSE).
#' @param plot_only Logical, whether to perform plotting only (default: FALSE).
#' @param self_plots Logical, whether to generate self plots (default: TRUE).
#' @param plot_xy_segmented Logical, whether to generate segmented XY plots (default: TRUE).
#' @param eval_th Evaluation threshold (default: 98).
#' @param depth Depth of the analysis (default: 3).
#' @param chunklen Chunk length or "default" (default: "default").
#' @param minlen Minimum length or "default" (default: "default").
#' @param compression Compression or "default" (default: "default").
#' @param max_size_col_plus_rows Maximum sum of rows and columns for alignment plots (default: 250).
#' @param max_n_alns Maximum number of alignments to process (default: 150).
#' @param use_paf_library Logical, whether to use the PAF library (default: FALSE).
#' @param conversionpaf_link Logical, whether to use the ConversionPAF link (default: FALSE).
#' @param xpad Padding for the X-axis (default: 1).
#' @param plot_minlen Minimum length for plotting (default: 350).
#' @param maxlen_refine Maximum length to refine (default: 1e10).
#' @param n_tests Number of tests to perform (default: 10).
#' @param n_max_testchunks Maximum number of test chunks (default: 5).
#' @param baseline_log_minsize_min Minimum log size for the baseline (default: 8).
#' @param baseline_log_minsize_max Maximum log size for the baseline or "default" (default: FALSE).
#' @param discovery_exact Logical, whether to use exact discovery mode (default: FALSE).
#' @param hltrack Logical, whether to use the high-level track (default: FALSE).
#' @param hllink Logical, whether to use the high-level link (default: FALSE).
#' @param aln_pad_factor Alignment padding factor (default: 1.0).
#' @param debug Logical, whether to enable debug mode (default: FALSE).
#' @param clean_after_yourself Logical, whether to delete all non-plot files, such as fastas (default: FALSE)
#' @param testrun_std Logical, whether to perform a test/demo run. Mut. Exclusive with testrun_fullfa. (default: FALSE)
#' @param testrun_fullfa Logical, whether to perform a test/demo run of full fastas. Mut. Exclusive with testrun_std (default: FALSE). 
#' @param minimap2_bin Path to minimap2 binary or "default" (default: "default").
#' @param bedtools_bin Path to bedtools binary or "default" (default: "default").
#' @return No explicit return value; function generates output files and visualizations.
#'
#'
#' @export
nahrwhals <- function(
                      genome_x_fa, genome_y_fa, seqname_x, start_x, end_x,
                      genome_y_fa_mmi = "default", anntrack = FALSE,
                      logfile = 'res/res.tsv', samplename_x = 'Fasta_x',
                      samplename_y = 'Fasta_y', compare_full_fastas = FALSE,
                      plot_only = FALSE, self_plots = TRUE,
                      plot_xy_segmented = TRUE, eval_th = 98, depth = 3,
                      chunklen = "default", minlen = "default", compression = "default",
                      max_size_col_plus_rows = 250, max_n_alns = 150,
                      use_paf_library = FALSE, conversionpaf_link = FALSE,
                      xpad = 1, plot_minlen = 350, maxlen_refine = 1e10,
                      n_tests = 10, n_max_testchunks = 5, baseline_log_minsize_min = 8,
                      baseline_log_minsize_max = FALSE, discovery_exact = FALSE,
                      hltrack = FALSE, hllink = FALSE, aln_pad_factor = 1.0,
                      debug = FALSE, clean_after_yourself = FALSE, testrun_std = FALSE,
                      testrun_fullfa = FALSE, minimap2_bin = "default", bedtools_bin = "default"
                      ) {



  # Get the function call
  call <- match.call()
  
  # Get the formal arguments with their default values
  defaults <- formals(nahrwhals)
  
  # Remove the first element (function name)
  call <- as.list(call)[-1]
  
  # Evaluate the call in the parent frame
  call_evaluated <- do.call(function(...) list(...), 
                            eval(call, envir = parent.frame()))
  
  # Merge the provided parameters and default values
  params <- modifyList(defaults, call_evaluated)
  
  if (testrun_std && testrun_fullfa) {
    stop("Parameters testrun_std and testrun_fullfa cannot both be TRUE.")
  }

  if (testrun_std){
    cat('Testmode! Running a standard testrun with sample data.')
    params$genome_x_fa = system.file("extdata/assemblies", "hg38_partial.fa", package="nahrwhals")
    params$genome_y_fa = system.file("extdata/assemblies", "assembly_partial.fa", package="nahrwhals")
    params$seqname_x = 'chr1_partial'
    params$start_x = 1700000
    params$end_x = 3300000
  } else if (testrun_fullfa){
    cat('Testmode! Running a "full fasta to full fasta" testrun with sample data.')
    params$genome_x_fa = system.file("extdata/extracted_fastas", "sequence1.fa", package="nahrwhals")
    params$genome_y_fa = system.file("extdata/extracted_fastas", "sequence2.fa", package="nahrwhals")
    params$compare_full_fastas = T
  }

  ###### Compute parameters that are functions of other parameters. #####
  ###### (unless those parameters have been specified explicitly; ) #####
  ###### i.e. they are not set to 'default' or similar.             #####
  default_param_values = c('default', 'Default', 'auto', 'Auto', '', NA, NULL, F)

  if (params$chunklen %in% default_param_values){
    params$chunklen = determine_chunklen(params$start_x, params$end_x)
  }
  if (params$minlen %in% default_param_values){
    params$minlen = determine_compression(params$start_x, params$end_x)
  }
  if (params$compression %in% default_param_values){
    params$compression = determine_compression(params$start_x, params$end_x)
  }
  if (params$genome_y_fa_mmi %in% default_param_values){
    params$genome_y_fa_mmi = paste0(params$genome_y_fa, '.mmi')
  }
  if (params$bedtools_bin %in% default_param_values){
    params$bedtools_bin = 'bedtools'
  }
  if (params$minimap2_bin %in% default_param_values){
    params$minimap2_bin = 'minimap2'
  }
  if (params$baseline_log_minsize_max %in% default_param_values){
    params$baseline_log_minsize_max = max(log2(20000), log2((params$end_x-params$start_x) / 10))
  }

  # For better reporting of results in fasta_direct
  if (params$compare_full_fastas == T){
    params$seqname_x = 'manual'
    params$start_x = 1
    params$end_x = Inf
  }
  # Print final parameter values 
  # for (param_name in names(params)) {
  #   cat(param_name, "=", params[[param_name]], "\n")
  # }

  ##### Run the main NAHRwhals wrapper function #####
  wrapper_aln_and_analyse(params)
    
  print('done!')

}