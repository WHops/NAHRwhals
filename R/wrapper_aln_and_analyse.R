#' wrapper_aln_and_analyse
#' @description The main wrapper function for running NAHRwhals.
#' This is where all the information flow together.
#' @param params a single list containing all input information.
#' @export
wrapper_aln_and_analyse <- function(params) {
  # Directly enter debug mode?
  if (params$debug) {
    print("Debug mode!")
  }
  # Start writing a log file.
  log_collection <<- init_log_with_def_values()
  log_collection[c(
    "chr",
    "start",
    "end",
    "xpad",
    "chunklen",
    "samplename",
    "depth"
  )] <<-
    c(
      params$seqname_x, params$start_x, params$end_x, params$xpad,
      params$chunklen, paste0(params$samplename_x, "_", params$samplename_y), params$depth
    )
  
  # Determine 'main' output name for this run
  sequence_name_output <- manufacture_output_res_name(
    params$outdir, params$seqname_x, params$start_x, params$end_x
  )
  
  # Create output folder tree
  make_output_folder_structure(params$outdir, sequence_name_output)
    
  # Define output files
  outlinks <- define_output_files(sequence_name_output, paste0(params$samplename_x, "_", params$samplename_y))
  while(!is.null(dev.list())){dev.off()}

  # start a pdf
  pdf(file = outlinks$outpdf_main)

  

  # Print only if params$silent is FALSE
  if (!params$silent) {
    print("Step 1: Identifying homologous regions in the y assembly")
  }

  if (params$compare_full_fastas == T) {
    print("Option pairwise_fasta_direct recognized. Comparing entire fasta files.")
    chr_start_end_pad <- list(chr = "seqname", start = "1", end = as.numeric(nchar(read.table(params$genome_x_fa))))
    run_silent(paste0("cp ", params$genome_x_fa, " ", outlinks$genome_x_fa_subseq))
    run_silent(paste0("cp ", params$genome_y_fa, " ", outlinks$genome_y_fa_subseq))


  } else if (params$use_paf_library == T) {
    chr_start_end_pad_params <- extract_sequence_wrapper(params, outlinks)
    chr_start_end_pad <- chr_start_end_pad_params[[1]]
    params <- chr_start_end_pad_params[[2]]

  } else if (params$use_paf_library == F) {
    params <- make_params_conversionpaf(params, outlinks)
    if (is.na(params$conversionpaf_link)){
      res_empty <- data.frame(eval = 0, mut1 = "ref")
      res_empty <- annotate_res_out_with_positions_lens(res_empty, NULL)
      write_results(res_empty, outlinks, params)
      while(!is.null(dev.list())){dev.off()}

      return()
    }
    chr_start_end_pad_params <- extract_sequence_wrapper(params, outlinks)
    
    if (is.null(chr_start_end_pad_params)){
      res_empty <- data.frame(eval = 0, mut1 = "ref")
      res_empty <- annotate_res_out_with_positions_lens(res_empty, NULL)
      write_results(res_empty, outlinks, params)
      while(!is.null(dev.list())){dev.off()}

      return()
    }
    chr_start_end_pad <- chr_start_end_pad_params[[1]]
    params <- chr_start_end_pad_params[[2]]
    run_silent(paste0("rm ", params$conversionpaf_link))
  }

  if (log_collection$cluttered_boundaries == T){
    write_results(data.frame(), outlinks, params)
    log_collection <<- init_log_with_def_values()
    while(!is.null(dev.list())){dev.off()}

    return('EARLY FINISH!')
  }

  if (!params$silent) {
  print("Step 2: Computing detailed pairwise alignment")
  }
  plot_x_y <- produce_pairwise_alignments_minimap2(params, outlinks, chr_start_end_pad)


  # Step 2.1: If plot_only, exit. No SV call
  if (params$plot_only) {
    print("Plots done. Not attemptying to produce SV calls (plot_only is set to TRUE).")
    if (params$clean_after_yourself) {
      clean_after_yourself(outlinks)
    }

    while(!is.null(dev.list())){dev.off()}
    return()
  }

  if (exists("log_collection")) {
    log_collection$exceeds_y <<- F
  }
  # Step 2.2: If alignment has a contig break, exit. No SV calls.
  if (log_collection$exceeds_y) {
    print("Assembly contig is broken. Not attempting to produce SV calls")
    res_empty <- data.frame(eval = 0, mut1 = "ref")
    res_empty <- annotate_res_out_with_positions_lens(res_empty, NULL)

    write_results(res_empty, outlinks, params)
    while(!is.null(dev.list())){dev.off()}

    return()
  }

  if (!params$silent) {
  print('Step 3: Segmenting pairwise alignment')
  }

  grid_xy <- wrapper_condense_paf(params, outlinks)
  # Step 3.1: If the alignment is cluttered, exit. No SV calls.
  if (is.null(grid_xy)) {
    print("Dotplot is cluttered. SV calculation is not proceeded.")
    res_empty <- data.frame(eval = 0, mut1 = "ref")
    res_empty <- annotate_res_out_with_positions_lens(res_empty, NULL)
    write_results(res_empty, outlinks, params)
    while(!is.null(dev.list())){dev.off()}

    return()
  } else {
    # What do we need here?
  }

  make_segmented_pairwise_plot(grid_xy, plot_x_y, outlinks)

  # Make in between another plot

  suppressMessages(suppressWarnings({
  p <- plot_matrix_ggplot_named(grid_xy[[3]], grid_xy[[1]], grid_xy[[2]])
  print(p)
  ggplot2::ggsave(
    p,
    file = outlinks$outfile_plot_grid,
    height = 15,
    width = 15,
    units = "cm",
    device = "pdf"
  )
  }))
  gridmatrix <- gridlist_to_gridmatrix(grid_xy)

  
  # reject bitloci with cut-off alignments at the borders ('clutter')
  if (is_cluttered(gridmatrix)){#} | log_collection$cluttered_boundaries==T) {
    print("Alignments overlap with borders. Make frame larger!")
    res_empty <- data.frame(eval = 0, mut1 = "ref")
    res_empty <- annotate_res_out_with_positions_lens(res_empty, NULL)
    write_results(res_empty, outlinks, params)
    while(!is.null(dev.list())){dev.off()}

    return()
  }
  if (!params$silent) {
  print('Step 4: BFS search for mutation chains.')
  }
  res <- solve_mutation(gridmatrix, maxdepth = 1, solve_th = params$eval_th, compression = params$compression, is_cluttered_already_paf = log_collection$cluttered_boundaries==T)
  if (max(res$eval) < 99.8){
    # Fahre schwere Geschuetze auf
    res_julia <- solve_mutation_julia_wrapper(params, gridmatrix, grid_xy[[1]], outlinks$solver_inmat_path, outlinks$solver_inlens_path, outlinks$solver_juliares_path, params$julia_solver_path)
    
    if (!all((dim(res_julia)) == c(1,3))){
      if (!'ref' %in% res_julia$mut1){
        res_julia = rbind(res_julia, c(res[res$mut1 == 'ref', 'eval'], 'ref', rep(NA,ncol(res_julia)-2)))
        res_julia$eval <- as.numeric(as.character(res_julia$eval))

        # Also wir nehmen die gewinner, und dann den gewinner mit der kleinsten depth und dann der höchsten eval? 
        res <- res_julia[order(-res_julia$eval, -rowSums(is.na(res_julia))), ]
      } else  {
        res = res_julia
      }
    }
  }



  make_modified_grid_plot(res, gridmatrix, outlinks)

  Sys.sleep(1)  
  # Step 5: Save
  write_results(res, outlinks, params)
  while(!is.null(dev.list())){dev.off()}

}




#' Determine desired padding value for the x-axis of a dot plot
#'
#' This function determines the optimal padding value for the x-axis of a dot plot based on the length of the input sequence. It uses a set of pre-defined length cutoffs and padding values to choose the best padding value for the input sequence.
#'
#' @param start The starting coordinate of the input sequence.
#' @param end The ending coordinate of the input sequence.
#'
#' @return The desired padding value for the x-axis of the dot plot.
#'
#' @export
determine_xpad <- function(start, end) {
  # Play with padding values
  if ((end - start) < 1000) {
    xpad <- 3
  } else if ((end - start) < 100000) {
    xpad <- 3
  } else if ((end - start) > 5000000) {
    xpad <- 1
  } else {
    xpad <- 2
  }
  return(xpad)
}

#' Determine optimal minimum length for plotting a dot plot
#'
#' This function determines the optimal minimum length for plotting a dot plot based on the length of the input sequence. It uses a set of pre-defined length cutoffs and minimum length values to choose the best minimum length for the input sequence.
#'
#' @param start The starting coordinate of the input sequence.
#' @param end The ending coordinate of the input sequence.
#'
#' @return The optimal minimum length for plotting the dot plot.
#'
#' @export
determine_plot_minlen <- function(start, end) {
  if ((end - start) > 100000) {
    minlen <- 500
  }
  if ((end - start) > 5000000) {
    minlen <- 1000
  } else if ((end - start) > 10000000) {
    minlen <- 2000
  } else if ((end - start) < 5000) {
    minlen <- 100
  } else {
    minlen <- 200
  }
  return(minlen)
}


#' Determine optimal chunk length for processing a sequence
#'
#' This function determines the optimal chunk length for processing a sequence based on the length of the input sequence. It uses a set of pre-defined upper length cutoffs and chunk length values to choose the best chunk length for the input sequence.
#'
#' @param start The starting coordinate of the input sequence.
#' @param end The ending coordinate of the input sequence.
#'
#' @return The optimal chunk length for processing the input sequence.
#'
#' @export
determine_chunklen <- function(start, end) {
  # This is a bit poorly encoded.
  length_upper_cutoffs <- c(10000, 500000, 1000000, 2000000, 5000000, 1e10)
  chunklens <- c(100, 1000, 10000, 10000, 10000, 50000)

  interval_len <- end - start

  dists <- length_upper_cutoffs - interval_len
  dists[dists < 0] <- Inf

  chunklen <- chunklens[which.min(dists[dists > 0])]

  return(chunklen)
}


#' Determine optimal compression parameter for segmenting a dot plot
#'
#' This function determines the optimal compression parameter for segmenting a dot plot based on the length of the input sequence. It uses a set of pre-defined upper length cutoffs and compression values to choose the best compression parameter for the input sequence.
#'
#' @param start The starting coordinate of the input sequence.
#' @param end The ending coordinate of the input sequence.
#'
#' @return The optimal compression parameter for segmenting the dot plot.
#'
#' @export
determine_compression <- function(start, end) {
  length_upper_cutoffs <- c(50000, 500000, 5000000, 1e10)
  compressions <- c(100, 1000, 10000, 20000)

  interval_len <- end - start

  dists <- length_upper_cutoffs - interval_len
  dists[dists < 0] <- Inf

  compression_parameter_bp <- compressions[which.min(dists[dists > 0])]

  return(compression_parameter_bp)
}


#' Save ggplot plot to file using custom settings
#'
#' This function saves a ggplot plot to a file using custom settings for the filename, device, width, height, and units. It is a quick helper function for saving ggplot plots.
#'
#' @param inplot The ggplot plot object to be saved.
#' @param filename The base name for the output file.
#' @param device The device for output file generation (e.g. "pdf", "png", etc.).
#' @param width The width of the plot in the specified units.
#' @param height The height of the plot in the specified units.
#' @param units The units for the plot dimensions (e.g. "cm", "in", etc.).
#'
#' @export
save_plot_custom <-
  function(inplot,
           filename,
           device,
           width = 20,
           height = 20,
           units = "cm") {
    ggplot2::ggsave(
      filename = paste0(filename, ".", device),
      plot = inplot,
      width = width,
      height = height,
      units = "cm",
      dpi = 300,
      device = device
    )
    #print("plot saved.")
  }



#' Create output file name for NAHRwhals results
#'
#' This function creates an output file name for the results of a NAHRwhals analysis, based on the input sequence name, start and end coordinates.
#'
#' @param seqname_x The name of the input sequence.
#' @param start_x The starting coordinate of the input sequence.
#' @param end_x The ending coordinate of the input sequence.
#'
#' @return The output file name for the NAHRwhals results.
#'
#' @export
manufacture_output_res_name <- function(out_dir, seqname_x, start_x, end_x) {
  # Manufacture the name
  sequence_name_output <- paste(
    paste0(out_dir, "/", seqname_x),
    format(start_x, scientific = F),
    format(end_x, scientific = F),
    sep = "-"
  )


  return(sequence_name_output)
}

#' Create empty output folders for NAHRwhals output files
#'
#' This function creates empty output folders for various intermediate and final files generated by NAHRwhals.
#'
#' @param sequence_name_output The output directory path for the files.
#'
#' @export
make_output_folder_structure <- function(outdir, sequence_name_output) {
  # Create a lot of subfolders

  subfolders_to_create = c(
    "self/pdf",
    "self/paf",
    "diff/pdf/grid",
    "diff/paf",
    "fasta"
  )

  for (subfolder in subfolders_to_create) {
    dir.create(paste0(sequence_name_output, "/", subfolder), showWarnings = F, recursive = T)
  }


}

#' Define output file paths for NAHRwhals output
#'
#' This function defines output file paths for various intermediate and final files generated by NAHRwhals.
#'
#' @param sequence_name_output The output directory path for the files.
#' @param samplename The name of the sample being analyzed.
#'
#' @return A list of file paths for the output files.
#'
#' @export
define_output_files <- function(sequence_name_output, samplename) {
  outlinks <- list()

  outlinks$outpdf_main <- paste0(sequence_name_output, "/", samplename, "_all.pdf")
  outlinks$outpaf_link_self_x <- paste0(
    sequence_name_output,
    "/self/paf/aln_ref",
    samplename,
    ".paf"
  )
  outlinks$outpaf_link_self_y <- paste0(sequence_name_output, "/self/paf/", samplename, "_y.paf")
  outlinks$outpaf_link_x_y <- paste0(sequence_name_output, "/diff/paf/", samplename, "_xy.paf")

  outlinks$res_table_xy <- paste0(sequence_name_output, "/diff/", samplename, "_res.tsv")

  outlinks$outfile_plot_self_x <- paste0(sequence_name_output, "/self/pdf/", samplename, "_x_self")
  outlinks$outfile_plot_self_y <- paste0(sequence_name_output, "/self/pdf/", samplename, "_y")
  outlinks$outfile_plot_x_y <- paste0(sequence_name_output, "/diff/pdf/", samplename, "_x_y")

  outlinks$outfile_plot_pre_grid <- paste0(
    sequence_name_output,
    "/diff/pdf/grid/",
    samplename,
    "_x_y_grid_pre.pdf"
  )
  outlinks$outfile_plot_grid <- paste0(
    sequence_name_output,
    "/diff/pdf/grid/",
    samplename,
    "_x_y_grid.pdf"
  )
  outlinks$outfile_colored_segment <- paste0(
    sequence_name_output,
    "/diff/pdf/grid/",
    samplename,
    "_x_y_colored.pdf"
  )
  outlinks$outfile_plot_grid_mut <- paste0(
    sequence_name_output,
    "/diff/pdf/grid/",
    samplename,
    "_x_y_grid_mut.pdf"
  )

  outlinks$solver_inmat_path <- paste0(sequence_name_output, "/diff/", samplename, "_x_y_grid_mut.tsv")
  outlinks$solver_inlens_path <- paste0(sequence_name_output, "/diff/", samplename, "_x_y_grid_mut_lens.tsv")
  outlinks$solver_juliares_path <- paste0(sequence_name_output, "/diff/", samplename, "_x_y_grid_mut_juliares.tsv")

  outlinks$genome_x_fa_subseq <- paste0(sequence_name_output, "/fasta/", samplename, "_x.fa")
  outlinks$genome_y_fa_subseq <- paste0(sequence_name_output, "/fasta/", samplename, "_y.fa")

  outlinks$genome_x_fa_altref_subseq <- paste0(sequence_name_output, "/fasta/", samplename, "altref_x.fa")
  outlinks$genome_y_fa_altref_subseq <- paste0(sequence_name_output, "/fasta/", samplename, "altref_y.fa")

  return(outlinks)
}

#' Delete intermediate files created by NAHRwhals
#'
#' This function deletes intermediate fasta sequences created by NAHRwhals. It takes a list of file paths as input and deletes the files if they exist and are not empty.
#'
#' @param outlinks A list of file paths for intermediate fasta sequences.
#'
#' @export
clean_after_yourself <- function(outlinks) {
  if (!is.na(file.size(outlinks$genome_x_fa_subseq))) {
    run_silent(paste0("rm ", outlinks$genome_x_fa_subseq))
  }
  if (!is.na(file.size(outlinks$genome_y_fa_subseq))) {
    run_silent(paste0("rm ", outlinks$genome_y_fa_subseq))
  }
  if (!is.na(file.size(paste0(outlinks$genome_x_fa_subseq, ".chunk.fa")))) {
    run_silent(paste0("rm ", outlinks$genome_x_fa_subseq, ".chunk.fa"))
  }
  if (!is.na(file.size(paste0(outlinks$genome_x_fa_subseq, ".chunk.fa")))) {
    run_silent(paste0("rm ", outlinks$genome_y_fa_subseq, ".chunk.fa"))
  }
}

#' Write tree search results to output file and logfile
#'
#' This function writes the results of a tree search to an output file and saves a log of the results to a logfile.
#'
#' @param res A table of search results to write to an output file.
#' @param outlinks A list of output links where the resulting files will be saved.
#' @param params A list of input parameters for the search.
#'
#' @export
write_results <- function(res, outlinks, params) {
  # Save res table
  write.table(
    res,
    file = outlinks$res_table_xy,
    col.names = T,
    row.names = F,
    quote = F,
    sep = "\t"
  )

  # Save to logfile
  save_to_logfile(get("log_collection", envir = globalenv()), res, params$logfile, params, alt_x = F)
}

#' Wrapper function for condensing PAF files into bit-locus format
#'
#' This function serves as a wrapper for the `wrapper_paf_to_bitlocus` function. It condenses PAF files into bit-locus format based on the input parameters and saves the resulting plot files to disk.
#'
#' @param params A list of input parameters used for the bit-locus plot generation.
#' @param outlinks A list of output links where the resulting plot files will be saved.
#'
#' @return A matrix containing the condensed bit-locus data.
#'
#' @export
wrapper_condense_paf <- function(params, outlinks) {
  # Make the condensation
  grid_xy <- wrapper_paf_to_bitlocus(
    outlinks$outpaf_link_x_y,
    params,
    gridplot_save = outlinks$outfile_plot_grid,
    pregridplot_save = outlinks$outfile_plot_pre_grid
  )

  # p = plot_matrix_ggplot_named(grid_list, gridlines.x, gridlines.y)
  # print(p)
  return(grid_xy)
}
