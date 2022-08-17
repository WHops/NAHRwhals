



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
    system(paste0('cp ', genome_x_fa, ' ', outlinks$genome_x_fa_subseq))
    system(paste0('cp ', genome_y_fa, ' ', outlinks$genome_y_fa_subseq))
    
    start_x_pad = params$start_x
    end_x_pad = params$end_x

    system(paste0('cp ', genome_x_fa, ' ', outlinks$genome_x_fa_subseq))
    system(paste0('cp ', genome_y_fa, ' ', outlinks$genome_y_fa_subseq))
    
  }
  
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




  if (!params$plot_only) {
    
    # If alignment has a contig break, we don't have to find SVs.
    if (log_collection$exceeds_y == F){
      # Make an xy grid
      #params$compression= 50000
      #params$minlen = 50000
      grid_xy = wrapper_paf_to_bitlocus(
        outlinks$outpaf_link_x_y,
        params,
        gridplot_save = outlinks$outfile_plot_grid,
        pregridplot_save = outlinks$outfile_plot_pre_grid
      )
      gridmatrix = gridlist_to_gridmatrix(grid_xy)
      res = solve_mutation(gridmatrix, depth = params$depth, discovery_exact = params$discovery_exact)
      res$eval = as.numeric(res$eval)
      

      
      res = res[order(res$eval, decreasing = T),]
      res = filter_res(res, threshold = params$eval_th)
      print(head(res))
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
          plot(grid_mut_plot)

      
    } else if (log_collection$exceeds_y == T){
      res = data.frame(eval = 0, mut1 = 'ref', mut2 = NA, mut3 = NA)
    }
    
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
    save_to_logfile(get('log_collection', envir = globalenv()), res, params$logfile)
    
  }
  if (params$clean_after_yourself) {
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
determine_chunklen_compression <- function(start, end) {
  if ((end - start) > 5000 * 1000) {
    chunklen = 100000
  }
  if ((end - start) > 500 * 1000) {
    chunklen = 10000
  } else if ((end-start) < 5000){
    chunklen = 500
  }
    else if (((end - start)) < 50 * 1000) {
    chunklen = 1000
  } else {
    chunklen = 1000
  }
  
  return(chunklen)
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

