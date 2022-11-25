

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
  } else if (((end - start)) < 50 * 1000) {
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
  save_to_logfile(get('log_collection', envir = globalenv()), res, params$logfile)
  
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
  
  if (length(xstart) > 433){
    print('Too many segments to make colored plot')
    return()
  }
  
  colors = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  
  n <- length(xstart)
  qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  
  plot_x_y_segmented = plot_x_y + 
    ggplot2::geom_rect(data=datx,
                       ggplot2::aes(xmin=xstart, xmax=xend, ymin=0, ymax=ymax, fill=col_vector[1:length(xstart)]),
                       alpha=0.5,) + 
    ggplot2::guides(fill = FALSE)#  +
  # geom_segment(data=daty,
  #              aes(x=0, xend=xmax, y=ystart, yend=ystart), color='grey')
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