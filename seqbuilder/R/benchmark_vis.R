#!/usr/local/bin/Rscript

#' Run a benchmark visualization. Thrown into a function but is really a script. 
#' 
#' @description TO BE REVIEWED / CORRECTED!
#' 
#' 
#' @author Wolfram Höps
#' @rdname evaluation
#' @export
benchmark_vis <- function(){
  benchmarkfile = "/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/resfile3.txt"
  
  #res = res[1:210,]
  
  
  res = read.table(benchmarkfile, sep='\t', header=T)
  ggplot2::ggplot(res[(res$sim == 0.90),]) + 
    ggplot2::geom_point(ggplot2::aes(x=sdlen, y=id, fill=as.character(chunklen))) + 
    ggplot2::geom_line(ggplot2::aes(x=sdlen, y=id, color=as.character(chunklen))) + 
    ggplot2::theme_bw() +
    ggplot2::scale_x_log10()
  
  ggplot2::ggplot(res[(res$sim == 0.99),]) + 
    ggplot2::geom_point(ggplot2::aes(x=chunklen, y=id, fill=as.character(sim))) + 
    ggplot2::geom_line(ggplot2::aes(x=chunklen, y=id, color=as.character(sdlen))) + 
    ggplot2::theme_bw() +
    ggplot2::scale_x_log10()
  
  
  plots = list()
  count = 1
  for (strand in c('+', '-')){
    for (sim in c(0.99, 0.95, 0.90)){
      resp = res[(res$strand==strand) & (res$sim == sim),]
      print(dim(resp))
      plots[[paste0('p', count)]] = 
        ggplot2::ggplot(resp) + 
        ggplot2::geom_tile(ggplot2::aes(x=sdlen, y=chunklen, fill=id)) + 
        ggplot2::scale_fill_viridis_c(limits=c(0.25,1)) + 
        ggplot2::coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
        ggplot2::scale_x_log10(breaks=2**(6:16)) + 
        ggplot2::scale_y_log10(breaks=2**(6:12)) + 
        ggplot2::theme_bw() +
        ggplot2::labs(x='SD length [bp]', y='Alignment chunk lengh [bp]')
      count = count + 1
    }
  }
  
  ggpubr::ggarrange(plots$p1, plots$p2, plots$p3,
                    plots$p4, plots$p5, plots$p6,
                    labels = c('0.99, +', '0.95, +', '0.90, +', '0.99, -', '0.95, -', '0.90, -'),
                    ncol = 3, nrow = 2)
  
}



