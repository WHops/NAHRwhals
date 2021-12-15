#!/usr/local/bin/Rscript

benchmarkfile = "/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/resfile3.txt"

#res = res[1:210,]
library(ggplot2)
library(ggpubr)

res = read.table(benchmarkfile, sep='\t', header=T)
ggplot(res[(res$sim == 0.90),]) + geom_point(aes(x=sdlen, y=id, fill=as.character(chunklen))) + 
  geom_line(aes(x=sdlen, y=id, color=as.character(chunklen))) + theme_bw() +
  scale_x_log10()

ggplot(res[(res$sim == 0.99),]) + geom_point(aes(x=chunklen, y=id, fill=as.character(sim))) + 
  geom_line(aes(x=chunklen, y=id, color=as.character(sdlen))) + theme_bw() +
  scale_x_log10()


plots = list()
count = 1
for (strand in c('+', '-')){
  for (sim in c(0.99, 0.95, 0.90)){
    resp = res[(res$strand==strand) & (res$sim == sim),]
    print(dim(resp))
    plots[[paste0('p', count)]] = 
      ggplot(resp) + geom_tile(aes(x=sdlen, y=chunklen, fill=id)) + 
      scale_fill_viridis_c(limits=c(0.25,1)) + 
      coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
      scale_x_log10(breaks=2**(6:16)) + 
      scale_y_log10(breaks=2**(6:12)) + 
      theme_bw() +
      labs(x='SD length [bp]', y='Alignment chunk lengh [bp]')
    count = count + 1
  }
}



ggarrange(plots$p1, plots$p2, plots$p3,
          plots$p4, plots$p5, plots$p6,
          labels = c('0.99, +', '0.95, +', '0.90, +', '0.99, -', '0.95, -', '0.90, -'),
          ncol = 3, nrow = 2)

