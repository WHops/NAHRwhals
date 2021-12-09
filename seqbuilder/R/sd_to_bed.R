#!/usr/local/bin/Rscript

# Functions

sd_to_bed <- function(sdlink, outbedfile=NULL){
  
  library(dplyr)
  sd = read.table(sdlink, sep='\t')
  
  colnames(sd) =  c(
    "chrom", "chromStart", "chromEnd", "name",
    "score","strand", "otherChrom","otherStart", "otherEnd",
    "otherSize", "uid", "posBasesHit", "testResult",
    "verdict", "chits","ccov","alignfile","alignL","indelN",
    "indelS","alignB","matchB","mismatchB", "transitionsB",
    "transversionsB","fracMatch","fracMatchIndel","jcK", "k2K"
  )
  
  sd_keep = sd[,c('chrom','chromStart','chromEnd', 'uid',
                 'otherChrom', 'otherStart', 'otherEnd', 
                 'strand', 'fracMatch')]
  
  print(sd_keep)
  sd_keep_uniq = sd_keep %>% group_by(uid) %>% slice(1)
  sd_keep_uniq_sort = sd_keep_uniq[order(sd_keep_uniq$chromStart),]
  
  if (is.null(outbedfile)){
    return(sd_keep_uniq_sort)
  } else if ( !is.null(outbedfile) ){
    write.table(sd_keep_uniq_sort, file=outbedfile, sep='\t', row.names=F, col.names=F, quote=F)
  }
}


