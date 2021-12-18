#!/usr/local/bin/Rscript

# Functions

#' Convert a (normal, double-entry-) tsv file to bed
#' 
#' @param sdlink [character/link] a link to the input sd file to be converted to bed.
#' @param outbedfile [character/link/NULL] a link to the output bedfile. If null, a 
#' bed-like dataframe is returned
#' @return a bed-file style dataframe if outbedfile == NULL
#' 
#' @author Wolfram Höps
#' @rdname format_sd_to_bed
#' @export
sd_to_bed <- function(sdlink, outbedfile=NULL){
  
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
  
  #sd_keep_uniq = sd_keep %>% group_by(uid) %>% slice(1)
  sd_keep_uniq = dplyr::slice(dplyr::group_by(sd_keep, uid), 1)
  
  sd_keep_uniq_sort = sd_keep_uniq[order(sd_keep_uniq$chromStart),]
  
  if (is.null(outbedfile)){
    return(sd_keep_uniq_sort)
  } else if ( !is.null(outbedfile) ){
    write.table(sd_keep_uniq_sort, file=outbedfile, sep='\t', row.names=F, col.names=F, quote=F)
  }
}


