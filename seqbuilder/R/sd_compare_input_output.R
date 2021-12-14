#!/usr/local/bin/Rscript

# Compare output paf with input paf


turn_intervals_to_vector_0_1 <- function(intervals, maxlen){
  # Inclusive. tstart and tend are both included. 
  
  sd_coverage_bool = 0*c(1:maxlen)
  for (row in 2:dim(intervals)[1]){
    sd_coverage_bool[intervals[row,1]:intervals[row,2]] = 1
  }
  return(sd_coverage_bool)
}

compare_sd_paf_input_output <- function(preslink, ptruelink, chunklen, runname, outfile){
  
  
  debug = F
  if (debug){
    preslink = "/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/res/paf/test_chunked.paf"
    ptruelink = "/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/res/test.tsv"
  }
  
  pres = read.table(preslink, sep='\t')
  ptrue = read.table(ptruelink, sep='\t')
  
  colnames_paf = c('qname','qlen','qstart','qend',
                   'strand','tname','tlen','tstart',
                   'tend','nmatch','alen','mapq')
  colnames_sd =  c(
    "chrom", "chromStart", "chromEnd", "name",
    "score","strand", "otherChrom","otherStart", "otherEnd",
    "otherSize", "uid", "posBasesHit", "testResult",
    "verdict", "chits","ccov","alignfile","alignL","indelN",
    "indelS","alignB","matchB","mismatchB", "transitionsB",
    "transversionsB","fracMatch","fracMatchIndel","jcK", "k2K"
  )
  
  colnames(pres) = colnames_paf
  colnames(ptrue) = colnames_sd
  
  maxlen = max(c(pres$tend, pres$tstart, ptrue$chromStart, ptrue$chromEnd))
  pres_vector = turn_intervals_to_vector_0_1(pres[,c('tstart', 'tend')], maxlen)
  ptrue_vector = turn_intervals_to_vector_0_1(ptrue[,c('chromStart', 'chromEnd')], maxlen)
  
  both = data.frame(res = pres_vector, true = ptrue_vector)
  both$identical = both$res == both$true
  
  out_line = data.frame(
    runname = runname,
    sdlen = ptrue[1,'chromEnd'] - ptrue[1,'chromStart'],
    sim = ptrue[1,'fracMatch'],
    chunklen = chunklen,
    id = round(sum(both$identical) / length(both$identical),3),
    tp = round((sum(both[both$true==T,]$identical) / 
                length(both[both$true==T,]$identical)), 3),
    tn = round((sum(both[both$true==F,]$identical) / 
            length(both[both$true==F,]$identical)), 3)
  )
  
  # Save
  write.table(x=out_line, file=outfile, sep = "\t", col.names = !file.exists(outfile), 
              append = T, quote = F, row.names=F)
}


# runs only when script is run by itself
if (sys.nframe() == 0){
  
  library(optparse)

    # INPUT
  option_list = list(
    make_option(c("-p", "--preslink"), type="character", default=NULL,
                help="Link to a results PAF file", metavar="character"),
    make_option(c("-q", "--ptruelink"), type="character", default=NULL,
                help="Link to an instructions SD file", metavar="character"),
    make_option(c("-c", "--chunklen"), type="character", default=0,
                help="chunklength used", metavar="character"),
    make_option(c("-n", "--runname"), type="character", default=NULL,
                help="Name that the result will apear in outfile", metavar="character"),
    make_option(c("-o", "--outfile"), type="character", default=NULL,
                help="Outfile summarizing resulte", metavar="character")
  )
  opt <- parse_args(OptionParser(option_list=option_list))
  
  preslink = opt$preslink
  ptruelink = opt$ptruelink
  chunklen = opt$chunklen
  runname = opt$runname
  outfile = opt$outfile
  
  compare_sd_paf_input_output(preslink, ptruelink, chunklen, runname, outfile)
  
  
}


