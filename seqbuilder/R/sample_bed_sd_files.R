#!/usr/local/bin/Rscript


library(optparse)

option_list = list(
  make_option(c("-t", "--seqlen"), type="numeric", default=NULL,
              help="Total sequence length [bp]", metavar="numeric"),
  make_option(c("-s", "--sdlen"), type="numeric", default=NULL,
              help="length of sds", metavar="numeric"),
  make_option(c("-d", "--distance"), type="numeric", default=NULL,
              help="distance between SDs", metavar="numeric"),
  make_option(c("-m", "--fracmatch"), type="numeric", default=NULL, 
              help="SD similarity", metavar="numeric"),
  make_option(c("-o", "--outfile"), type="character", default=NULL, 
              help="Length of chunks to use for minimap2", metavar="character")
  
)


debug=F
if (debug){
  opt = list()
  opt$seqlen = 10000
  opt$distance = 100
  opt$sdlen = 500
  opt$inter_sd_len = 500
}

options(error=traceback)

parser <- OptionParser(usage = "%prog -i alignments.coords -o out [options]",option_list=option_list)
opt = parse_args(parser)
print(opt)
colnames_bed = c('chrom','chromStart',
                      'chromEnd', 'uid', 'otherChrom',
                      'otherStart','otherEnd', 'strand', 'fracMatch')

opt$inter_sd_len = 100
n_sds = floor(opt$seqlen / ((opt$sdlen*3) + opt$distance))
print(n_sds)
contig = "simcontig"
print((0:(n_sds-1)))
print(((opt$sdlen * 2) + opt$inter_sd_len + opt$dist))
chromStart = (0:(n_sds-1)) * ((opt$sdlen * 2) + opt$inter_sd_len + opt$dist)
chromEnd = chromStart + opt$sdlen
otherStart = chromEnd + opt$inter_sd_len
otherEnd = otherStart + opt$sdlen
uid = paste0("SD", 1:n_sds)

print(chromStart)
sds = data.frame('chrom' = contig, 
                     'chromStart' = chromStart,
                     'chromEnd' = chromEnd, 
                     'uid' = uid, 
                     'otherChrom' = contig,
                     'otherStart' = otherStart, 
                     'otherEnd' = otherEnd, 
                     'strand' = '+',
                     'fracmatch' = opt$fracmatch)

write.table(sds, file=opt$outfile, sep='\t', row.names=F, col.names=F, quote=F)
