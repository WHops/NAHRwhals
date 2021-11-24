#!/bin/Rscript

source_here <- function(x, ...) {
  dir <- "."
  if(sys.nframe()>0) {
    frame <- sys.frame(1)
    if (!is.null(frame$ofile)) {
      dir <- dirname(frame$ofile)
    }
  }
  source(file.path(dir, x), ...)
}

library(optparse)
#sourcedir = ("/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/R/")
#source(here::here('seqbuilder.R'))

source_here('seqbuilder.R')
source_here('seqbuilder_functions.R')
source_here('seqplotter.R')



# sourcefiles = list.files(sourcedir)
# sourcefiles = sourcefiles[sourcefiles != 'seqbuilder_wrapper.R']
# sapply(paste0(sourcedir,sourcefiles), source)

# INPUT
option_list = list(
  make_option(c("-l", "--seqlen"), type="numeric", default=NULL,
              help="Total sequence length [bp]", metavar="numeric"),
  make_option(c("-s", "--sdfile"), type="character", default=NULL,
              help="Bedfile with desired SDs", metavar="character"),
  make_option(c("-o", "--outfasta"), type="character", default="./outputcorr/",
              help="Output fasta", metavar="character"),
  make_option(c("-p", "--outplot"), type="character", default="./outputcorr/",
              help="Output simple dotplot", metavar="character")
  )

opt <- parse_args(OptionParser(option_list=option_list))

seqlen = opt$seqlen
sdfile = opt$sdfile
outfasta = opt$outfasta
outplot = opt$outplot



simulate_seq(seqlen, sdfile, outfasta)

# Make an exact dotplot, in case this is feasible. 
# Otherwise (or additionally?) we are going to make a minimap2 dotplot. 
if (seqlen <= 25000){
  print('making exact plot.')
  make_dotplot(outfasta, outfasta, 15, outplot)
} else if (seqlen > 25000){
  print('Sequence1 is too long for exact dotplot (>25 kbp)')
}

