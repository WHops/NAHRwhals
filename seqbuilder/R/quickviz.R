#!/usr/bin/env Rscript

# Quickviz
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

suppressMessages(source_here('seqbuilder.R'))
suppressMessages(source_here('seqbuilder_functions.R'))
suppressMessages(source_here('seqplotter.R'))
suppressMessages(source_here('pafCoordsDotPlotly.R'))
suppressMessages(source_here('seqmodder.R'))
suppressMessages(source_here('compress_paf.R'))
suppressMessages(source_here('paf_to_bed.R'))
suppressMessages(source_here('sd_to_bed.R'))

# INPUT
option_list = list(
  make_option(c("-t", "--targetfasta"), type="character", default=NULL,
              help="Link to a results PAF file", metavar="character"),
  make_option(c("-q", "--queryfasta"), type="character", default=NULL,
              help="Link to an instructions SD file", metavar="character"),
  make_option(c("-c", "--chunklen"), type="numeric", default=1000,
              help="Link to an instructions SD file", metavar="numeric"),
  make_option(c("-o", "--outpaf"), type="character", default=0,
              help="chunklength used", metavar="character"),
  make_option(c("-p", "--outplot"), type="character", default=NULL,
              help="Outfile summarizing resulte", metavar="character")
)
opt <- parse_args(OptionParser(option_list=option_list))

print(opt)
make_chunked_minimap_alnment(opt$targetfasta, opt$queryfasta, 
                             opt$outpaf, opt$outplot)