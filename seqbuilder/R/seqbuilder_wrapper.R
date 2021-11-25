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

source_here('seqbuilder.R')
source_here('seqbuilder_functions.R')
source_here('seqplotter.R')
source_here('pafCoordsDotPlotly.R')
source_here('seqmodder.R')


# sourcefiles = list.files(sourcedir)
# sourcefiles = sourcefiles[sourcefiles != 'seqbuilder_wrapper.R']
# sapply(paste0(sourcedir,sourcefiles), source)

# INPUT
option_list = list(
  make_option(c("-l", "--seqlen"), type="numeric", default=NULL,
              help="Total sequence length [bp]", metavar="numeric"),
  make_option(c("-s", "--sdfile"), type="character", default=NULL,
              help="Bedfile with desired SDs", metavar="character"),
  make_option(c("-v", "--SVfile"), type="character", default=NULL,
              help="Textfile with SVs to model", metavar="character"),
  make_option(c("-o", "--outprefix"), type="character", default="./outputcorr/", 
              help="Output fasta", metavar="character")
  
  )

opt <- parse_args(OptionParser(option_list=option_list))

seqlen = opt$seqlen
sdfile = opt$sdfile
svfile = opt$SVfile
outprefix = opt$outprefix


# Prepare the output directories
fadir = "../res/fa/"
pafdir = "../res/paf/"
plotdir = "../res/plot/"

dir.create(fadir)
dir.create(pafdir)
dir.create(plotdir)

outfasta = paste0(fadir, outprefix, '.fa')
outpaf = paste0(pafdir, outprefix, '.paf')
outplot1 = paste0(plotdir, outprefix, '.pdf')
outplot2 = paste0(plotdir, outprefix, 'minimap')
outmutfasta = paste0(fadir, outprefix, '_mut.fa')
outmutsd = paste0(fadir, outprefix, '_mut.bed')
outmutpaf = paste0(pafdir, outprefix, '_mut.paf')
outplot3 = paste0(plotdir, outprefix, '_mut.pdf')
outplot4 = paste0(plotdir, outprefix, 'minimap_mut')


## MAKE THE SEQUENCE
simulate_seq(seqlen, sdfile, outfasta)

# Make an exact dotplot, in case this is feasible.
# Otherwise (or additionally?) we are going to make a minimap2 dotplot.
if (seqlen <= 25000){
  print('making exact plot.')
  make_dotplot(outfasta, outfasta, 15, outplot1)
} else if (seqlen > 25000){
  print('Sequence1 is too long for exact dotplot (>25 kbp)')
}

## MAKE A MINIMAP2+DOTPLOTLY dotplot
system(paste0("minimap2 -x asm20 -c -z400,50 -s 0 -M 0.2 -N 100 -P --hard-mask-level ", outfasta, " ", outfasta, " > ", outpaf))
system(paste0("./pafCoordsDotPlotly.R -i ", outpaf,
                                     " -o ", outplot2,
                                     " -m 50 -p 10 -q 10"))


## MUTATE THE SEQUENCE
mutate_seq(outfasta, sdfile, svfile, outmutfasta, outmutsd)

# Make an exact dotplot, in case this is feasible.
# Otherwise (or additionally?) we are going to make a minimap2 dotplot.
if (seqlen <= 25000){
  print('making exact plot.')
  make_dotplot(outfasta, outmutfasta, 15, outplot3)
} else if (seqlen > 25000){
  print('Sequence1 is too long for exact dotplot (>25 kbp)')
}

## MAKE A MINIMAP2+DOTPLOTLY dotplot
system(paste0("minimap2 -x asm20 -c -z400,50 -s 0 -M 0.2 -N 100 -P --hard-mask-level ", outfasta, " ", outmutfasta, " > ", outmutpaf))
system(paste0("./pafCoordsDotPlotly.R -i ", outmutpaf,
              " -o ", outplot4,
              " -m 50 -p 10 -q 10"))


