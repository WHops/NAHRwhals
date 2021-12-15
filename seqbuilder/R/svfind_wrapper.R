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

# runs only when script is run by itself
if (sys.nframe() == 0){
  
  library(optparse)
  
  suppressMessages(source_here('seqbuilder.R'))
  suppressMessages(source_here('seqbuilder_functions.R'))
  suppressMessages(source_here('seqplotter.R'))
  suppressMessages(source_here('pafCoordsDotPlotly.R'))
  suppressMessages(source_here('seqmodder.R'))
  suppressMessages(source_here('compress_paf.R'))
  suppressMessages(source_here('paf_to_bed.R'))
  suppressMessages(source_here('sd_to_bed.R'))
  
  
  # INPUT

  debug=T
  if (debug){
    opt = list()
    opt$input_fasta = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/res/fa/simple_inv_mut.fa'
    opt$reference = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/res/fa/simple_inv.fa'
    
  } else {
    opt <- parse_args(OptionParser(option_list=option_list))
    option_list = list(
      make_option(c("-i", "--input_fasta"), type="character", default=NULL,
                  help="input fasta", metavar="character"),
      make_option(c("-r", "--reference"), type="character", default='hg38',
                  help="Reference sequence. Can be 'hg38', 'chm13' or a custom input fasta.", metavar="character"),
      make_option(c("-o", "--outfile"), type="character", default=NULL,
                  help="Outputfile with SVs", metavar="character"),
      make_option(c("-c", "--chunklen"), type="numeric", default=512,
                  help="Chunklen", metavar="numeric"),
    )
    
  }
  
  # Find coordinates matching best the input. This has to be coded still. For now 
  # We give these manually. 
  #coords_start_end = find_matching_coordinates_dummy(opt$input_fasta, opt$reference)
  
  coords_start_end = c(0, 10000)
  
  outpaf = 'test/blub.paf'
  outplot = 'test/plot.pdf'
  outbed = 'test/blub.bed'
  # Find SD possibilities
  #make_chunked_minimap_alnment(opt$input_fasta, opt$input_fasta, outpaf, outplot, 
  #                             chunklen = 100)
  
  # Write bedfile
  paf_write_bed(outpaf, outbed)
  
  
}
