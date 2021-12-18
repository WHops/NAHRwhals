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
    optparse::make_option(c("-l", "--seqlen"), type="numeric", default=NULL,
                help="Total sequence length [bp]", metavar="numeric"),
    optparse::make_option(c("-s", "--sdfile"), type="character", default=NULL,
                help="TSV file with desired SDs", metavar="character"),
    optparse::make_option(c("-v", "--SVfile"), type="character", default=NULL,
                help="Textfile with SVs to model", metavar="character"),
    optparse::make_option(c("-o", "--outprefix"), type="character", default="./outputcorr/", 
                help="Output fasta", metavar="character"),
    optparse::make_option(c("-c", "--chunklen"), type="numeric", default=1000, 
                help="Length of chunks to use for minimap2", metavar="character")
    
    )
  
  debug=F
  if (debug){
    opt = list()
    seqlen = 10000
    sdfile = "../data/sds10y.tsv"
    svfile = "../data/svs10.txt"
    outprefix = "debug"
    chunklen = 1000
  } else {
    opt <- parse_args(OptionParser(option_list=option_list))
    
    seqlen = opt$seqlen
    sdfile = opt$sdfile
    svfile = opt$SVfile
    chunklen = opt$chunklen
    outprefix = opt$outprefix
  }
  
  # Prepare the output directories

  fadir = "../res/fa/"
  pafdir = "../res/paf/"
  plotdir = "../res/plot/"
  beddir = "../res/bed/"
  
  for (dir in c('fa', 'paf', 'plot', 'bed')){
    dir.create(paste0("../res/",dir,"/"), showWarnings=F)
  }


  # Files containing simulated sequence before mut
  outfasta = paste0(fadir, outprefix, '.fa')
  outpaf = paste0(pafdir, outprefix, '_chunked.paf')
  outbed = paste0(beddir, outprefix, '.bed')
  
  outplot1 = paste0(plotdir, outprefix, '.pdf')
  outplot2 = paste0(plotdir, outprefix, 'minimap')
  
  # Files containing simulated sequence after mut
  outmutfasta = paste0(fadir, outprefix, '_mut.fa')
  outmutsd = paste0(pafdir, outprefix, '_mut_self_gt.paf')
  outmutpaf = paste0(pafdir, outprefix, '_mut.paf')
  outmutpaf_self = paste0(pafdir, outprefix, '_mut_self.paf')
  outmutbed = paste0(beddir, outprefix, '_mut.bed')
  outmutbedgt = paste0(beddir, outprefix, '_mut_gt.bed')
  
  outplot3 = paste0(plotdir, outprefix, '_mut.pdf')
  outplot4 = paste0(plotdir, outprefix, 'minimap_mut')
  outplot5 = paste0(plotdir, outprefix, '_mut.self.pdf')
  outplot6 = paste0(plotdir, outprefix, '_mut.self.gt.pdf')
  
  
  ## MAKE THE SEQUENCE
  simulate_seq(seqlen, sdfile, outfasta)
  
  # Make an exact dotplot, in case this is feasible.
  # Otherwise (..additionally) we are going to make a minimap2 dotplot.
  if (seqlen <= 25000){
    make_dotplot(outfasta, outfasta, 15, outplot1)
  }
  
  # Make minimap2 alignment
  make_chunked_minimap_alnment(outfasta, outfasta, outpaf, outplot2, 
                               chunklen = chunklen,
                               hllink = opt$sdfile, hltype = 'sd')
  
  # Write bedfile
  paf_write_bed(outpaf, outbed)
  
  
  if (is.null(svfile)){
    print("No Mutation input. Only returning primary sequences")
  } else {
      
    ## MUTATE THE SEQUENCE
    mutate_seq(outfasta, sdfile, svfile, outmutfasta, outmutsd)
    sd_to_bed(outmutsd, outmutbedgt)
    
    
    # Make an exact dotplot, in case this is feasible.
    # Otherwise (or additionally?) we are going to make a minimap2 dotplot.
    if (seqlen <= 25000){
      make_dotplot(outfasta, outmutfasta, 15, outplot3)
      make_dotplot(outmutfasta, outmutfasta, 15, outplot6)
      
    }
    
    ## MAKE A MINIMAP2+DOTPLOTLY dotplot
    make_chunked_minimap_alnment(outfasta, outmutfasta, outmutpaf, outplot4, chunklen = chunklen)
    
    ## MAKE A MINIMAP2+DOTPLOTLY dotplot with self
    print('here')
    make_chunked_minimap_alnment(outmutfasta, outmutfasta, outmutpaf_self, outplot5, 
                                 chunklen = chunklen, hllink = outmutsd, hltype = 'paf')
    # Write bedfile
    paf_write_bed(outmutpaf_self, outmutbed)
  }
}

