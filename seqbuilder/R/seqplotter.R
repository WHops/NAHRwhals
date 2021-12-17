library(Rcpp)
library(dotplot)
library(Biostrings)
library(ggplot2)

#' A helperfunction for making an exact dotplot, wrapping a functionality of dotplot package. 
#' Should not be used for sequences above 25 kbp. 
#' 
#' @desciption Give me two sequences (character of DNAString I think), 
#' and I will return to you a ggplot object with an exact (bp-precise)
#' dotplot. This relies on the mkDotPlotDataFrame function from 'dotplot' package. 
#'   
#' @param seq1 A DNA sequence, as character or DNAString
#' @param seq2 A DNA sequence, as character or DNAString
#' @return a ggplot2 object: exact dotplot between the two seqs. 
#' 
#' 
#' @author Wolfram Höps
#' @rdname plotting
#' @export
dotPlotr <- function(seq1, seq2, wsize = 5, wstep = 1, nmatch = -1){

  if (length(seq1[1]) > 1)
    stop("seq1 should be provided as a single string")
  if (length(seq2[1]) > 1)
    stop("seq2 should be provided as a single string")
  if (wsize < 1)
    stop("non allowed value for wsize")
  if (wstep < 1)
    stop("non allowed value for wstep")
  if (nmatch < 1)
    nmatch = wsize * 0.5
  if (nmatch > wsize)
    stop("nmatch > wsize is not allowed")
  
  
  wsize = 15
  nmatch = 15
  wstep = 1

  seq2r = as.character(reverseComplement(DNAString(seq2)))
  
  p1 = mkDotPlotDataFrame(seq1, paste0(seq2), wsize, wstep, nmatch)
  p2 = mkDotPlotDataFrame(seq1, seq2r, wsize, wstep, nmatch)

  p2$y = nchar(seq2r) - p2$y
  
  p = ggplot() + geom_point(data=p1, aes(x=x,y=y), shape=15, size=0.5)  + 
    geom_point(data=p2, aes(x=x,y=y), shape=15, size=0.5) + 
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
    theme_bw()

  return(p)
  
}

#' Main wrapper function for making an exact dotplot. Should not be used for sequences
#' above 25 kbp. 
#' 
#' @desciption Give me two sequences (character of DNAString I think), 
#' and I will return to you a ggplot object with an exact (bp-precise)
#' dotplot. This relies on the mkDotPlotDataFrame function from 'dotplot' package. 
#'   
#' @param seq1link [character/link] A link to a target fasta file, as character or DNAString
#' @param seq2link [character/link] A link to a query fasta file, as character or DNAString
#' @param save [bool] if T, save plot to a file. Otherwise, return.
#' @return a ggplot2 object or nothing, depending on the 'save' parameter. 
#' 
#' 
#' @author Wolfram Höps
#' @rdname plotting
#' @export
make_dotplot <- function(targetfa, queryfa, wsize, outfile, save=T, debugmode=F){
  
  
  if (debugmode){
    targetfa = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/res/blub.fa'
    queryfa = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/res/blub.fa'
    outfile = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/res/blub.pdf'
    wsize=10
  }
  
  # Load both
  seq1f = readDNAStringSet(targetfa)
  seq1name = names(seq1f)
  seq1seq = as.character(seq1f)
  
  # Load both
  seq2f = readDNAStringSet(queryfa)
  seq2name = names(seq2f)
  seq2seq = as.character(seq2f)
  
  
  print('sequences loaded')
  p = dotPlotr(seq1seq, seq2seq, wsize)
  p = p + labs(x=seq1name, y=seq2name)
  
    if (save==T){
      ggsave(outfile, plot=p, device='pdf', height=5, width=8, units = 'in')
      print(paste0('Done. Dotplot for simulated sequence written to: ', outfile))
      
  } else if (save==F){
      return(p)
  }
  
}




# runs only when script is run by itself
if (sys.nframe() == 0){
  
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
  
  source_here('seqbuilder_functions.R')
  library(optparse)

  # INPUT
  option_list = list(
    make_option(c("-a", "--seqfile_1"), type="character", default=NULL,
                help="Fastafile 1", metavar="character"),
    make_option(c("-b", "--seqfile_2"), type="character", default=NULL,
                help="Fastafile 2", metavar="character"),
    make_option(c("-o", "--outplot"), type="character", default=NULL,
                help="output plot", metavar="character"),
    make_option(c("-w", "--wsize"), type="numeric", default=10,
                help="Windowsize", metavar="numeric")
  )
  opt <- parse_args(OptionParser(option_list=option_list))
  
  seqfile_1 = opt$seqfile_1
  seqfile_2 = opt$seqfile_2
  outplot = opt$outplot
  wsize = opt$wsize
  
  # Test seqlength
  seqlen1 = as.numeric(nchar(as.character(readDNAStringSet(seqfile_1))))
  seqlen2 = as.numeric(nchar(as.character(readDNAStringSet(seqfile_2))))
  print(seqlen1)
  dotplot_size = seqlen1 * seqlen2
  
  if (dotplot_size < (25000**2)){
    print('making exact plot.')
    make_dotplot(seqfile_1, seqfile_2, 10, outplot)
  } else if (dotplot_size > (25000**2)){
    print('Sequences are too long for exact dotplot (dotplotlimit: 625.000.000 pixels, e.g. 2x25kbp sequences')
  }
  

}
