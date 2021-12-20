
#' A helperfunction for making an exact dotplot, wrapping a functionality of dotplot package. 
#' Should not be used for sequences above 25 kbp. 
#' 
#' @description Give me two sequences (character of DNAString I think), 
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
  
  # Yes these shouldn't be loaded directly. But anything else crashes the code, 
  # So we leave it like that FOR NOW. 
  library(Rcpp)
  library(dotplot)
  
  wsize = 15
  nmatch = 15
  wstep = 1

  seq2r = as.character(Biostrings::reverseComplement(Biostrings::DNAString(seq2)))

  
  p1 = mkDotPlotDataFrame(seq1, paste0(seq2), wsize, wstep, nmatch)
  p2 = dotplot::mkDotPlotDataFrame(seq1, seq2r, wsize, wstep, nmatch)

  p2$y = nchar(seq2r) - p2$y
  
  p = ggplot2::ggplot() + 
        ggplot2::geom_point(data=p1, ggplot2::aes(x=x,y=y), shape=15, size=0.5)  + 
        ggplot2::geom_point(data=p2, ggplot2::aes(x=x,y=y), shape=15, size=0.5) + 
        ggplot2::coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
        ggplot2::theme_bw()

  return(p)
  
}

#' Main wrapper function for making an exact dotplot. Should not be used for sequences
#' above 25 kbp. 
#' 
#' @description Give me two sequences (character of DNAString I think), 
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
  seq1f = Biostrings::readDNAStringSet(targetfa)
  seq1name = names(seq1f)
  seq1seq = as.character(seq1f)
  
  # Load both
  seq2f = Biostrings::readDNAStringSet(queryfa)
  seq2name = names(seq2f)
  seq2seq = as.character(seq2f)
  
  
  print('sequences loaded')
  p = dotPlotr(seq1seq, seq2seq, wsize)
  p = p + ggplot2::labs(x=seq1name, y=seq2name)
  
    if (save==T){
      ggplot2::ggsave(outfile, plot=p, device='pdf', height=5, width=8, units = 'in')
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
  
  #source_here('seqbuilder_functions.R')

  # INPUT
  option_list = list(
    optparse::make_option(c("-a", "--seqfile_1"), type="character", default=NULL,
                help="Fastafile 1", metavar="character"),
    optparse::make_option(c("-b", "--seqfile_2"), type="character", default=NULL,
                help="Fastafile 2", metavar="character"),
    optparse::make_option(c("-o", "--outplot"), type="character", default=NULL,
                help="output plot", metavar="character"),
    optparse::make_option(c("-w", "--wsize"), type="numeric", default=10,
                help="Windowsize", metavar="numeric")
  )
  opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
  
  seqfile_1 = opt$seqfile_1
  seqfile_2 = opt$seqfile_2
  outplot = opt$outplot
  wsize = opt$wsize
  
  # Test seqlength
  seqlen1 = as.numeric(nchar(as.character(Biostrings::readDNAStringSet(seqfile_1))))
  seqlen2 = as.numeric(nchar(as.character(Biostrings::readDNAStringSet(seqfile_2))))
  print(seqlen1)
  dotplot_size = seqlen1 * seqlen2
  
  if (dotplot_size < (25000**2)){
    print('making exact plot.')
    make_dotplot(seqfile_1, seqfile_2, 10, outplot)
  } else if (dotplot_size > (25000**2)){
    print('Sequences are too long for exact dotplot (dotplotlimit: 625.000.000 pixels, e.g. 2x25kbp sequences')
  }
  

}
