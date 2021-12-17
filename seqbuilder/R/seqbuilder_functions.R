#' A core wrapper function. Give me two fasta files and I'll give you
#' an output paf and a dotplot pdf. 
#' 
#' @description This is sitting at the core of the nahrchainer module. It ties
#' together core modules, such as query-sequence chunking, minimap2 alignment, 
#' paf recombination and plotting. s
#' 
#' @param targetfasta [character/link] link to the 'target' single-sequence fasta (sometimes reference, e.g. chm13.)
#' @param queryfasta [character/link] link to the 'query' single-sequence fasta.
#' @param chunklen [numeric] length of sequence chunks to split the query to (default = 1000 bp.)
#' @param hllink [character/link] link to an SD annotation file (bed, tsv or paf) to include as 
#' highlights in the plot.
#' @param hltype [character] filetype of hllink. Can be 'NULL', 'bed', 'tsv', 'paf'.
#' @param outpaf [character/link] Path to the output paffile to be written.
#' @param outplot [character/link] Path to the output plot to be written.
#' 
#' @param keep_ref [numeric] Number of alignments to keep or something. Plotting parameter. 
#' @param plot_size [numeric] Plot size in inch. 
#' @param outsd [character/link] Outputfile for the modified SD tsv file. 
#' @return nothing. But output files written. 
#' 
#' @author Wolfram Höps
#' @rdname alignment
#' @export
make_chunked_minimap_alnment <- function(targetfasta, queryfasta, outpaf, outplot, 
                                         chunklen = 1000, keep_ref = 10000, 
                                         plot_size = 10, 
                                         hllink = hllink,
                                         hltype = NULL){
  
  # Define intermediate files
  queryfasta_chunk = paste0(queryfasta, ".chunk.fa")
  outpaf_chunk = paste0(outpaf, '.chunk')
  outpaf_awk = paste0(outpaf, '.awked')
  
  # Run a series of chunking, aligning and merging functions/scripts
  
  # Single-sequence query fasta gets chopped into pieces.
  shred_seq(queryfasta, queryfasta_chunk, chunklen)
  
  # Self explanatory
  run_minimap2(targetfasta, queryfasta_chunk, outpaf_chunk)
  
  # Awk is used to correct the seuqence names. This is because I know only there
  # how to use regex...
  awk_edit_paf(outpaf_chunk, outpaf_awk)
  
  # paf of fragmented paf gets put back together. 
  compress_paf_fnct(outpaf_awk, outpaf)
  
  # Make a dotplot of that final paf (and with sd highlighting). 
  pafdotplot_make(outpaf, outplot, keep_ref=keep_ref, plot_size=plot_size,
                  hllink = hllink, hltype = hltype)
  
}

#' Chunkify query fasta
#' 
#' @description This is a helperfunction calling an external script to 
#' chop a query sequence into chunks. 
#' 
#' @param infasta [character/link] single-seq fasta to be chopped
#' @param outfasta_chunk [character/link] output chopped multi-seq fasta.
#' @param chunklen [numeric] length of sequence chunks in bp
#' @param scriptloc [character/link] link to shred.ss from bbmap.
#' @return nothing. Only output files written. 
#' 
#' @author Wolfram Höps
#' @rdname alignment
#' @export
shred_seq <- function(infasta, outfasta_chunk, chunklen, scriptloc='../../../bbmap/shred.sh'){
  print(paste0(scriptloc," in=", infasta, " out=", outfasta_chunk, " length=", chunklen))
  system(paste0(scriptloc," in=", infasta, " out=", outfasta_chunk, " length=", chunklen))
}

#' Submit a system command to run minimap2
#' 
#' @description This is a helperfunction to run minimap2
#' 
#' @param targetfasta [character/link] link to the 'target' single-sequence fasta (sometimes reference, e.g. chm13.)
#' @param queryfasta [character/link] link to the 'query' fasta. Can be single or multi-fasta
#' @param outpaf [character/link] Path to the output paffile to be written.
#' @param minimap2loc [character/link] link to minimap2 binary.

#' @return nothing. Only output files written. 
#' 
#' @author Wolfram Höps
#' @rdname alignment
#' @export
run_minimap2 <- function(targetfasta, queryfasta, outpaf, minimap2loc = "/Users/hoeps/opt/anaconda3/bin/minimap2"){
  #system(paste0(minimap2loc," -x asm20 -c -z400,50 -s 0 -M 0.2 -N 100 -P --hard-mask-level ", fastatarget, " ", fastaquery, " > ", outpaf))
  system(paste0(minimap2loc," -x asm20 -P -c -s 0 -M 0.2 ", targetfasta, " ", queryfasta, " > ", outpaf))
  
}


#' Submit a system command to run awk to change sequence chunk names
#' 
#' @description This is a helperfunction to run awk
#' 
#' @param inpaf [character/link] input paf 
#' @param outpaf [character/link] output paf
#' 
#' @return nothing. Only output files written. 
#' 
#' @author Wolfram Höps
#' @rdname alignment
#' @export 
awk_edit_paf <- function(inpaf, outpaf){
  system(paste0("../scripts/awk_on_paf.sh ", inpaf, " ", outpaf))
}

#' Helperfunction to save a fasta file. 
#' 
#' @description a simple function that takes a data.frame that has a column name
#' and seq and writes a fasta file from it. Taken from
#' https://bootstrappers.umassmed.edu/guides/main/r_writeFasta.html.
#' 
#' @param data [character/link] data frame with column 'name' and 'seq'
#' @param filename [character/link] output filename. 
#' 
#' @return nothing. Only output files written. 
#' 
#' @author Nicholas Hathaway
#' @rdname alignment
#' @export 
writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

#' https://rdrr.io/cran/insect/src/R/complement.R
#' Reverse complement DNA in character string format.
#'
#' This function reverse complements a DNA sequence or vector of DNA
#'   sequences that are stored as character strings.
#'
#' @param z a vector of DNA sequences in upper case character string format.
#' @return a vector of DNA sequences as upper case character strings.
#' @details This function accepts only DNA sequences in concatenated character
#'   string format, see \code{\link[ape]{complement}} in the \code{\link[ape]{ape}}
#'   package for "DNAbin" input objects, and \code{\link[seqinr]{comp}} in the
#'   \code{\link[seqinr]{seqinr}} package for when the input object is a character
#'   vector.
#' @author Shaun Wilkinson
#' @examples rc("TATTG")
################################################################################
rc <- function(z){
  rc1 <- function(zz){
    s <- strsplit(zz, split = "")[[1]]
    s <- rev(s)
    dchars <- strsplit("ACGTMRWSYKVHDBNI", split = "")[[1]]
    comps <- strsplit("TGCAKYWSRMBDHVNI", split = "")[[1]]
    s <- s[s %in% dchars] # remove spaces etc
    s <- dchars[match(s, comps)]
    s <- paste0(s, collapse = "")
    return(s)
  }
  z <- toupper(z)
  tmpnames <- names(z)
  res <- unname(sapply(z, rc1))
  if(!is.null(attr(z, "quality"))){
    strev <- function(x) sapply(lapply(lapply(unname(x), charToRaw), rev), rawToChar)
    attr(res, "quality") <- unname(sapply(attr(z, "quality"), strev))
  }
  names(res) <- tmpnames
  return(res)
}
################################################################################

#' Create a simple random sequence. 
#' 
#' @description Simple function to create a random string of ACGT. 
#' 
#' @param n [numeric] length of desired sequence (bp)
#' @param gcfreq [character/link] desired GC frequency. 
#' @return A character vector of a random DNA sequence. 
#' 
#' @author Wolfram Höps
#' @rdname seq_modeling
#' @export 
randDNASeq <- function(n, gcfreq, seed=1234){
  bases= c('A','C','G','T')
  
  set.seed(seed)
  seq = sample(bases, n, replace=T, 
               prob = c((1-gcfreq)/2, gcfreq/2, gcfreq/2, (1-gcfreq)/2) )
  return(paste(seq, collapse=''))
}