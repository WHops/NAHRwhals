# 

#from https://bootstrappers.umassmed.edu/guides/main/r_writeFasta.html
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


randDNASeq <- function(n, gcfreq, seed=1234){
  bases= c('A','C','G','T')
  
  set.seed(seed)
  seq = sample(bases, n, replace=T, 
               prob = c((1-gcfreq)/2, gcfreq/2, gcfreq/2, (1-gcfreq)/2) )
  return(paste(seq, collapse=''))
}