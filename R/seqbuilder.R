# Whoeps, 24th Nov 2021

#' Introduce repeats to a simulated sequence
#' 
#' This function takes a DNA sequence and an SD instruction (.tsv format), 
#' and introduces this SD. 
#' Also calls add_snps function, to introduce snps in the SD pair, 
#' to create the similarity specified in the input tsv. 
#' 
#' @param seq_base A DNA sequence, as character
#' @param sds A table including columns 'chromStart', 'chromEnd', 'uid', 'strand'
#' @param uid Name of the SD pair to be created. Must appear (exactly) twice in sds.
#' 
#' @examples 
#' 
#' Some example function call or so.
#' 
#' @author Wolfram Höps
#' @rdname seq_modeling
#' @export
model_sds <- function(seq_base, sds, uid){
  
  sds_sub = sds[sds$uid == uid,]
  partner_upstr = substr(seq_base,sds_sub[1,'chromStart'], sds_sub[1,'chromEnd'])
  partner_upstr_snpped = add_snps(partner_upstr, sds_sub[1,'fracMatch'])

  # Copy SD over
  if (sds_sub[2,'strand'] == '+'){
    substr(seq_base,sds_sub[2,'chromStart'], sds_sub[2,'chromEnd']) = partner_upstr_snpped
  } else if (sds_sub[2,'strand'] == '-'){
    substr(seq_base,sds_sub[2,'chromStart'], sds_sub[2,'chromEnd']) = rc(partner_upstr_snpped)
  } else {
    stop("Error in simulating SDs. Check your input SD bedfile.")
  }
  return(seq_base)
}

#' Add snps to a simulated SD pair
#' 
#' Mutate a proportion of bases in a DNA sequence. Can be used to 
#' simulate evolutionary divergence of two SD pairs.
#' 
#' @param seq DNA sequence, as 'character'
#' @param similarity How similar should the result be to the input (fraction of bases, [0-1])
#' 
#' @author Wolfram Höps
#' @rdname seq_modeling
#' @export
add_snps <- function(seq, similarity, seed=1234){
  
  bases= c('A','C','G','T')
  set.seed(seed)
  
  # According to similarity, choose bases to mutate. 
  # Bases can mutate to self, so we add one third to bases to mutate.
  # Going from 'similarity/100' to 'similarity' because now we work with
  # similarity between 0 and 1, not 0 and 100 in the sd file. 
  idx_to_change = sample(nchar(seq), (1-similarity) * (4/3) * nchar(seq))
  for (index in idx_to_change){
    substr(seq, index, index) = sample(bases, 1)
    
  }
  return(seq)
  
}


#' Create a simulated DNA sequence, containing SDs of desired similarity
#' 
#' @param seqlen [numeric] length of simulated sequence in bp
#' @param sdfile [character] A link to an SD file of the form of SD_hg38.tsv tables. 
#' @param outfasta [character] a link to the outputfasta that will be saved
#' 
#' @author Wolfram Höps
#' @rdname seq_modeling
#' @export
simulate_seq <- function(seqlen, sdfile, outfasta, debugmode=F){
  
  # Input: is a tsv file? 
  stopifnot("Error: Input to simulate_seq has to be a .tsv file (29 columns). 
            If you have a .bed file as input, convert it to .tsv with 
            nahrtoolkit::convert_bed_to_tsv inbed outtsv" = endsWith(sdfile, '.tsv'))
  
  
  if (debugmode){
    seqlen = 10000
    sdfile = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/data/sds10y.tsv'
    outfasta = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/res/seq.fa'
  }
  
  # Load and QC SD input file
  sds_colnames = c(
    "chrom", "chromStart", "chromEnd", "name",
    "score","strand", "otherChrom","otherStart", "otherEnd",
    "otherSize", "uid", "posBasesHit", "testResult",
    "verdict", "chits","ccov","alignfile","alignL","indelN",
    "indelS","alignB","matchB","mismatchB", "transitionsB",
    "transversionsB","fracMatch","fracMatchIndel","jcK", "k2K"
  )
  sds = read.table(sdfile); colnames(sds) = sds_colnames

  
  stopifnot("Number of SDs is not a multiple of two" = dim(sds)[1] %% 2 == 0)
  stopifnot("SD coordinates exceed sequence length" = max(sds[,c('chromStart', 'chromEnd', 'otherStart', 'otherEnd')]) < seqlen)
  
  # Get a random sequence of desired length
  seq_base = randDNASeq(seqlen, 0.46)

  # Add the SDs
  seq_modified = seq_base
  for (uid in unique(sds$uid)){
    seq_modified = model_sds(seq_modified, sds, uid)
  }
  
  writeFasta(data.frame(name='sim-sequence', seq=seq_modified),
                        filename=outfasta)
  
  print(paste0('Done. Simulated sequence written to: ', outfasta))
}




