# Whoeps, 24th Nov 2021


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


#' Create a simulated sequence
#' @export
simulate_seq <- function(seqlen, sdfile, outfasta, debugmode=F){
  
  if (debugmode){
    seqlen = 10000
    sdfile = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/data/sds10y.tsv'
    outfasta = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/res/seq.fa'
  }
  
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


#sdfile_large = "/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/data/sds_large.bed"
#outfasta = "/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/res/seq.fa"
#simulate_seq(1000000, sdfile_large, outfasta)




