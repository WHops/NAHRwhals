# Whoeps, 24th Nov 2021


model_sds <- function(seq_base, sds, sds_name){
  
  sds_sub = sds[sds$sdname == sds_name,]
  partner_upstr = substr(seq_base,sds_sub[1,'start'], sds_sub[1,'end'])
  partner_upstr_snpped = add_snps(partner_upstr, sds[1,'similarity'])

  # Copy SD over
  if (sds_sub[2,'strand'] == '+'){
    substr(seq_base,sds_sub[2,'start'], sds_sub[2,'end']) = partner_upstr
  } else if (sds_sub[2,'strand'] == '-'){
    substr(seq_base,sds_sub[2,'start'], sds_sub[2,'end']) = rc(partner_upstr)
  } else {
    stop("Error in simulating SDs. Check your input SD bedfile.")
  }
  return(seq_base)
}

add_snps <- function(seq, similarity){
  
  bases= c('A','C','G','T')
  
  # According to similarity, choose bases to mutate. 
  # Bases can mutate to self, so we add one third to bases to mutate.
  idx_to_change = sample(nchar(seq), (1-(similarity/100)) * (4/3) * nchar(seq))
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
    sdfile = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/data/sds.bed'
    outfasta = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/res/seq.fa'
  }
  
  bed_colnames = c('seqname','start','end','sdname','similarity','strand')
  sds = read.table(sdfile); colnames(sds) = bed_colnames

  stopifnot("Number of SDs is not a multiple of two" = dim(sds)[1] %% 2 == 0)
  stopifnot("SD coordinates exceed sequence length" = max(sds[,c('start', 'end')]) < seqlen)
  
  # Get a random sequence of desired length
  seq_base = randDNASeq(seqlen, 0.46)

  # Add the SDs
  seq_modified = seq_base
  for (sdname in unique(sds$sdname)){
    seq_modified = model_sds(seq_modified, sds, sdname)
  }
  
  writeFasta(data.frame(name='sim-sequence', seq=seq_modified),
                        filename=outfasta)
  
  print(paste0('Done. Simulated sequence written to: ', outfasta))
}


#sdfile_large = "/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/data/sds_large.bed"
#outfasta = "/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/res/seq.fa"
#simulate_seq(1000000, sdfile_large, outfasta)




