# Whoeps, 25th Nov 2021
library(Biostrings)

debug = F

mutate_seq <- function(seq_fasta, sds_bed, sv_instr_txt, outfasta, outsd){
  
  if (debug){
    seq_fasta = "../res/fa/s10.fa"
    sds_bed = "../data/sds10x.bed"
    sv_instr_txt = "../data/svs10.txt"
  }
  
  # Load both sequence
  seqf = readDNAStringSet(seq_fasta)
  seqname = names(seqf)
  seq = as.character(seqf)
  
  # Load sd info
  bed_colnames = c('seqname','start','end','sdname','similarity','strand')
  sds = read.table(sds_bed); colnames(sds) = bed_colnames
  
  # Load Mutation instructions
  svs = read.table(sv_instr_txt, sep='\t'); colnames(svs) = c('SD', 'SV')
  
  # Make every SV
  for (n_sv in 1:dim(svs)[1]){
    print(n_sv)
    if (svs[n_sv,]$SV == 'INV'){
      seq_sds = carry_out_inv(seq, sds, svs[n_sv,])
      seq = seq_sds[[1]]
      sds = seq_sds[[2]]
    } else if (svs[n_sv,]$SV == 'DEL'){
      seq_sds = carry_out_del(seq, sds, svs[n_sv,])
      seq = seq_sds[[1]]
      sds = seq_sds[[2]]
    } else {
      stop(paste0("Unknown mutation type: ",  svs[n_sv,]$SV))
    }
  }
  
  writeFasta(data.frame(name='sim-sequence-mutated', seq=seq),
             filename=outfasta)
  
  write.table(sds, file=outsd, sep='\t', col.names=F, row.names = F, quote = F)
  
  print('Success! Mutated sequence and sd written.')
}


carry_out_del <- function(seq, sds, sv){
  
  #### A) Delete sequence 
  sds_backup = sds
  # Make helpers for sequence deletion
  sds$sd_middle = sds$start + ((sds$end - sds$start) / 2)
  sds_specific = sds[sds$sdname == sv$SD,]
  
  # Check if SD is possible given SD orientation
  stopifnot("Error: DEL can not be simulated by inversely oriented SDs." = 
              all(sds_specific$strand == '+'))
  
  # Remove sequence
  seq_mutated = paste0(
    substr(seq, 1, sds_specific[1,]$sd_middle),
    substr(seq, sds_specific[2,]$sd_middle, nchar(seq))
  )

  #### B) Update SDs table information
  
  # First, remove the downstream copy of the mediator SD
  # This is not elegantly coded. 
  sds[sds$sdname == sv$SD,] = sds[sds$sdname == sv$SD,][1,]
  sds = unique(sds)
  
  # Which SDs are buried inside the inversion of interest?
  sds$buried_in_inv =   sds$start >= sds_specific[1,'end'] & 
    sds$end <= sds_specific[2,'start']
  
  # ... remove those
  sds = sds[!sds$buried_in_inv,]
  
  # Everything else gets shifted upstream by the amount of deletion. 
  delsize = sds_specific[2,]$end - sds_specific[1,]$end
  
  sds[sds$start > sds_specific[1,]$end]$start = 
    sds[sds$start > sds_specific[1,]$end]$start - delsize
  
  sds[sds$start > sds_specific[1,]$end]$end = 
    sds[sds$start > sds_specific[1,]$end]$end- delsize
  
  return(list(seq_mutated, sds))
}








carry_out_inv <- function(seq, sds, sv){
  
  # TODO: What if a repeat is (partially) contained in another repeat? 
  
  #### A) Invert sequence 
  
  # Make helpers for sequence inversion
  sds$sd_middle = sds$start + ((sds$end - sds$start) / 2)
  sds_specific = sds[sds$sdname == sv$SD,]
  
  # Check if SD is possible given SD orientation
  stopifnot("Error: INV can not be simulated by directly oriented SDs." = 
              all(sds_specific$strand == '-'))
  
  substr(seq, sds_specific[1,]$sd_middle, sds_specific[2,]$sd_middle) = 
    as.character(reverseComplement(DNAString(substr(seq, sds_specific[1,]$sd_middle, sds_specific[2,]$sd_middle))))
  
  #### B) Update SDs table information
  sds_revisited = sds_document_inversion(sds, sv)

  
  return(list(seq, sds_revisited))
}

sds_document_inversion <- function(sds, sv){
  
  # Calculating again here. 
  sds_specific = sds[sds$sdname == sv$SD,]
  
  # Which SDs are buried inside the inversion of interest?
  sds$buried_in_inv =   sds$start >= sds_specific[1,'end'] & 
    sds$end <= sds_specific[2,'start']
  
  # Are there any SD positions to update?
  if (sum(sds$buried_in_inv) > 0){
    
      # A) Update location of affected SDs. #
      # Bad style ('start_new', 'end_new'), should be rewritten. 
      sds$start_new = 0
      sds$end_new = 0
      
      sds[sds$buried_in_inv==T,]$end_new = 
        sds_specific[2,'start'] - (sds[sds$buried_in_inv==T,]$start - sds_specific[1,'end'])
      sds[sds$buried_in_inv==T,]$start_new = 
        sds_specific[2,'start'] - (sds[sds$buried_in_inv==T,]$end - sds_specific[1,'end'])
      
      sds[sds$buried_in_inv==T,]$end = sds[sds$buried_in_inv==T,]$end_new
      sds[sds$buried_in_inv==T,]$start = sds[sds$buried_in_inv==T,]$start_new
      
      sds$end_new = NULL
      sds$start_new = NULL
      
      
      # B) Update orientation state #
      flipped_pairs = 
        unique((sds %>% group_by(sdname) %>% filter(sum(buried_in_inv)==1))$sdname)
      
      # Bad style, could be written better. 
      if (dim(sds[(sds$sdname %in% flipped_pairs) & (sds$strand == '+'),])[1] > 0){
        sds[(sds$sdname %in% flipped_pairs) & (sds$strand == '+'),]$strand = "+flip"
      }
      
      if (dim(sds[(sds$sdname %in% flipped_pairs) & (sds$strand == '-'),])[1] > 0){
        sds[(sds$sdname %in% flipped_pairs) & (sds$strand == '-'),]$strand = "+"
      }
      
      if (dim(sds[(sds$sdname %in% flipped_pairs) & (sds$strand == '+flip'),])[1] > 0){
        sds[sds$strand == '+flip',]$strand = '-'
      }
      
      # Remove helpercolumns
      sds$sd_middle = NULL
      sds$buried_in_inv = NULL
  }
  return(sds)
}
