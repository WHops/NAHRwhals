# Whoeps, 25th Nov 2021
library(Biostrings)

debug = F


carry_out_inv <- function(seq, sds, sv){
  
  # TODO: What if a repeat is (partially) contained in another repeat? 
  
  #### A) Invert sequence 
  
  # Make helpers for sequence inversion
  sds_specific = sds[sds$uid == sv$SD,]
  
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
  sds_specific = sds[sds$uid == sv$SD,]
  
  # Which SDs are buried inside the inversion of interest?
  sds$buried_in_inv = (  
    sds$chromStart >= sds_specific[1,'chromEnd'] & 
      sds$chromEnd <= sds_specific[2,'chromStart']
  )
  
  # Are there any SD positions to update?
  if (sum(sds$buried_in_inv) > 0){
    
    # A.1) Update location of affected SDs. #
    sds[sds$buried_in_inv==T,] = transform(
      sds[sds$buried_in_inv==T,],
      chromEnd = sds_specific[2,'chromStart'] - (chromStart - sds_specific[1,'chromEnd']),
      chromStart = sds_specific[2,'chromStart'] - (chromEnd - sds_specific[1,'chromEnd'])
    )
    
    # A.2) Everyone whose partner has moved. This is not well written...
    for (transformed_uid in unique(sds[sds$buried_in_inv==T,])){
      sds_ui = sds[sds$uid == transformed_uid,]
      flipped_id = which(sds_ui$buried_in_inv)
      partner_of_flipped_id = 3 - flipped_id # this transforms 2 -> 1 and 1 -> 2
      
      sds[sds$uid == transformed_uid,][partner_of_flipped_id,] = 
        transform(
          sds[sds$uid == transformed_uid,][partner_of_flipped_id,],
          otherStart = sds[sds$uid == transformed_uid,][flipped_id,'chromStart'],
          otherEnd = sds[sds$uid == transformed_uid,][flipped_id,'chromEnd']
        )
      
    }
    
    # B) Update orientation state #
    flipped_pairs = 
      unique((sds %>% group_by(uid) %>% filter(sum(buried_in_inv)==1))$uid)
    
    sds[(sds$uid %in% flipped_pairs),] = transform(
      sds[(sds$uid %in% flipped_pairs),],
      strand = ifelse(strand == '+', '-', '+')
    )
    
    # Remove helpercolumns
    sds$sd_middle = NULL
    sds$buried_in_inv = NULL
  }
  return(sds)
}




carry_out_del <- function(seq, sds, sv){
  
  #### A) Delete sequence 

  # Make helpers for sequence deletion
  sds_specific = sds[sds$uid == sv$SD,]
  
  # Check if SD is possible given SD orientation
  stopifnot("Error: DEL can not be simulated by inversely oriented SDs." = 
              all(sds_specific$strand == '+'))
  
  # Remove sequence
  seq_mutated = paste0(
    substr(seq, 1, sds_specific[1,]$sd_middle),
    substr(seq, sds_specific[2,]$sd_middle + 1, nchar(seq))
  )

  #### B) Update SDs table information
  
  sds_revisited = sds_document_deletion(sds, sv)
  
  return(list(seq_mutated, sds))
}
sds_document_deletion <- function(sds, sv){
  
  # First, remove the downstream copy of the mediator SD
  # This is not elegantly coded. 
  sds[sds$uid == sv$SD,] = sds[sds$uid == sv$SD,][1,]
  sds = unique(sds)
  
  # Which SDs are buried inside the inversion of interest?
  sds$buried_in_del = (  
    sds$chromStart >= sds_specific[1,'chromEnd'] & 
    sds$chromEnd <= sds_specific[2,'chromStart']
  )
  del_pairs = unique(sds[sds$buried_in_del == T,]$uid)
  
  # ... remove those
  sds = sds[!sds$uid %in% del_pairs,]
  
  # Everything else gets shifted upstream by the amount of deletion. 
  delsize = sds_specific[2,]$chromEnd - sds_specific[1,]$chromEnd
  
  sds[sds$chromStart > sds_specific[1,]$chromEnd] = transform(
    sds[sds$chromStart > sds_specific[1,]$chromEnd],
    chromStart = chromStart - delsize,
    chromEnd   = chromEnd - delsize
  )
  
  sds[sds$otherStart > sds_specific[1,]$otherEnd] = transform(
    sds[sds$otherStart > sds_specific[1,]$otherEnd],
    otherStart = otherStart - delsize,
    otherEnd   = otherEnd - delsize
  )
  
  # Remove helpercolumns
  sds$sd_middle = NULL
  sds$buried_in_inv = NULL
  
  return(sds)
}


carry_out_dup <- function(seq, sds, sv){
  
  
  #### A) Duplicate sequence 
  
  # Make helpers for sequence duplication
  sds_specific = sds[sds$uid == sv$SD,]
  
  # Check if SV is possible given SD orientation
  stopifnot("Error: DUP can not be simulated by inversely oriented SDs." = 
              all(sds_specific$strand == '+'))
  
  # Add sequence
  seq_mutated = paste0(
    substr(seq, 1, sds_specific[2,]$sd_middle), # Seq up to upstream bp
    substr(seq, sds_specific[1,]$sd_middle + 1, sds_specific[2,]$sd_middle),
    substr(seq, sds_specific[2,]$sd_middle + 1, nchar(seq))
  )
  
  #### B) Update SDs table information
  
  sds_revisited = sds_document_duplication(sds, sv)
  
  return(list(seq_mutated, sds))
}
sds_document_duplication <- function(sds,sv){
  
  # Everything else gets shifted downstream by the amount of duplication 
  dupsize = sds_specific[2,]$chromEnd - sds_specific[1,]$chromEnd
  
  orig = sds[sds$uid == svs[1,'SD'],]
  copy1 = sds[sds$uid == svs[1,'SD'],]

  copy1[1,] = transform(
    copy1[1,],
    otherStart = otherStart + dupsize,
    otherEnd = otherEnd + dupsize
  )
  copy1[2,] = transform(
    copy1[2,],
    chromStart = chromStart + dupsize,
    chromEnd = chromEnd + dupsize
  )
  
  copy2 = copy1
  
  copy2[1,] = transform(
    copy2[1,],
    chromStart = chromStart + dupsize,
    chromEnd = chromEnd + dupsize

  )
  copy2[2,] = transform(
    copy2[2,],
    otherStart = otherStart + dupsize,
    otherEnd = otherEnd + dupsize
  )
  
  copy1$uid = paste0(orig$uid, 'a')
  copy2$uid = paste0(orig$uid, 'b')
  
  all = rbind(rbind(orig, copy1), copy2)
  
  #sds[sds$uid == sv$SD,] = sds[sds$uid == sv$SD,][1,]
  #sds = unique(sds)
  
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









mutate_seq <- function(seq_fasta, sds_tsv, sv_instr_txt, outfasta, outsd){
  
  if (debug){
    seq_fasta = "../res/fa/s10.fa"
    sds_tsv = "../data/sds10y.tsv"
    sv_instr_txt = "../data/svs10.txt"
  }
  
  # Load both sequence
  seqf = readDNAStringSet(seq_fasta)
  seqname = names(seqf)
  seq = as.character(seqf)
  
  # Load sd info
  sds_colnames = c(
    "chrom", "chromStart", "chromEnd", "name",
    "score","strand", "otherChrom","otherStart", "otherEnd",
    "otherSize", "uid", "posBasesHit", "testResult",
    "verdict", "chits","ccov","alignfile","alignL","indelN",
    "indelS","alignB","matchB","mismatchB", "transitionsB",
    "transversionsB","fracMatch","fracMatchIndel","jcK", "k2K"
  )
  sds = read.table(sds_tsv); colnames(sds) = sds_colnames
  
  # Load Mutation instructions
  svs = read.table(sv_instr_txt, sep='\t'); colnames(svs) = c('SD', 'SV')
  
  # We assume the breakpoint is always in the middle of the repeat. 
  # If we want the breakpoints different at some point, this here is
  # the place to change this. 
  sds$sd_middle = sds$chromStart + ((sds$chromEnd - sds$chromStart) / 2)
  
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
