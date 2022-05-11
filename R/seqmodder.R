# Whoeps, 25th Nov 2021

#' Introduce mutations to a DNA sequence. 
#' 
#' This is sitting at the core of the nahrchainer module. 
#' 
#' @param seq_fasta [character/link] A single-sequence fasta file containing the
#' original sequence, which is to be modified. Should contain SDs. 
#' @param sds_tsv [character/link] A tsv file DESCRIBING the seq_fasta. This information will
#' be used to update presumed SD position post modification (not working well yet).
#' @param sv_instr_txt [character/link] A tab-separated file with colums [SD-uid, SV].
#' SD-uid has to be existing in sds_tsv. SV can be [INV, DEL, DUP]. s
#' @param outfasta [character/link] Outputfile for the modified fasta
#' @param outsd [character/link] Outputfile for the modified SD tsv file. 
#' @return nothing. But output files written. 
#' 
#' @author Wolfram Höps
#' @export
mutate_seq <- function(seq_fasta, sds_tsv, sv_instr_txt, outfasta, outsd, debug=F){
  
  debug=F
  if (debug){
    seq_fasta = "../res/fa/invs.fa"
    sds_tsv = "../data/sds10y.tsv"
    sv_instr_txt = "../data/svs10.txt"
  }
  
  # Load both sequence
  seqf = Biostrings::readDNAStringSet(seq_fasta)
  seqname = names(seqf) #should be Biostrings too?
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
  
  
  # Make every SV
  for (n_sv in 1:dim(svs)[1]){
    if (svs[n_sv,]$SV == 'INV'){
      seq_sds = carry_out_inv(seq, sds, svs[n_sv,])
      seq = seq_sds[[1]]
      sds = seq_sds[[2]]
    } else if (svs[n_sv,]$SV == 'DEL'){
      seq_sds = carry_out_del(seq, sds, svs[n_sv,])
      seq = seq_sds[[1]]
      sds = seq_sds[[2]]
    } else if (svs[n_sv,]$SV == 'DUP'){
      seq_sds = carry_out_dup(seq, sds, svs[n_sv,])
      seq = seq_sds[[1]]
      sds = seq_sds[[2]]
    } else {
      stop(paste0("Unknown mutation type: ",  svs[n_sv,]$SV))
    }
  }
  
  # Write output files
  writeFasta(data.frame(name='sim-sequence-mutated', seq=seq),
             filename=outfasta)
  
  write.table(sds, file=outsd, sep='\t', col.names=F, row.names = F, quote = F)
  
  print('Success! Mutated sequence and sd written.')
}


#' Introduce an INV to a DNA sequence
#' 
#' @description This function is typically called by the mutate_seq function
#' 
#' @param seq [character] DNA sequence as character vector
#' @param sds [data.frame] A dataframe with tsv SD information
#' @param sv [data.frame] A data framewith columns [SD-uid, SV].
#' SD-uid has to be existing in sds_tsv. SV can be [INV, DEL, DUP]
#' @return a list object containing the mutated sequence and sds dataframe. 
#' 
#' @author Wolfram Höps
#' @export
carry_out_inv <- function(seq, sds, sv){
  
  # We assume the breakpoint is always in the middle of the repeat. 
  # If we want the breakpoints different at some point, this here is
  # the place to change this. 
  sds$sd_middle = sds$chromStart + ((sds$chromEnd - sds$chromStart) / 2)
  
  # TODO: What if a repeat is (partially) contained in another repeat? 
  
  #### A) Invert sequence 
  
  # Make helpers for sequence inversion
  sds_specific = sds[sds$uid == sv$SD,]
  
  # Check if SD is possible given SD orientation
  stopifnot("Error: INV can not be simulated by directly oriented SDs." = 
              all(sds_specific$strand == '-'))
  test = Biostrings::DNAString(substr(seq, sds_specific[1,]$sd_middle, sds_specific[2,]$sd_middle))
  
  
  substr(seq, sds_specific[1,]$sd_middle, sds_specific[2,]$sd_middle) = 
    as.character(Biostrings::reverseComplement(Biostrings::DNAString(substr(seq, sds_specific[1,]$sd_middle, sds_specific[2,]$sd_middle))))
  
  #### B) Update SDs table information
  sds_revisited = sds_document_inversion(sds, sv)
  
  # W, 14th March 2022. We need to sort here. 
  # So that '1,end' and '2,start' are what they
  # are supposed to be. 
  sds_revisited = sds_revisited[
    order(sds_revisited$uid, sds_revisited$chromStart),
  ]
  
  return(list(seq, sds_revisited))
}


#' Document the introduction of an INV in the SD tsv dataframe. 
#' 
#' @description This function is typically called by the carry_out_inv function.
#' ALSO the function is highly experimental and probably not working.
#' We recommend not rely on it at this point. 
#' 
#' @param sds [data.frame] A dataframe with tsv SD information
#' @param sv [data.frame] A data framewith columns [SD-uid, SV].
#' SD-uid has to be existing in sds_tsv. SV can be [INV, DEL, DUP]
#' @returns a mutated sds dataframe. 
#' 
#' @author Wolfram Höps
#' @export
sds_document_inversion <- function(sds, sv){
  
  # W, 14th March 2022. We need to sort here. 
  # So that '1,end' and '2,start' are what they
  # are supposed to be. 
  sds = sds[
    order(sds$uid, sds$chromStart),
  ]
  
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
    for (transformed_uid in unique(sds[sds$buried_in_inv==T,'uid'])){
      
      
      sds_ui = sds[sds$uid == transformed_uid,]
      flipped_id = which(sds_ui$buried_in_inv)
      if (length(flipped_id) == 1){
      partner_of_flipped_id = 3 - flipped_id # this transforms 2 -> 1 and 1 -> 2
      
      sds[sds$uid == transformed_uid,][partner_of_flipped_id,] = 
        transform(
          sds[sds$uid == transformed_uid,][partner_of_flipped_id,],
          otherStart = sds[sds$uid == transformed_uid,][flipped_id,'chromStart'],
          otherEnd = sds[sds$uid == transformed_uid,][flipped_id,'chromEnd']
        )
      }
      
    }
    
    # B) Update orientation state #
    flipped_pairs = 
      unique((sds %>% group_by(uid) %>% filter(sum(buried_in_inv)==1))$uid)
    
    sds[(sds$uid %in% flipped_pairs),] = transform(
      sds[(sds$uid %in% flipped_pairs),],
      strand = ifelse(strand == '+', '-', '+')
    )
    
    sds$buried_in_inv = NULL
  }
  
  sds$sd_middle = NULL
  return(sds)
}

#' Introduce a DEL to a DNA sequence
#' 
#' @description This function is typically called by the mutate_seq function
#' 
#' 
#' @param seq [character] DNA sequence as character vector
#' @param sds [data.frame] A dataframe with tsv SD information
#' @param sv [data.frame] A data framewith columns [SD-uid, SV].
#' SD-uid has to be existing in sds_tsv. SV can be [INV, DEL, DUP]
#' @return a list object containing the mutated sequence and sds dataframe. 
#' 
#' @author Wolfram Höps
#' @export
carry_out_del <- function(seq, sds, sv){
  
  # We assume the breakpoint is always in the middle of the repeat. 
  # If we want the breakpoints different at some point, this here is
  # the place to change this. 
  sds$sd_middle = sds$chromStart + ((sds$chromEnd - sds$chromStart) / 2)
  
  #### A) Delete sequence 
  #sv = svs[n_sv,]
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
  
  # W, 14th March 2022. We need to sort here. 
  # So that '1,end' and '2,start' are what they
  # are supposed to be. 
  sds_revisited = sds_revisited[
    order(sds_revisited$uid, sds_revisited$chromStart),
  ]
  
  return(list(seq_mutated, sds_revisited))
}

#' Document the introduction of a DEL in the SD tsv dataframe. 
#' 
#' @description This function is typically called by the carry_out_del function.
#' ALSO the function is highly experimental and probably not working.
#' We recommend not rely on it at this point. 
#' 
#' @param sds [data.frame] A dataframe with tsv SD information
#' @param sv [data.frame] A data framewith columns [SD-uid, SV].
#' SD-uid has to be existing in sds_tsv. SV can be [INV, DEL, DUP]
#' @returns a mutated sds dataframe. 
#' 
#' @author Wolfram Höps
#' @export
sds_document_deletion <- function(sds, sv){
  
  sds_specific = sds[sds$uid == sv$SD,]
  
  # Remove the mediator pair
  sds = sds[!sds$uid == sv$SD,]
  
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
  
  sds[sds$chromStart > sds_specific[1,]$chromEnd,] = transform(
    sds[sds$chromStart > sds_specific[1,]$chromEnd,],
    chromStart = chromStart - delsize,
    chromEnd   = chromEnd - delsize
  )
  
  sds[sds$otherStart > sds_specific[1,]$otherEnd,] = transform(
    sds[sds$otherStart > sds_specific[1,]$otherEnd,],
    otherStart = otherStart - delsize,
    otherEnd   = otherEnd - delsize
  )
  
  # Remove helpercolumns
  sds$buried_in_inv = NULL
  
  return(sds)
}

#' Helperfunction 1/2 for carry_out_dup
#' 
#' @author Wolfram Höps
#' @export
expand_one_of_pair <- function(sds, uid, n_of_pair, dupsize){
  
  # We put in the whole sds first. 
  
  
  orig = sds[sds$uid == uid,]
  copy1 = sds[sds$uid == uid,]
  
  if (n_of_pair == 1){
    # 1-1.5
    copy1[1,] = transform(
      copy1[1,],
      otherStart = chromStart + dupsize,
      otherEnd = chromEnd + dupsize
    )
    
    # 1.5-1
    copy1[2,] = transform(
      copy1[2,],
      chromStart = otherStart + dupsize,
      chromEnd = otherEnd + dupsize
    )
    
    copy2 = orig
    
    # 1.5-2
    copy2[1,] = transform(
      copy2[1,],
      chromStart = chromStart + dupsize,
      chromEnd = chromEnd + dupsize
    )
    
    # 2-1.5
    copy2[2,] = transform(
      copy2[2,],
      otherStart = otherStart + dupsize,
      otherEnd = otherEnd + dupsize
    )
  } else if (n_of_pair == 2){
    # 1-3
    copy1[1,] = transform(
      copy1[1,],
      otherStart = otherStart + dupsize,
      otherEnd = otherEnd + dupsize
    )
    
    # 3-1
    copy1[2,] = transform(
      copy1[2,],
      chromStart = chromStart + dupsize,
      chromEnd = chromEnd + dupsize
    )
    
    copy2 = copy1
    
    # 2-3
    copy2[1,] = transform(
      copy2[1,],
      chromStart = otherStart - dupsize,
      chromEnd = otherEnd - dupsize
      
    )
    
    # 3-2
    copy2[2,] = transform(
      copy2[2,],
      otherStart = chromStart - dupsize,
      otherEnd = chromEnd - dupsize
    )
  }
  copy1$uid = paste0(orig$uid, 'a', n_of_pair)
  copy2$uid = paste0(orig$uid, 'b', n_of_pair)
  
  all = rbind(rbind(sds, copy1), copy2)
  
  return(all)
}

#' Helperfunction 2/2 for carry_out_dup
#' 
#' @author Wolfram Höps
#' @export
expand_whole_pair <- function(sds, uid, dupsize){
  
  sds = sds_bu
  sds1 = expand_one_of_pair(sds, uid, 1, dupsize)
  sds2 = expand_one_of_pair(sds1, uid, 2, dupsize)
  
  # Add the last connection (between the two new SDs)
  # Manually. 
  to_cp = sds2[sds2$uid == uid,]
  to_cp[,c('chromStart', 'chromEnd', 'otherStart', 'otherEnd')] = 
    to_cp[,c('chromStart', 'chromEnd', 'otherStart', 'otherEnd')] + dupsize
  to_cp$uid = paste0(to_co$uid, '3')
  sds3 = rbind(sds2, to_cp)
  
  sds3 = sds3[order(sds3$chromStart),]
  
  return(sds3)
}

#' Introduce a DUP to a DNA sequence
#' 
#' @description This function is typically called by the mutate_seq function
#' 
#' @param seq [character] DNA sequence as character vector
#' @param sds [data.frame] A dataframe with tsv SD information
#' @param sv [data.frame] A data framewith columns [SD-uid, SV].
#' SD-uid has to be existing in sds_tsv. SV can be [INV, DEL, DUP]
#' @return a list object containing the mutated sequence and sds dataframe. 
#' 
#' @author Wolfram Höps
#' @export
carry_out_dup <- function(seq, sds, sv){
  
  # We assume the breakpoint is always in the middle of the repeat. 
  # If we want the breakpoints different at some point, this here is
  # the place to change this. 
  sds$sd_middle = sds$chromStart + ((sds$chromEnd - sds$chromStart) / 2)
  
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
  
  # W, 16th Feb 2022:
  # [Something is wrong in this functin. Commenting it out for now.]
  sds_revisited = sds# sds_document_duplication(sds, sv)
  
  # W, 14th March 2022. We need to sort here. 
  # So that '1,end' and '2,start' are what they
  # are supposed to be. 
  sds_revisited = sds_revisited[
    order(sds_revisited$uid, sds_revisited$chromStart),
  ]
  
  return(list(seq_mutated, sds_revisited))
}


#' Document the introduction of a DUP in the SD tsv dataframe. 
#' 
#' @description This function is typically called by the carry_out_dup function.
#' ALSO the function is highly experimental and probably not working.
#' We recommend not rely on it at this point. 
#' 
#' @param sds [data.frame] A dataframe with tsv SD information
#' @param sv [data.frame] A data framewith columns [SD-uid, SV].
#' SD-uid has to be existing in sds_tsv. SV can be [INV, DEL, DUP]
#' @returns a mutated sds dataframe. 
#' 
#' @author Wolfram Höps
#' @export
sds_document_duplication <- function(sds,sv){
  
  sds_specific = sds[sds$uid == sv$SD,]
  #sds_bu = sds
  #sv_bu = sv
  sds = sds[sds$uid %in% c("SDA", "SDB", "SDC"),]

  # Remember this for later
  uid_orig = unique(sds$uid)
  
  dupsize = sds_specific[2,]$chromEnd - sds_specific[1,]$chromEnd
  
  sds = expand_one_of_pair(sds, "SDC", 2,  dupsize)
  
  # Which SDs are buried inside the duplication of interest?
  sds$buried_in_dup =   sds$chromStart >= sds_specific[1,'chromEnd'] & 
    sds$chromEnd <= sds_specific[2,'chromStart']
   
  pairs_one_duplicated = 
    unique((sds %>% group_by(uid) %>% filter(sum(buried_in_dup)==1))$uid)
  
  pairs_two_duplicated = 
    unique((sds %>% group_by(uid) %>% filter(sum(buried_in_dup)==2))$uid)
  
  for (uid in pairs_one_duplicated){
    sds = expand_one_of_pair(sds, uid, 2,  dupsize)
  }
  
  for (uid in pairs_two_duplicated){
    sds = expand_whole_pair(sds, uid, dupsize)
  }
  
  
  
  # Everything else gets shifted upstream by the amount of deletion. 

  # Chrom starts, ends
  sds_predup = sds[sds$uid %in% uid_orig,]
  
  
  sds_predup[sds_predup$chromStart > sds_specific[1,]$otherEnd,] = transform(
    sds_predup[sds_predup$chromStart > sds_specific[1,]$otherEnd,],
    chromStart = chromStart + dupsize,
    chromEnd   = chromEnd + dupsize
  )
  
  sds_predup[sds_predup$otherStart > sds_specific[1,]$otherEnd,] = transform(
    sds_predup[sds_predup$otherStart > sds_specific[1,]$otherEnd,],
    otherStart = otherStart + dupsize,
    otherEnd   = otherEnd + dupsize
  )
  
  # Adjust starts and ends of everything downstream
  sds_predup[sds_predup$chromStart > (sds_specific[1,]$otherEnd),]$chromStart = 
    sds_predup[sds_predup$chromStart > sds_specific[1,]$otherEnd,]$chromStart + dupsize
  
  sds_predup[sds_predup$chromEnd > sds_specific[1,]$otherEnd,]$chromEnd = 
    sds_predup[sds_predup$chromEnd > sds_specific[1,]$otherEnd,]$chromEnd + dupsize
  
  
  # Add the duplication itself 
  
  # print(sds_specific)
  sds_dup = sds_specific
  sds_dup[1,] = transform(
    sds_dup[1,],
    chromEnd = otherEnd,
    otherEnd = otherEnd + dupsize
  )
  sds_dup[2,] = transform(
    sds_dup[2,],
    chromEnd = chromEnd + dupsize,
    otherEnd = chromEnd
  )
  sds_dup$uid = paste0(sds_dup$uid, 'dup')
  sds_dup$buried_in_dup = F
  # print('@@########################')
  # print(sds_dup)
  # print('@@########################')
  # print(sds_predup)
  # Re-merge
  sds_return = rbind(rbind(sds_predup, sds[!sds$uid %in% uid_orig,]), sds_dup)
  #sds_return = rbind(sds_predup, sds[!sds$uid %in% uid_orig,])
  
  sds_return = sds_return[order(sds_return$uid),]
  
  return(sds_return)
}


#' carry_out_compressed_sv
#'
#' @description Carry out an SV IN GRIDSPACE
#'
#' @param bitl  an nxm matrix (n=gridpoints_x, m=gridpoints_y)
#' @param input_ins a vector with instruction: p1, p2, sv.
#'
#' @author Wolfram Höps
#' @export
carry_out_compressed_sv <- function(bitl, input_ins) {
  
  pair = as.numeric(input_ins[1:2])
  action = input_ins[3]
  
  # Modify the matrix accordingly
  if (action == 'del') {
    bitl_mut = bitl[, -c(as.numeric(pair[1]):(as.numeric(pair[2]) - 1))]
  } else if (action == 'dup') {
    
    bitl_mut = cbind(bitl[, 1:as.numeric(pair[2])],
                     bitl[, as.numeric(pair[1] + 1):dim(bitl)[2]])
    
    if (pair[2] - pair[1] == 1){
      if (pair[1] == 1){
        # W, 10th March 2022, Changed this a bit to account for cases
        # with only 2 columns. 
        colnames(bitl_mut)[dim(bitl_mut)[2]] = colnames(bitl_mut)[2] 
      } else if (pair[2] == dim(bitl)[2]){
        colnames(bitl_mut)[dim(bitl_mut)[2]] = colnames(bitl_mut)[dim(bitl_mut)[2]-1]
      }
    }
  } else if (action == 'inv') {
    
    # W, 8th March 2022. 
    # Noticed a problem with one-column inversions. Here, the expression 
    # -bitl[, (as.numeric(pair[2]) - 1):(as.numeric(pair[1]) + 1)])
    # Returns a one-line dataframe, which loses its column name. 
    # Therefore, we now handle one-column inversions differently now. 
    if ((pair[2] - pair[1]) == 1){
      bitl_mut = bitl
    } else if (pair[2] - pair[1] > 1){
      
        # A lot of bordercase handling...
        if ((pair[1] == 1) & (pair[2] != dim(bitl)[2])){
          bitl_mut = cbind( -bitl[, (as.numeric(pair[2])):(as.numeric(pair[1]))],
                             bitl[, as.numeric(pair[2]+1):dim(bitl)[2]])
          # Manually edit the last entry of the header. This is necessary if pair[2] == dim-1.
          # This is in response to an error ('error #2, may11th 2022')
          colnames(bitl_mut)[dim(bitl_mut)[2]] = colnames(bitl)[dim(bitl_mut)[2]]
          
          # This black was added not for an error, but i realized that I think it should be in (may11th 2022)
        } else if ((pair[1] > 1) & (pair[2] == dim(bitl)[2])){
          bitl_mut = cbind( bitl[, 1:(as.numeric(pair[1])-1)],
                            -bitl[, (as.numeric(pair[2])):(as.numeric(pair[1]))])
          # Manually edit the last entry of the header. This is necessary if pair[1] == 2
          # This is in response to an error ('error #2, may11th 2022')
          colnames(bitl_mut)[1] = colnames(bitl)[1]
          
        } else if ((pair[1] == 1) & (pair[2] == dim(bitl)[2])){
          bitl_mut = -bitl[, (as.numeric(pair[2])):(as.numeric(pair[1]))]
          
        } else {
        bitl_mut = cbind(cbind(bitl[, 1:as.numeric(pair[1])],
                               -bitl[, (as.numeric(pair[2]) - 1):(as.numeric(pair[1]) + 1)]),
                         bitl[, as.numeric(pair[2]):dim(bitl)[2]])
        } 
        # Border case handling. If inversion is exactly one column, 
        # (i.e. pair 3-5), then we assign the colnames manually. 
        if (pair[2] - pair[1] == 2){
          colnames(bitl_mut) = colnames(bitl)
        } else if (pair[2] == dim(bitl)[2]) {
          colnames(bitl_mut)[pair[2]] = colnames(bitl_mut)[pair[1]]
        }
    }
    } else if (action == 'ref'){
      bitl_mut = bitl
    }
  
  # It can happen that only one column is left. In that case...
  
  if (any((class(bitl_mut) != c('matrix', 'array')))){
    bitl_mut = as.matrix(bitl_mut)
    row.names(bitl_mut) = row.names(bitl)
    colnames(bitl_mut) = colnames(bitl)[1]
  }
  
  # if ("" %in% colnames(bitl_mut)){
  #   browser()
  # }
  # W, 7th March 2022. Excluding this line since bitl now has 
  # meaningful rownames and colnames. 
  #colnames(bitl_mut) = 1:dim(bitl_mut)[2]

  # Turn minus-zero into zero.
  bitl_mut[bitl_mut == 0] = 0
  return(bitl_mut)
  
}

#' modify_gridmatrix
#' 
#' Introduce a chain of mutations to gridmatrix. 
#' Is used after mutation simulation to verify/visualize
#' the result. 
#' @param gmx  GridMatriX. Row/Colnames don't matter. 
#' @param r1 mutation instruction. Not sure about format right now. TODO. 
#' @export
modify_gridmatrix <- function(gmx, r1){
  
  # If there is no mutation to be done, return matrix as it is. 
  if (r1[1,'mut1'] == 'ref'){
    return(gmx)
  }
  
  # How many mutations do we have to run?
  nmut = length(r1[, colSums(is.na(r1)) == 0]) - 1
  
  # Run each mutation. 
  for (i in 1:nmut){
    instr = r1[,i+1]
    start = as.numeric(strsplit(instr, '_')[[1]][1])
    end = as.numeric(strsplit(instr, '_')[[1]][2])
    action = strsplit(instr, '_')[[1]][3]
    
    # Run ONE mutation
    gmx = carry_out_compressed_sv(gmx, c(start, end, action))
  }
  
  # Not ENTIRELY sure, but for some reason we have to transpose here to get the 
  # correct result back apparently? Might be wrong. Keep an eye on it. 
  gmx = t(gmx)
  
  # Simple row and colnames for easy plotting. 
  colnames(gmx) = 1:ncol(gmx)
  row.names(gmx) = 1:nrow(gmx)
  
  return(gmx)
}


