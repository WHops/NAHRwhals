# Fasta extraction and asssembly tools

#' extract_subseq
#' helperfunction to extract a subfasta by coordinates
#'
#' Is awfully slow unfortunately. We should find a way to only read
#' the chromosome of interest.
#' @author Wolfram Hoeps
#' @export
extract_subseq <- function(infasta, seqname, start, end, outfasta) {
  # Read the whole infasta
  seq = Biostrings::readDNAStringSet(infasta)
  
  # Subset it to the sequence name of interest
  if (!is.null(seqname)) {
    subseq = Biostrings::subseq(seq[[seqname]], start = start, end = end)
  } else {
    subseq = Biostrings::subseq(seq, start = start, end = end)
  }
  # Write fasta. Imported function from seqbilder_functions.R
  writeFasta(data.frame(name = 'seq', seq = as.character(subseq)), outfasta)
  
}



# To document
#' filter_paf_to_main_region
#' Script to filter down a whole-genome paf to the 'best matching' region.
#' query should be short sequence fragments (1kb, 10kb, ...) of the roi.

#' @author Wolfram Hoeps
#' @export
filter_paf_to_main_region <- function(paflink, outpaflink) {
  # Read in paf
  paf = read.table(paflink,
                   sep = '\t',
                   fill = T,
                   row.names = NULL)
  
  colnames_paf = c(
    'qname',
    'qlen',
    'qstart',
    'qend',
    'strand',
    'tname',
    'tlen',
    'tstart',
    'tend',
    'nmatch',
    'alen',
    'mapq'
  )
  colnames(paf)[1:length(colnames_paf)] = colnames_paf
  
  # Sum of lengths of the sniplets. Group by sniplets (qname), take one per group
  # sum their qlen.
  insequence_len = sum(na.omit(as.numeric(
    dplyr::slice(dplyr::group_by(paf, qname), 1)$qlen
  )))
  
  # For every snipplet, find the alignment with the largest number of matches.
  paf_winners = dplyr::filter(dplyr::group_by(paf, qname), nmatch == max(nmatch))
  
  # Of these winning alignments, which is the most frequent target contig/sequence/chr?
  chr_winners = tail(names(sort(table(
    paf_winners$tname
  ))), 1)
  # And those snipplets matching to the winner chromosome - where do they fall (median)?
  middle_median = median(((
    paf_winners[paf_winners$tname == chr_winners,]$tend +
      paf_winners[paf_winners$tname == chr_winners,]$tstart
  ) / 2))
  
  # Now simply extend from median towards front and back.
  start_winners = middle_median - (1.5 * (insequence_len / 2))
  end_winners = middle_median + (1.5 * (insequence_len / 2))
  
  # Here we cut, but I am not sure if this is really the way to go?
  paf_cut = paf[((paf$tname == chr_winners) &
                   (paf$tstart >= start_winners) &
                   (paf$tend <= end_winners)),]
  
  
  aln_containing = cpaf_i[(cpaf_i$qstart <= start) &
                            (cpaf_i$qend >= end),]
  
  
  
  write.table(
    paf_cut,
    file = outpaflink,
    sep = '\t',
    row.names = F,
    col.names = F,
    quote = F
  )
  
  
}

#' @export
flip_query_target <- function(inpaf, outpaf) {
  # Read in paf
  paf = read.table(inpaf,
                   sep = '\t',
                   fill = T,
                   row.names = NULL)
  
  colnames_paf = c(
    'qname',
    'qlen',
    'qstart',
    'qend',
    'strand',
    'tname',
    'tlen',
    'tstart',
    'tend',
    'nmatch',
    'alen',
    'mapq'
  )
  colnames(paf)[1:length(colnames_paf)] = colnames_paf
  
  # Transform
  paf_transform = transform(
    paf,
    qname = tname,
    tname = qname,
    qlen = tlen,
    tlen = qlen,
    qstart = tstart,
    tstart = qstart,
    qend = tend,
    tend = qend
  )
  
  write.table(
    paf_transform,
    file = outpaf,
    sep = '\t',
    row.names = F,
    col.names = F,
    quote = F
  )
  
}



# TRASH #

#' extract_subseq_bedtools

#' @author Wolfram Hoeps
#' @export
extract_subseq_bedtools <-
  function(infasta, seqname, start, end, outfasta) {
    # Needed for parallel runs
    
    bedtoolsloc = query_config("bedtools")
    random_tag = as.character(runif(1, 1e10, 1e11))
    tmp_bedfile = paste0('region2_', random_tag, '.bed')
    region = paste0(seqname,
                    "\t",
                    format(start, scientific = F),
                    "\t",
                    format(end, scientific = F))
    system(paste0('echo "', region, '" > ', tmp_bedfile))
    system(
      paste0(
        bedtoolsloc,
        " getfasta -fi ",
        infasta,
        " -bed ",
        tmp_bedfile,
        " > ",
        outfasta
      )
    )
    system(paste0('rm ', tmp_bedfile))
    
    print(paste0('Subsequence extracted and saved to ', outfasta))
    
  }



#' @export
find_punctual_liftover <- function(cpaf, pointcoordinate, chrom) {
  # Find all alignments overlapping with the pointcoordinate
  overlappers = cpaf[((cpaf$qname == chrom) &
                        (cpaf$qstart <= pointcoordinate) &
                        (cpaf$qend >= pointcoordinate)
  ),]
  
  # If no alignment overlaps, we return an empty value.
  if (dim(overlappers)[1] == 0) {
    return(NULL)
  }
  
  # Find the alignment with the hightest nmatch value.
  best_aln = overlappers[overlappers$nmatch == max(overlappers$nmatch),]
  
  # Liftover the point with that highest alignment.
  if (best_aln$strand == '+') {
    coord = best_aln$tstart + (pointcoordinate - best_aln$qstart)
  } else {
    coord = best_aln$tend - (pointcoordinate - best_aln$qstart)
  }
  
  return(data.frame(seqname = best_aln$tname, liftover_coord = coord))
}


#' @export
liftover_coarse <-
  function(seqname,
           start,
           end,
           conversionpaf_link,
           n_probes = 100,
           lenfactor = 1.5) {
    # Load conversionpaf
    cpaf = read.table(
      conversionpaf_link,
      sep = '\t',
      fill = T,
      row.names = NULL
    )
    
    colnames_paf = c(
      'qname',
      'qlen',
      'qstart',
      'qend',
      'strand',
      'tname',
      'tlen',
      'tstart',
      'tend',
      'nmatch',
      'alen',
      'mapq'
    )
    colnames(cpaf)[1:length(colnames_paf)] = colnames_paf
    
    # We take n_probes single points from the alignment, see where they fall, and take their median
    # as the median of our alignment.
    
    # Define a vector of probe positions (equally spaced between start and end)
    pos_probes = as.integer(start + (((end - start) / (n_probes - 1)) * (0:(n_probes -
                                                                              1))))
    
    # Liftover every probe
    liftover_coords <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(liftover_coords) <- c('seqname', 'liftover_coord')
    for (pointcoord in pos_probes) {
      liftover_coords = rbind(liftover_coords,
                              find_punctual_liftover(cpaf, pointcoord, seqname))
    }
    
    # What is the majority vote for the target chromosome?
    winner_chr = names(sort(table(liftover_coords$seqname), decreasing = TRUE)[1])
    
    # And those snipplets matching to the winner chromosome - where do they fall (median)?
    middle_median = median(liftover_coords[liftover_coords$seqname == winner_chr,]$liftover_coord)
    
    # Now  extend from median towards front and back.
    insequence_len = end - start
    start_winners = middle_median - (lenfactor * (insequence_len / 2))
    end_winners =   middle_median +   (lenfactor * (insequence_len / 2))
    
    # Make sure we don't exceed chromosome boundaries
    start_winners_cutoff = as.integer(max(0, start_winners))
    end_winners_cutoff =   as.integer(min(cpaf[cpaf$tname == winner_chr,][1, 'tlen'], end_winners))
    
    if ((start_winners < 0) |
        (end_winners > cpaf[cpaf$tname == winner_chr,][1, 'tlen'])) {
      print('Warning! Reaching the end of alignment!')
    }
    
    return(
      list(
        lift_contig = winner_chr,
        lift_start = start_winners_cutoff,
        lift_end = end_winners_cutoff
      )
    )
    
  }
