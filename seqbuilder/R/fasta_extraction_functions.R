# Fasta extraction and asssembly tools

#' extract_subseq
#' helperfunction to extract a subfasta by coordinates
#'
#' Is awfully slow unfortunately. We should find a way to only read
#' the chromosome of interest.
#' @author Wolfram Hoeps
#' @rdname alignment
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


#' @description A core function. Give me a fasta with a sequence of interest (e.g. from hg38),
#' and give me an assembly fasta, and I will return the fasta of a corresponding
#' sequence in the assembly. sas
#'
#' @param targetfasta [character/link] link to the 'target' single-sequence fasta (sometimes reference, e.g. chm13.)
#' @param queryfasta [character/link] link to the 'query' single-sequence fasta.
#' @return nothing. But output files written.
#'
#' @author Wolfram Höps
#' @rdname alignment
#' @export
# liftover_segment_direct <- function(targetfasta, queryfasta){
#
#   # Queryfasta: e.g. a piece of hg38 sequence.
#   # Targetfasta: typically a de-novo alignment.
#
#   queryfasta = outfasta
#   targetfasta = "/Users/hoeps/PhD/projects/huminvs/genomes/hifi-asm/HG00512_hgsvc_pbsq2-ccs_1000-hifiasm.h1-un.fasta"
#   outpaf_tmp = "./tpm_out.paf"
#   # Make an alignment with mmap2.
#   # We align the sequence we want to find to the assembly (targetfasta)
#   run_minimap2(queryfasta, targetfasta, outpaf_tmp)
#
#
# }

# outpaf_tmp = '/Users/hoeps/phd/projects/nahrcall/nahrchainer/seqbuilder/res/oa1.chunk'
# outpaf_tmp2 = '/Users/hoeps/phd/projects/nahrcall/nahrchainer/seqbuilder/res/oa1.chunk.filter'
#
# filter_paf_to_main_region(outpaf_tmp, outpaf_tmp2, inlen = 200000)

# To document



#' filter_paf_to_main_region
#' Script to filter down a whole-genome paf to the 'best matching' region.
#' query should be short sequence fragments (1kb, 10kb, ...) of the roi.

#' @author Wolfram Hoeps
#' @rdname paf_filter
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
    paf_winners[paf_winners$tname == chr_winners, ]$tend +
      paf_winners[paf_winners$tname == chr_winners, ]$tstart
  ) / 2))
  
  # Now simply extend from median towards front and back.
  start_winners = middle_median - (1.5 * (insequence_len / 2))
  end_winners = middle_median + (1.5 * (insequence_len / 2))
  
  # Here we cut, but I am not sure if this is really the way to go?
  paf_cut = paf[((paf$tname == chr_winners) &
                   (paf$tstart >= start_winners) &
                   (paf$tend <= end_winners)), ]
  
  
  aln_containing = cpaf_i[(cpaf_i$qstart <= start) &
                            (cpaf_i$qend >= end), ]
  
  
  
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

# hg38fa = "/Users/hoeps/PhD/projects/huminvs/genomes/hg38/hg38.fa"
# outfasta = "~/Desktop/trash.fa"
# seqname = 'chr11'
# start = 13104252
# end = 13122521
#
# extract_subseq_bedtools(hg38fa, seqname, start, end, outfasta)
#


# TRASH #

#' extract_subseq_bedtools

#' @author Wolfram Hoeps
#' @rdname alignment
#' @export
extract_subseq_bedtools <-
  function(infasta, seqname, start, end, outfasta) {
    # Needed for parallel runs
    
    random_tag = as.character(runif(1, 1e10, 1e11))
    tmp_bedfile = paste0('region2_', random_tag, '.bed')
    region = paste0(seqname, "\t", start, "\t", end)
    system(paste0('echo "', region, '" > ', tmp_bedfile))
    system(
      paste0(
        "/usr/local/bin/bedtools getfasta -fi ",
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

# determine_case <- function(cpaf_intersect, aln_containing){
#
#   # Case A: at least 1 unbroken alignment
#   if (dim(aln_containing)[1] >= 1){
#     return("A")
#   }
#
#
#   else if ()
# }

#' liftover_coarse


#' liftover_coarse <-
#'   function(seqname,
#'            start,
#'            end,
#'            conversionpaf_link,
#'            factor = 0.5) {
#'     start = 74769950
#'     end = 76058098
#'     
#'     start = 75000001
#'     end = 75010001
#'     cpaf = read.table(
#'       conversionpaf_link,
#'       sep = '\t',
#'       fill = T,
#'       row.names = NULL
#'     )
#'     
#'     colnames_paf = c(
#'       'qname',
#'       'qlen',
#'       'qstart',
#'       'qend',
#'       'strand',
#'       'tname',
#'       'tlen',
#'       'tstart',
#'       'tend',
#'       'nmatch',
#'       'alen',
#'       'mapq'
#'     )
#'     colnames(cpaf)[1:length(colnames_paf)] = colnames_paf
#'     
#'     # Query: from (..hg38..)
#'     # Target: to ( ..assembly.. )
#'     
#'     cpaf_i = cpaf[(cpaf$qname == seqname), ]
#'     
#'     
#'     # Get alignments that
#'     cpaf_intersect = cpaf_i[(cpaf_i$qend >= start) &
#'                               (cpaf_i$qstart <= end), ]
#'     
#'     aln_containing = cpaf_i[(cpaf_i$qstart <= start) &
#'                               (cpaf_i$qend >= end), ]
#'     
#'     # determine case
#'     case = determine_case(cpaf_intersect, aln_containing)
#'     
#'     
#'     # Start coord
#'     alns_containing_my_start = cpaf_i[(cpaf_i$qend >= start) &
#'                                         (cpaf_i$qstart <= start), ]
#'     longest_alns_containing_my_start = alns_containing_my_start[order(alns_containing_my_start$nmatch, decreasing =
#'                                                                         F), ][1, ]
#'     
#'     # End coord
#'     alns_containing_my_end = cpaf_i[(cpaf_i$qend >= end) &
#'                                       (cpaf_i$qstart <= end), ]
#'     longest_aln_containing_my_end = alns_containing_my_end[order(alns_containing_my_end$nmatch, decreasing =
#'                                                                    F), ][1, ]
#'     
#'     if (longest_alns_containing_my_start$qname != longest_aln_containing_my_end$qname) {
#'       if (longest_alns_containing_my_start$nmatch > longest_aln_containing_my_end$nmatch) {
#'         longest_aln_containing_my_end = longest_alns_containing_my_start
#'       } else {
#'         longest_aln_containing_my_start = longest_alns_containing_my_end
#'       }
#'     }
#'     
#'     
#'     if (longest_alns_containing_my_start$strand == '+') {
#'       liftover_start = longest_alns_containing_my_start$tstart + (start - longest_alns_containing_my_start$qstart)
#'     } else {
#'       liftover_start = longest_alns_containing_my_start$tend - (start - longest_alns_containing_my_start$qstart)
#'     }
#'     
#'     
#'     
#'     if (longest_aln_containing_my_end$strand == '+') {
#'       liftover_end = longest_aln_containing_my_end$tend - (longest_aln_containing_my_end$qend - end)
#'     } else {
#'       liftover_end = longest_aln_containing_my_end$tstart + (longest_aln_containing_my_end$qend - end)
#'     }
#'     
#'     
#'     liftover_start_real = as.integer(min(liftover_start, liftover_end) - ((end -
#'                                                                              start) * factor))
#'     liftover_end_real = as.integer(max(liftover_start, liftover_end) + ((end -
#'                                                                            start) * factor))
#'     
#'     
#'     to_return = list(
#'       lift_contig = alns_containing_my_start$tname,
#'       lift_start = max(liftover_start_real, 0),
#'       lift_end = min(liftover_end_real, alns_containing_my_start$tlen)
#'     )
#'     
#'     return(to_return)
#'   }

#' @export
find_punctual_liftover <- function(cpaf, pointcoordinate) {
  
  # Find all alignments overlapping with the pointcoordinate
  overlappers = cpaf[((cpaf$qname == seqname) &
                        (cpaf$qstart <= pointcoordinate) &
                        (cpaf$qend >= pointcoordinate)
  ), ]
  
  # If no alignment overlaps, we return an empty value. 
  if (dim(overlappers)[1] == 0){
    return(NULL)
  }
  
  # Find the alignment with the hightest nmatch value. 
  best_aln = overlappers[overlappers$nmatch == max(overlappers$nmatch), ]
  
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
                              find_punctual_liftover(cpaf, pointcoord))
    }
    
    # What is the majority vote for the target chromosome?
    winner_chr = names(sort(table(liftover_coords$seqname), decreasing = TRUE)[1])
    
    # And those snipplets matching to the winner chromosome - where do they fall (median)?
    middle_median = median(liftover_coords[liftover_coords$seqname == winner_chr, ]$liftover_coord)
    
    # Now  extend from median towards front and back.
    insequence_len = end - start
    start_winners = middle_median - (lenfactor * (insequence_len / 2))
    end_winners =   middle_median +   (lenfactor * (insequence_len / 2))
    
    # Make sure we don't exceed chromosome boundaries
    start_winners_cutoff = as.integer(max(0, start_winners))
    end_winners_cutoff =   as.integer(min(cpaf[cpaf$tname == winner_chr, ][1, 'tlen'], end_winners))
    
    if( (start_winners < 0) | (end_winners > cpaf[cpaf$tname == winner_chr, ][1, 'tlen'])){
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
# infasta = "/Users/hoeps/PhD/projects/huminvs/genomes/hg38/hg38.fa"
# outfasta = "~/Desktop/trash.fa"
# seqname = 'chr11'
# start = 13104252
# end = 13122521
#
# extract_subseq(infasta, seqname, start, end, outfasta)
