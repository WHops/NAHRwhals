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

#' flip_query_target
#'
#' Not entirely sure this is still in use. Might be obsolete.
#' @author Wolfram Höps
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



#' extract_subseq_bedtools
#' Give me an input fasta link, chromosome coordinates and an
#' output destionation file. I will extract the sequence of that
#' coordinates from the fasta, and output it into the outfasta.
#' @author Wolfram Hoeps
#' @export
extract_subseq_bedtools <-
  function(infasta, seqname, start, end, outfasta) {
    # Where is bedtools?
    bedtoolsloc = query_config("bedtools")
    
    # we need to create a temporary bedfile that will be deleted in a few lines.
    random_tag = as.character(runif(1, 1e10, 1e11))
    tmp_bedfile = paste0('region2_', random_tag, '.bed')
    
    # Write coordinates into temporary bedfile
    region = paste0(seqname,
                    "\t",
                    format(start, scientific = F),
                    "\t",
                    format(end, scientific = F))
    system(paste0('echo "', region, '" > ', tmp_bedfile))
    
    # Run bedtools
    system(paste0(
      bedtoolsloc,
      " getfasta -fi ",
      infasta,
      " -bed ",
      tmp_bedfile,
      " > ",
      outfasta
    ))
    
    # Delete tmp file
    system(paste0('rm ', tmp_bedfile))
    
    # Check if run was successful.
    stopifnot(
      "Error: Bedtools was unable to extract sequence. Please make sure a) Your input fasta paths are correct, b) target and queryfasta are different and c) your input fastas and the conversionpaf are matching." =
        file.size(outfasta) > 0
    )
    
    # Report your success :)
    print(paste0('Subsequence extracted and saved to ', outfasta))
    
  }


#' find_punctual_liftover
#' Translates a single point from one assembly to another.
#' Uses the pre-computed conversion-paf.
#' @export
find_punctual_liftover <- function(cpaf, pointcoordinate, chrom) {
  # Find all alignments overlapping with the pointcoordinate
  overlappers = cpaf[((cpaf$qname == chrom) &
                        (cpaf$qstart <= pointcoordinate) &
                        (cpaf$qend >= pointcoordinate)
  ),]
  
  # If no alignment overlaps, we return an empty value.
  if (dim(overlappers)[1] == 0) {
    return(c(NA, NA))
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

#' liftover_coarse
#'
#' Wrapperfunction to find a liftover sequence from one assembly to the other.
#' Better description TBD.
#'
#' @author Wolfram Höps
#' @export
liftover_coarse <-
  function(seqname,
           start,
           end,
           conversionpaf_link,
           n_probes = 100,
           lenfactor = 1.2, #Typically 1.0, because we want to get symmetrical dotplots.
           whole_chr = F) {
    

    cpaf = read_and_prep_paf(conversionpaf_link)
    # Assert that the input coordinates are within chromosome boundaries
    # stopifnot("Error: Input coordinates exceed input chromosome length." =
    #             ((end <= cpaf[cpaf == seqname,][1, 'qlen']) &
    #                (start >= 0)))
    
    
    # We take n_probes single points from the alignment, see where they fall, and take their median
    # as the median of our alignment.
    
    # Define a vector of probe positions (equally spaced between start and end)
    pos_probes = as.integer(start + (((end - start) / (n_probes - 1)) * (0:(n_probes -
                                                                              1))))
    
    # Liftover every probe
    liftover_coords <- data.frame(matrix(ncol = 3, nrow = length(pos_probes)))
    colnames(liftover_coords) <- c('pos_probe','seqname', 'liftover_coord')
    liftover_coords$pos_probe = pos_probes
    counter = 1
    #browser()
    for (pointcoord in pos_probes) {
      liftover_coords[counter, c('seqname','liftover_coord')] = 
                              find_punctual_liftover(cpaf, pointcoord, seqname)
      counter = counter + 1
    }
    

    #browser()
    liftover_coords = na.omit(liftover_coords)
    # Debug Tent
    library(ggplot2)
    ggplot(liftover_coords) + geom_point(aes(x=pos_probe, y=liftover_coord))
    
    # Make sure we got any results.
    stopifnot(
      "Error: Unsuccessful liftover of the input sequence: Sequence not found on query" =
        dim(liftover_coords)[1] > 0
    )
    
    
    # What is the majority vote for the target chromosome?
    winner_chr = names(sort(table(liftover_coords$seqname), decreasing = TRUE)[1])
    
    extrapolation_mode = T
    # If we want to get a whole chromosome (tested for chrY only), we start at 1, and
    # go as far as probes land (... thusly avoiding heterochromatin.)
    if (whole_chr){
      start_winners = 1
      end_winners = max(liftover_coords[liftover_coords$seqname == winner_chr,]$liftover_coord)
      
    } else if (extrapolation_mode ==F) {
      # And those snipplets matching to the winner chromosome - where do they fall (median)?
      middle_median = median(liftover_coords[liftover_coords$seqname == winner_chr,]$liftover_coord)
      
      # Now  extend from median towards front and back.
      insequence_len = end - start
      start_winners = middle_median - (lenfactor * (insequence_len / 2))
      end_winners =   middle_median +   (lenfactor * (insequence_len / 2))
      
      # Warn if we are exceeding chromosome boundaries in the query.
      if ((start_winners < 0) |
          (end_winners > cpaf[cpaf$tname == winner_chr,][1, 'tlen'])) {
        print('Warning! Reaching the end of alignment!')
    }
      # Make an entry to the output logfile #
      if (exists('log_collection')){
        log_collection$exceeds_y <<- T
      }
      # Log file entry done #
      
    } else if (extrapolation_mode) {
      
      # And those snipplets matching to the winner chromosome - where do they fall (median)?
      liftover_coords_maxseq = liftover_coords[liftover_coords$seqname == winner_chr,]
      
      # Probes are sorted. 
      mapping_direction = sign(sum(sign(diff(liftover_coords_maxseq$liftover_coord))))
      
      middle_median = median(liftover_coords_maxseq$liftover_coord)
      mapped_x_region = c(min(liftover_coords_maxseq$pos_probe), max(liftover_coords_maxseq$pos_probe))
      mapped_x_region_mid = median(mapped_x_region)
      
      extend_median_bwd = mapped_x_region_mid - start
      extend_median_fwd = end - mapped_x_region_mid
      
      # Make sure we got any results.
      stopifnot(
        "Error: No mapping direction recognized." =
          mapping_direction %in% c(1,-1)
      )
      
      if (mapping_direction == 1){
        start_winners = middle_median - extend_median_bwd
        end_winners = middle_median + extend_median_fwd
      } else if (mapping_direction == -1){
        start_winners = middle_median - extend_median_fwd
        end_winners = middle_median + extend_median_bwd
      } else {

      }
      
      # Warn if we are exceeding chromosome boundaries in the query.
      if ((start_winners < 0) |
          (end_winners > cpaf[cpaf$tname == winner_chr,][1, 'tlen'])) {
        print('Warning! Reaching the end of alignment!')
      }
      # Make an entry to the output logfile #
      if (exists('log_collection')){
        log_collection$exceeds_y <<- T
      }
      # Log file entry done #
      
    }
    
    # Make sure we don't exceed chromosome boundaries in the query.
    start_winners_cutoff = as.integer(max(0, start_winners))
    end_winners_cutoff =   as.integer(min(cpaf[cpaf$tname == winner_chr,][1, 'tlen'], end_winners))
    
    
    
    return(
      list(
        lift_contig = winner_chr,
        lift_start = start_winners_cutoff,
        lift_end = end_winners_cutoff
      )
    )
    
  }


#' enlarge_interval_by_factor
#'
#' @author Wolfram Höps
#' @export
enlarge_interval_by_factor <-
  function(start,
           end,
           factor,
           seqname_f = NULL,
           conversionpaf_f = NULL,
           log_collection = NULL) {
    
    # Some border cases #
    
    # end should be larger than start
    stopifnot("Error: interval end must be larger than start." = (end > start))
    
    # Factors smaller 0 make no sense
    stopifnot("Error: xpad factor must be larger than 0." = factor > 0)
    
    # Factors smaller 1 make little sense
    if (factor < 1) {
      warning(
        "xpad factors below 1 are typically not a good idea. (You are shorening your input sequence). xpad=1 means no padding."
      )
    }
    
    # If 1, don't do anything.
    if (factor == 1) {
      return(c(start, end))
    }
    
    # How long is interval
    interval_len = end - start
    desired_len = interval_len * factor
    
    # DO THE PADDING
    start_pad = round(start - (abs(interval_len - desired_len) / 2))
    end_pad =   round(end + (abs(interval_len - desired_len) / 2))
    
    
    # From here on: ensure that the padding was good :)
    if (start_pad < 0) {
      print('Warning: Start coordinate after padding exceeds chromosome boundary. Setting start to 0')
      start_pad = 0
      
      if (exists('log_collection')){
        log_collection$exceeds_x <<- T
      }
      
    }
    
    # IF conversionpaf is given (recommended), it is used to check how long
    # that chromosome actually is - and if our padding exceeds boundaries.
    if (!is.null(conversionpaf_f)) {
      stopifnot("Sequence padding: chromosome name is required for QC" = !is.na(seqname_f))
      
      paf = read_and_prep_paf(conversionpaf_f)
      
      chrlen = paf[paf$qname == seqname_f, ][1, 'qlen']
      if (end_pad > chrlen) {
        print(
          paste0(
            'Warning: End coordinate after padding (',
            end_pad,
            ') exceedis chromosome boundary. Setting end to chr end (',
            chrlen,
            ') instead.'
          )
        )
        end_pad = chrlen
        
        # Make an entry to the output logfile
        if (exists('log_collection')){
          log_collection$exceeds_x <<- T
        }
        
      }
    }
    
    # Make an entry to the output logfile #
    if (exists('log_collection')){
      log_collection[c('start_pad', 'end_pad')] <<- c(start_pad, end_pad)
    }
    # Log file entry done #
    
    return(c(start_pad, end_pad))
  }

#' read_and_prep_paf
#'
#' @author Wolfram Hoeps
#' @export
read_and_prep_paf <- function(paflink) {
  # Load conversionpaf
  cpaf = read.table(
    paflink,
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
  
  # Assert that the input coordinates are within chromosome boundaries
  return(cpaf)
}

