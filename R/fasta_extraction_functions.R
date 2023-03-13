#' extract_subseq_bedtools
#' Give me an input fasta link, chromosome coordinates and an
#' output destionation file. I will extract the sequence of that
#' coordinates from the fasta, and output it into the outfasta.
#' @author Wolfram Hoeps
#' @export
extract_subseq_bedtools <-
  function(infasta, seqname, start, end, outfasta) {

    print(infasta)
    print(seqname)
    print(outfasta)
    # Where is bedtools?
    bedtoolsloc = query_config("bedtools")
    print(bedtoolsloc)
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

    print(paste0(
      bedtoolsloc,
      " getfasta -fi ",
      infasta,
      " -bed ",
      tmp_bedfile,
      " > ",
      outfasta
    ))

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
  
  if(dim(best_aln)[1] > 1){
    best_aln = best_aln[1,]
  }
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
           whole_chr = F,
           search_mode = 'mad',
           refine_runnr = 0) {
    

    cpaf = read_and_prep_paf(conversionpaf_link)

    if (refine_runnr == 2){

      colnames_orig = colnames(cpaf)
      cpaf = cpaf[,c(6,7,8,9,5,1,2,3,4,10,11,12)]
      colnames(cpaf) = colnames_orig
      
      seqname = cpaf[1,'qname']
      start = 0
      end = cpaf[1,'qlen']
      
      # Somewhat risky here...
      aln_start = as.numeric(names(sort(table(sub(".*:([0-9]+).*", '\\1', cpaf$tname)),decreasing=TRUE)[1]))
      cpaf$tname = names(sort(table(sub(":[0-9]+-[0-9].*", "", cpaf$tname)),decreasing=TRUE)[1])
      cpaf$tstart = cpaf$tstart + aln_start
      cpaf$tend = cpaf$tend + aln_start
      
    } 
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
    for (pointcoord in pos_probes) {
      liftover_coords[counter, c('seqname','liftover_coord')] = 
                              find_punctual_liftover(cpaf, pointcoord, seqname)
      counter = counter + 1
    }
    
    liftover_coords = na.omit(liftover_coords)
    
    
    # Make sure we got any results.
    if (nrow(liftover_coords) <= 1){
      print("Error: Unsuccessful liftover of the input sequence: Sequence not found on query.")
      return(NULL)
    }
    
    # What is the majority vote for the target chromosome?
    winner_chr = names(sort(table(liftover_coords$seqname), decreasing = TRUE)[1])
    
    # If we want to get a whole chromosome (tested for chrY only), we start at 1, and
    # go as far as probes land (... thusly avoiding heterochromatin.)
    # whole_chr = F

    # if (whole_chr){
    #   start_winners = 1
    #   end_winners = min(28500000,max(liftover_coords[liftover_coords$seqname == winner_chr,]$liftover_coord))
      
    
    
  if (search_mode == 'extrapolation') {
      startend = find_coords_extrapolated(liftover_coords, cpaf, winner_chr, start, end)
      
    } else if (search_mode == 'mad'){ # mad = mean absolute deviation.
      startend = find_coords_mad(liftover_coords, cpaf, winner_chr, start, end, refine_runnr!=1)
      if (is.null(startend)){
        print('MAD failed. Going for extrapolation instead.')
        startend = find_coords_extrapolated(liftover_coords, cpaf, winner_chr, start, end, refine_runnr!=1)
        
      }
      
    }
    
    if (refine_runnr == 1){
      length = startend[2] - startend[1]
      if (length < 50000){
        shift_fact = 1.5
      } else if (length < 200000){
        shift_fact = 1.2
      } else {
        shift_fact = 1.1
      }
      start_winners = startend[1] - (length*((shift_fact-1) * 0.5))
      end_winners = startend[2] + (length*((shift_fact-1) * 0.5))
      # Make sure we don't exceed chromosome boundaries in the query.
      start_winners_cutoff = as.integer(max(0, start_winners))
      end_winners_cutoff =   as.integer(min(cpaf[cpaf$tname == winner_chr,][1, 'tlen'] - 1, end_winners - 1))
    } else if (refine_runnr == 2) {
      start_winners_cutoff = as.integer(max(0, startend[1]))
      end_winners_cutoff = as.integer(min(max(cpaf[cpaf$tname == winner_chr,'tend']), startend[2] - 1))
    } else {
      start_winners_cutoff = as.integer(max(0, startend[1]))
      end_winners_cutoff = as.integer(min(cpaf[cpaf$tname == winner_chr,][1, 'qlen'] + startend[1] - 1, startend[2] - 1))
    }
    
    
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
#' input start and end coordinates and a factor by which to enlarge. E.g. 40-50 (len 10) and factor 2 will result in 35-55 (len 20).
#' Makes sure no negative values are created. 
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
    if ((!is.null(conversionpaf_f)) & (!is.na(conversionpaf_f))) {
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


#' find_coords_mad
#'
#' @author Wolfram Hoeps
#' @export
find_coords_mad <- function(liftover_coords, cpaf, winner_chr, start, end, liftover_run){

  # Hardcode success thresholds
  th_x = 0.8
  th_y_low = 0.5
  th_y_high = 2
  start_map_limit = 10
  end_map_limit = 90

  # And those snipplets matching to the winner chromosome - where do they fall (median)?
  liftover_coords_maxseq = liftover_coords[liftover_coords$seqname == winner_chr,]
  
  # Filter outliers
  liftover_coords_maxseq = liftover_coords_maxseq[mad_mask_outliers(liftover_coords_maxseq$liftover_coord),]
  
  # The region is now the region
  
  start_winners = min(liftover_coords_maxseq$liftover_coord)
  end_winners = max(liftover_coords_maxseq$liftover_coord)
  
  # library(ggplot2)
  # p = ggplot(liftover_coords) + geom_point(aes(x=pos_probe, y=liftover_coord)) +
  # coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
  #print(p)
  mapped_x_region_frac = abs(max(liftover_coords_maxseq$pos_probe) - min(liftover_coords_maxseq$pos_probe)) / abs(end-start)
  
  mapped_y_region_frac = abs(end_winners - start_winners) / abs(end-start)
  
  start_pointers_arrived = min(as.numeric(row.names(liftover_coords_maxseq))) < start_map_limit
  end_pointers_arrived = max(as.numeric(row.names(liftover_coords_maxseq))) > end_map_limit
  success = (mapped_x_region_frac > th_x) & (mapped_y_region_frac > th_y_low) & (mapped_y_region_frac < th_y_high) &
    (start_pointers_arrived) & (end_pointers_arrived)
  
  if ((success == F) & (liftover_run == F)){
    return(NULL)
  }
  
  if (liftover_run == F){
    
    dist_between_probes = min(unique(diff(liftover_coords$pos_probe)))
    n_probes_distance = 5
    # If the break is very close (5_probes distance)
    # Warn if we are exceeding chromosome boundaries in the query.
    if ((start_winners < (dist_between_probes*n_probes_distance)) |
        (end_winners + (dist_between_probes*n_probes_distance)) > (cpaf[cpaf$tname == winner_chr,][1, 'tlen'])) {
      print('Warning! Reaching the end of alignment!')
      
      # Make an entry to the output logfile #
      if (exists('log_collection')){
        log_collection$exceeds_y <<- T
      }
    }
  }


  return(c(start_winners, end_winners))

  
}


#' find_coords_extrapolated
#'
#' @author Wolfram Hoeps
#' @export
find_coords_extrapolated <- function(liftover_coords, cpaf, winner_chr, start, end, liftover_run){
  
  # And those snipplets matching to the winner chromosome - where do they fall (median)?
  liftover_coords_maxseq = liftover_coords[liftover_coords$seqname == winner_chr,]
  
  # If <10% of probes are matching, there is no point in continuing this discussion
  stopifnot("No sequence homolog found in assembly. Consider trying a larger locus window to gap potential deletions." = 
              dim(liftover_coords_maxseq)[1] > 10)

  # Probes are sorted. The juicer makes direction positive if in doubt. (The juicer is what you want it to be :P)
  the_juicer = 1e-5
  
  start_exists = all(1:6 %in% as.numeric(row.names(liftover_coords_maxseq)))
  end_exists = all(94:100 %in% as.numeric(row.names(liftover_coords_maxseq)))
  
  if (start_exists){
    mapping_direction = sign(sum(sign(diff(liftover_coords_maxseq[as.character(seq(1,6)),]$liftover_coord))) + the_juicer)
  } else if (end_exists){
    mapping_direction = sign(sum(sign(diff(liftover_coords_maxseq[as.character(seq(94,100)),]$liftover_coord))) + the_juicer)
  } else {
    mapping_direction = sign(sum(sign(diff(liftover_coords_maxseq$liftover_coord))) + the_juicer)
  }
  middle_median = median(liftover_coords_maxseq$liftover_coord)
  
  #ggplot(liftover_coords) + geom_point(aes(x=pos_probe, y=liftover_coord))# +
  #coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
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
  
  pointspace = min(diff(liftover_coords_maxseq$pos_probe))
  points_from_border_to_consider_ok = 5
  if (liftover_run == F){   
    # Warn if we are exceeding chromosome boundaries in the query.
    if ((start_winners < (pointspace * points_from_border_to_consider_ok)) |
       ( (end_winners + (pointspace * points_from_border_to_consider_ok))  > cpaf[cpaf$tname == winner_chr,][1, 'tlen'])) {
      print('Warning! Reaching the end of alignment!')
  
      # Make an entry to the output logfile #
      if (exists('log_collection')){
        log_collection$exceeds_y <<- T
      }
    }
  }
  # Log file entry done #
  
  return(c(start_winners, end_winners))
}

#' mad_mask_outliers
#'
#' @author Wolfram Hoeps
#' @export
mad_mask_outliers <- function(obs, th=4){

  obs_median = median(obs)
  
  dist_from_median = ((obs - obs_median)**2)
  dist_from_median = sqrt(dist_from_median)
  med_abs_deviation = median(dist_from_median)
  
  modified_z_score = 0.6745 * dist_from_median / med_abs_deviation
  
return(modified_z_score < th)

}
