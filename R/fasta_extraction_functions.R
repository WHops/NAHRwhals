#' extract_subseq_bedtools
#'
#' Extracts a subsequence from a fasta file using bedtools getfasta.
#'
#' @param infasta (character): A link to the input fasta file.
#' @param seqname (character): The name of the sequence/chromosome to extract from the fasta file.
#' @param start (integer): The starting position (1-based) of the sequence to extract.
#' @param end (integer): The ending position (1-based) of the sequence to extract.
#' @param outfasta (character): The destination file to write the extracted subsequence in fasta format.
#' @param params (list): A list containing the path to the bedtools binary (bedtools_bin).
#'
#' @author Wolfram Hoeps
#' @export
extract_subseq_bedtools <- function(infasta, seqname, start, end, outfasta, params) {
  # Where is bedtools?
  bedtoolsloc <- params$bedtools_bin
  
  # Generate a random tag for temporary files to avoid conflicts
  random_tag <- as.character(runif(1, 1e10, 1e11))
  
  # Temporary files for BED and sequence length information
  tmp_bedfile <- paste0("region2_", random_tag, ".bed")
  tmp_seq_length_file <- paste0("seq_length_", random_tag, ".txt")
  
  # Get sequence lengths using bedtools faidx and store in a temporary file
  system(paste0("samtools faidx ", infasta))
  system(paste0("cut -f1,2 ", infasta, ".fai > ", tmp_seq_length_file))
  
  # Read the sequence lengths into R
  seq_lengths <- read.table(tmp_seq_length_file, stringsAsFactors = FALSE, col.names = c("seqname", "length"))
  
  # Adjust the 'end' position if it exceeds the contig length
  contig_length <- seq_lengths[seq_lengths$seqname == seqname, "length"]
  if (!is.na(contig_length) && contig_length < end) {
    end <- contig_length
  }
  
  # Write coordinates into the temporary BED file
  write.table(data.frame(seqname, start, end), file = tmp_bedfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  # Run bedtools getfasta
  system(paste0(bedtoolsloc, " getfasta -fi ", infasta, " -bed ", tmp_bedfile, " -fo ", outfasta))
  
  # Delete temporary files
  file.remove(tmp_bedfile, tmp_seq_length_file)
  
  # Check if run was successful by confirming the output file size is greater than 0
  if (file.size(outfasta) <= 0) {
    stop("Error: Bedtools was unable to extract sequence. Please make sure a) Your input fasta paths are correct, b) target and queryfasta are different, and c) your input fastas and the conversionpaf are matching.")
  }
}

#' find_punctual_liftover
#'
#' Translates a single point from one genome assembly to another using a pre-computed conversion-paf.
#'
#' @param cpaf (data.frame): A data frame containing the conversion-paf.
#' @param pointcoordinate (integer): The coordinate of the point to be lifted over.
#' @param chrom (character): The chromosome of the point to be lifted over.
#'
#' @return A data frame with two columns: seqname and liftover_coord.
#' seqname gives the name of the sequence to which the point was lifted over, and liftover_coord
#' gives the lifted-over coordinate.
#'
#' @author Wolfram Hoeps
#' @export
find_punctual_liftover <- function(cpaf, pointcoordinate, chrom) {
  # Find all alignments overlapping with the pointcoordinate
  overlappers <- cpaf[((cpaf$qname == chrom) &
    (cpaf$qstart <= pointcoordinate) &
    (cpaf$qend >= pointcoordinate)
  ), ]

  # If no alignment overlaps, we return an empty value.
  if (dim(overlappers)[1] == 0) {
    return(c(NA, NA))
  }

  # Find the alignment with the hightest nmatch value.
  best_aln <- overlappers[overlappers$nmatch == max(overlappers$nmatch), ]

  if (dim(best_aln)[1] > 1) {
    best_aln <- best_aln[1, ]
  }
  # Liftover the point with that highest alignment.
  if (best_aln$strand == "+") {
    coord <- best_aln$tstart + (pointcoordinate - best_aln$qstart)
  } else {
    coord <- best_aln$tend - (pointcoordinate - best_aln$qstart)
  }

  return(data.frame(seqname = best_aln$tname, liftover_coord = coord))
}


#' liftover_coarse
#'
#' This function lifts over a genomic region from one assembly to another.
#' The function extracts subregions of the input sequence, and then
#' calculates the lift-over positions by aligning these subregions
#' to the target assembly. By default, the function uses a mean absolute
#' deviation (MAD) based approach to calculate the lifted-over region.
#' However, it can also be set to use an extrapolation-based approach
#' to calculate the lifted-over region. The result of the function is a
#' list containing the lifted-over contig name, the start position, and the
#' end position.
#'
#' @param seqname (character): The name of the contig/chromosome to lift-over.
#' @param start (integer): The starting position (1-based) of the region to lift-over.
#' @param end (integer): The ending position (1-based) of the region to lift-over.
#' @param conversionpaf_link (character): A link to the precomputed PAF file containing
#' the conversion information between the two genome assemblies.
#' @param n_probes (integer): The number of subregions to use for calculating the lifted-over
#' position. Defaults to 100.
#' @param lenfactor (numeric): The length factor used to calculate the length of the subregions.
#' A higher lenfactor will result in longer subregions, which may improve the quality of the lifted-over
#' position. Defaults to 1.2.
#' @param whole_chr (logical): Whether to lift-over the entire chromosome/contig. By default, this is set to
#' FALSE, and only the region specified by the start and end coordinates is lifted-over.
#' @param search_mode (character): The search mode to use for finding the lifted-over position.
#' Can be set to "mad" (default), or "extrapolation".
#' @param refine_runnr (integer): The number of times to refine the lifted-over position. Can be set to
#' 0 (default), 1, or 2. A higher number will result in a more accurate lifted-over position, but will
#' also increase the runtime.
#'
#' @return A list containing the lifted-over contig name, the start position, and the end position.
#'
#' @author Wolfram Hoeps
#' @export
liftover_coarse <-
  function(seqname,
           start,
           end,
           conversionpaf_link,
           external_paf_bool,
           n_probes = 100,
           lenfactor = 1.2, # Typically 1.0, because we want to get symmetrical dotplots.
           whole_chr = F,
           search_mode = "mad",
           refine_runnr = 0
           ) {
    cpaf <- read_and_prep_paf(conversionpaf_link)
    
    if (external_paf_bool){
      newnames = c("tname"  ,"tlen"  , "tstart" ,"tend" ,  
                   "strand" ,"qname",  "qlen" ,  "qstart" ,
                   "qend" ,  "nmatch" ,"alen" ,  "mapq")
      colnames(cpaf) = newnames
      
      cpaf$tname = sub("_[0-9]+-[0-9]+$", "", cpaf$tname)
      
    }
    
    if (refine_runnr == 2) {
      colnames_orig <- colnames(cpaf)
      cpaf <- cpaf[, c(6, 7, 8, 9, 5, 1, 2, 3, 4, 10, 11, 12)]
      colnames(cpaf) <- colnames_orig
      text = cpaf[1,'tname']
      pattern <- ":([0-9]+-[0-9]+)_"
      matches <- regexpr(pattern, text)
      result <- regmatches(text, matches)
      
      if (attr(matches, "match.length") > 0) {
        extracted <- substring(result, 2, nchar(result)-1)
      } else {
        extracted <- character(0)  # No match found
      }
      
      
      start_here = as.numeric(strsplit(extracted, '-')[[1]][1])
      end_here = as.numeric(strsplit(extracted, '-')[[1]][2])
      len_here = end_here - start_here
      len_orig = len_here * (1/1.04)
      start_orig = len_orig * 0.02
      end_orig = start_orig + len_orig
      
      seqname <- cpaf[1, "qname"]
      start <- 0
      end <- cpaf[1, "qlen"]
      
      # Somewhat risky here...
      aln_start <- as.numeric(names(sort(table(sub(".*:([0-9]+).*", "\\1", cpaf$tname)), decreasing = TRUE)[1]))
      cpaf$tname <- names(sort(table(sub(":[0-9]+-[0-9].*", "", cpaf$tname)), decreasing = TRUE)[1])
      cpaf$tstart <- cpaf$tstart + aln_start
      cpaf$tend <- cpaf$tend + aln_start
    }
    # We take n_probes single points from the alignment, see where they fall, and take their median
    # as the median of our alignment.
    # Define a vector of probe positions (equally spaced between start and end)
    pos_probes <- as.integer(start + (((end - start) / (n_probes - 1)) * (0:(n_probes -
      1))))
    # Liftover every probe

    liftover_coords <- data.frame(matrix(ncol = 3, nrow = length(pos_probes)))
    colnames(liftover_coords) <- c("pos_probe", "seqname", "liftover_coord")
    liftover_coords$pos_probe <- pos_probes
    counter <- 1
    #browser()
    for (pointcoord in pos_probes) {
      liftover_coords[counter, c("seqname", "liftover_coord")] <-
        find_punctual_liftover(cpaf, pointcoord, seqname)
      counter <- counter + 1
    }

    liftover_coords <- na.omit(liftover_coords)

    # Make sure we got any results.
    if (nrow(liftover_coords) <= 1) {
      print("Error: Unsuccessful liftover of the input sequence: Sequence not found on query.")
      return(NULL)
    }

    # What is the majority vote for the target chromosome?
    winner_chr <- names(sort(table(liftover_coords$seqname), decreasing = TRUE)[1])

    # If we want to get a whole chromosome (tested for chrY only), we start at 1, and
    # go as far as probes land (... thusly avoiding heterochromatin.)
    # whole_chr = F

    # if (whole_chr){
    #   start_winners = 1
    #   end_winners = min(28500000,max(liftover_coords[liftover_coords$seqname == winner_chr,]$liftover_coord))
    if (search_mode == "extrapolation") {
      startend <- find_coords_extrapolated(liftover_coords, cpaf, winner_chr, start, end)
    } else if (search_mode == "mad") { # mad = mean absolute deviation.
      startend <- find_coords_mad(liftover_coords, cpaf, winner_chr, start, end, refine_runnr != 1)
      if (is.null(startend)) {
        startend <- find_coords_extrapolated(liftover_coords, cpaf, winner_chr, start, end, refine_runnr != 1)
      }
    }

    if (is.null(startend)){
      return(NULL)
    }
    if (external_paf_bool){
      start_winners_cutoff <- as.integer(max(0, startend[1]))
      end_winners_cutoff <- as.integer(min(max(cpaf[cpaf$tname == winner_chr, "tend"]), startend[2] - 1))
      
      return(
        list(
          lift_contig = winner_chr,
          lift_start = start_winners_cutoff,
          lift_end = end_winners_cutoff
        )
      )
    }
    
    if (refine_runnr == 1) {
      length <- startend[2] - startend[1]
      if (length < 50000) {
        shift_fact <- 1.5
      } else if (length < 200000) {
        shift_fact <- 1.2
      } else {
        shift_fact <- 1.1
      }
      
      start_winners <- startend[1] - (length * ((shift_fact - 1) * 0.5))
      end_winners <- startend[2] + (length * ((shift_fact - 1) * 0.5))

      # Make sure we don't exceed chromosome boundaries in the query.
      start_winners_cutoff <- as.integer(max(0, start_winners))
      end_winners_cutoff <- as.integer(min(cpaf[cpaf$tname == winner_chr, ][1, "tlen"] - 1, end_winners - 1))


    } else if (refine_runnr == 2) {
      start_winners_cutoff <- as.integer(max(0, startend[1]))
      end_winners_cutoff <- as.integer(min(max(cpaf[cpaf$tname == winner_chr, "tend"]), startend[2] - 1))

      

    } else {
      start_winners_cutoff <- as.integer(max(0, startend[1]))
      end_winners_cutoff <- as.integer(min(cpaf[cpaf$tname == winner_chr, ][1, "qlen"] + startend[1] - 1, startend[2] - 1))
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
#' Given a start coordinate, an end coordinate and a factor, the function returns new start and end coordinates after padding the original interval. The padding factor can be greater than 1, in which case the function increases the length of the input interval by that factor. The padding is applied symmetrically on both ends. The function makes sure that negative values are not created and that the output interval does not exceed the length of the chromosome. A log file can also be created to keep track of the operation.
#'
#' @param start start coordinate of the input interval
#' @param end end coordinate of the input interval
#' @param factor factor by which to increase the length of the input interval
#' @param seqname_f (optional) chromosome name, used for quality control
#' @param conversionpaf_f (optional) path to PAF file, used for quality control
#' @param log_collection (optional) log file to keep track of the operation
#' @return A vector of length 2, containing the new start and end coordinates.
#' @author Wolfram Hoeps
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
    interval_len <- end - start
    desired_len <- interval_len * factor

    # DO THE PADDING
    start_pad <- round(start - (abs(interval_len - desired_len) / 2))
    end_pad <- round(end + (abs(interval_len - desired_len) / 2))


    # From here on: ensure that the padding was good :)
    if (start_pad < 0) {
      #print("Warning: Start coordinate after padding exceeds chromosome boundary. Setting start to 0")
      start_pad <- 0

      if (exists("log_collection")) {
        log_collection$exceeds_x <<- T
      }
    }

    # IF conversionpaf is given (recommended), it is used to check how long
    # that chromosome actually is - and if our padding exceeds boundaries.
    if ((!is.null(conversionpaf_f)) & (!is.na(conversionpaf_f))) {
      stopifnot("Sequence padding: chromosome name is required for QC" = !is.na(seqname_f))

      paf <- read_and_prep_paf(conversionpaf_f)

      chrlen <- paf[paf$qname == seqname_f, ][1, "qlen"]
      if (end_pad > chrlen) {
        #print(
        #  paste0(
        #    "Warning: End coordinate after padding (",
        #    end_pad,
        #    ") exceedis chromosome boundary. Setting end to chr end (",
        #    chrlen,
        #    ") instead."
        #  )
        #)
        end_pad <- chrlen

        # Make an entry to the output logfile
        if (exists("log_collection")) {
          log_collection$exceeds_x <<- T
        }
      }
    }

    # Make an entry to the output logfile #
    if (exists("log_collection")) {
      log_collection[c("start_pad", "end_pad")] <<- c(start_pad, end_pad)
    }
    # Log file entry done #

    return(c(start_pad, end_pad))
  }

#' read_and_prep_paf
#'
#' Reads in a PAF file and renames the columns to match the PAF standard column names.
#' @param paflink The path to the PAF file to read.
#' @author Wolfram Hoeps
#' @export
read_and_prep_paf <- function(paflink) {
  # Load conversionpaf
  cpaf <- read.table(
    paflink,
    sep = "\t",
    fill = T,
    row.names = NULL
  )

  colnames_paf <- c(
    "qname",
    "qlen",
    "qstart",
    "qend",
    "strand",
    "tname",
    "tlen",
    "tstart",
    "tend",
    "nmatch",
    "alen",
    "mapq"
  )
  colnames(cpaf)[1:length(colnames_paf)] <- colnames_paf

  # Assert that the input coordinates are within chromosome boundaries
  return(cpaf)
}

#' find_coords_sane

#' find_coords_mad
#' This function finds the coordinates of a region given a set of liftover coordinates and other parameters.
#' @param liftover_coords A data frame containing liftover coordinates
#' @param cpaf A data frame containing cpaf information
#' @param winner_chr The chromosome with the highest number of matching sequences
#' @param start The start position of the region
#' @param end The end position of the region
#' @param liftover_run A logical value indicating whether this is a liftover run
#' @return A numeric vector containing the start and end coordinates of the region
#' @author Wolfram Hoeps
#' @export
find_coords_mad <- function(liftover_coords, cpaf, winner_chr, start, end, liftover_run) {
  # Hardcode success thresholds
  th_x <- 0.8
  th_y_low <- 1/3
  th_y_high <- 3
  start_map_limit <- 10
  end_map_limit <- 90

  # And those snipplets matching to the winner chromosome - where do they fall?
  liftover_coords_maxseq <- liftover_coords[liftover_coords$seqname == winner_chr, ]

  # Filter outliers
  # Here is that part...
  liftover_coords_maxseq <- liftover_coords_maxseq[mad_mask_outliers(liftover_coords_maxseq$liftover_coord), ]

  direction = sign(sum(sign(diff(liftover_coords_maxseq$liftover_coord))) + 1e-5)

  steplength = abs(min(abs(diff(liftover_coords_maxseq$liftover_coord - min(abs(liftover_coords_maxseq$liftover_coord))))))


  start_winners <- liftover_coords_maxseq[1,'liftover_coord']# min(liftover_coords_maxseq$liftover_coord)
  end_winners <- liftover_coords_maxseq[nrow(liftover_coords_maxseq),'liftover_coord']#max(liftover_coords_maxseq$liftover_coord)

  # p = ggplot(liftover_coords) + geom_point(aes(x=pos_probe, y=liftover_coord)) +
  # coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
  mapped_x_region_frac <- abs(max(liftover_coords_maxseq$pos_probe) - min(liftover_coords_maxseq$pos_probe)) / abs(end - start)

  mapped_y_region_frac <- abs(end_winners - start_winners) / abs(end - start)

  start_pointers_arrived <- min(as.numeric(row.names(liftover_coords_maxseq))) < start_map_limit
  end_pointers_arrived <- max(as.numeric(row.names(liftover_coords_maxseq))) > end_map_limit
  
  start_is_in_beginning <- which.min(liftover_coords$liftover_coord) < start_map_limit
  end_is_in_end <- which.max(liftover_coords$liftover_coord) > end_map_limit
  
  start_is_in_end <- which.min(liftover_coords$liftover_coord) > end_map_limit
  end_is_in_start <- which.max(liftover_coords$liftover_coord) < start_map_limit
  
  start_end_ok = (start_is_in_beginning && end_is_in_end) | (start_is_in_end && end_is_in_start)
  
  success <- (mapped_x_region_frac > th_x) & (mapped_y_region_frac > th_y_low) & (mapped_y_region_frac < th_y_high) &
    (start_pointers_arrived) & (end_pointers_arrived) #& start_end_ok

  if ((success == F) & (liftover_run == F)) {
    print("exiting")
    return(NULL)
  }


  mode = 'letstry'
  if (mode == 'pre_december'){


  if (direction == 1){
    start_winners = liftover_coords_maxseq[1,'liftover_coord'] - 
      ((min(as.numeric(rownames(liftover_coords_maxseq)))-1) * steplength)

    end_winners = liftover_coords_maxseq[nrow(liftover_coords_maxseq),'liftover_coord'] + 
      ((100 - max(as.numeric(rownames(liftover_coords_maxseq)))) *  steplength)
      
  } else if (direction == -1){

    start_winners = liftover_coords_maxseq[nrow(liftover_coords_maxseq),'liftover_coord'] - 
    ((min(as.numeric(rownames(liftover_coords_maxseq)))-1) * steplength)

    end_winners = liftover_coords_maxseq[1,'liftover_coord'] + 
    ((100 - max(as.numeric(rownames(liftover_coords_maxseq)))) *  steplength)

  }
  
  } else {


      start_winners = min(liftover_coords_maxseq[,'liftover_coord'])# - 
       # ((min(as.numeric(rownames(liftover_coords_maxseq)))-1) * steplength)

      end_winners = max(liftover_coords_maxseq[,'liftover_coord'])# + 
        #((100 - max(as.numeric(rownames(liftover_coords_maxseq)))) *  steplength)
        
  } 

  
  
  if (liftover_run == F) {
    dist_between_probes <- min(unique(diff(liftover_coords$pos_probe)))
    n_probes_distance <- 1

    # If less than 90% of probes fall on one target, we get suspicious
    if (max(as.numeric(table(liftover_coords$seqname))) < 90){
      # If the break is very close (5_probes distance)
      # Warn if we are exceeding chromosome boundaries in the query.
      if ((start_winners < (dist_between_probes * n_probes_distance)) |
        (end_winners + (dist_between_probes * n_probes_distance)) > (cpaf[cpaf$tname == winner_chr, ][1, "tlen"])) {
        #print("Warning! Reaching the end of alignment!")
  
        # Make an entry to the output logfile #
        if (exists("log_collection")) {
          log_collection$exceeds_y <<- T
        }
      }
    }
  }



  return(c(start_winners, end_winners))
}


#' find_coords_extrapolated
#' This function finds the coordinates of a region given a set of liftover coordinates and other parameters.
#' @param liftover_coords A data frame containing liftover coordinates
#' @param cpaf A data frame containing cpaf information
#' @param winner_chr The chromosome with the highest number of matching sequences
#' @param start The start position of the region
#' @param end The end position of the region
#' @param liftover_run A logical value indicating whether this is a liftover run
#' @return A numeric vector containing the start and end coordinates of the region
#' @author Wolfram Hoeps
#' @export
find_coords_extrapolated <- function(liftover_coords, cpaf, winner_chr, start, end, liftover_run) {
  # And those snipplets matching to the winner chromosome - where do they fall (median)?
  liftover_coords_maxseq <- liftover_coords[liftover_coords$seqname == winner_chr, ]

  # If <10% of probes are matching, there is no point in continuing this discussion
  if (dim(liftover_coords_maxseq)[1] < 10){
    return(NULL)
  }

  # Probes are sorted. The juicer makes direction positive if in doubt. (The juicer is what you want it to be :P)
  the_juicer <- 1e-5

  start_exists <- all(1:6 %in% as.numeric(row.names(liftover_coords_maxseq)))
  end_exists <- all(94:100 %in% as.numeric(row.names(liftover_coords_maxseq)))

  if (start_exists) {
    mapping_direction <- sign(sum(sign(diff(liftover_coords_maxseq[as.character(seq(1, 6)), ]$liftover_coord))) + the_juicer)
  } else if (end_exists) {
    mapping_direction <- sign(sum(sign(diff(liftover_coords_maxseq[as.character(seq(94, 100)), ]$liftover_coord))) + the_juicer)
  } else {
    mapping_direction <- sign(sum(sign(diff(liftover_coords_maxseq$liftover_coord))) + the_juicer)
  }
  middle_median <- median(liftover_coords_maxseq$liftover_coord)

  # ggplot(liftover_coords) + geom_point(aes(x=pos_probe, y=liftover_coord))# +
  # coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
  mapped_x_region <- c(min(liftover_coords_maxseq$pos_probe), max(liftover_coords_maxseq$pos_probe))
  mapped_x_region_mid <- median(mapped_x_region)

  extend_median_bwd <- mapped_x_region_mid - start
  extend_median_fwd <- end - mapped_x_region_mid
  # Make sure we got any results.
  stopifnot(
    "Error: No mapping direction recognized." =
      mapping_direction %in% c(1, -1)
  )

  if (mapping_direction == 1) {
    start_winners <- middle_median - extend_median_bwd
    end_winners <- middle_median + extend_median_fwd
  } else if (mapping_direction == -1) {
    start_winners <- middle_median - extend_median_fwd
    end_winners <- middle_median + extend_median_bwd
  } else {

  }

  dist_between_probes <- min(unique(diff(liftover_coords$pos_probe)))
  n_probes_distance <- 1
  if (max(as.numeric(table(liftover_coords$seqname))) < 90){
    # If the break is very close (5_probes distance)
    # Warn if we are exceeding chromosome boundaries in the query.
    if ((start_winners < (dist_between_probes * n_probes_distance)) |
        (end_winners + (dist_between_probes * n_probes_distance)) > (cpaf[cpaf$tname == winner_chr, ][1, "tlen"])) {
      #print("Warning! Reaching the end of alignment!")
      
      # Make an entry to the output logfile #
      if (exists("log_collection")) {
        log_collection$exceeds_y <<- T
      }
    }
  }
  
  # Log file entry done #

  return(c(start_winners, end_winners))
}

#' mad_mask_outliers
#' This function masks outliers in a vector using the median absolute deviation (MAD) method.
#' @param obs A numeric vector
#' @param th A numeric threshold (default 4)
#' @return A logical vector indicating which values in obs are not outliers
#' @author Wolfram Hoeps
#' @export
mad_mask_outliers <- function(obs, th = 4) {
  obs_median <- median(obs)

  dist_from_median <- ((obs - obs_median)**2)
  dist_from_median <- sqrt(dist_from_median)
  med_abs_deviation <- median(dist_from_median)

  modified_z_score <- 0.6745 * dist_from_median / med_abs_deviation

  return(modified_z_score < th)
}
