#!/usr/local/bin/Rscript

#' Merge paired alignment rows in PAF file
#'
#' @description Merges paired alignment rows in the provided PAF file.
#'
#' @param paffile PAF file data frame with column names
#' @param pairs Data frame with columns representing row numbers of paired alignments
#' @return PAF file data frame with merged rows
#'
#' @author Wolfram Hoeps
#' @export
merge_rows <- function(paffile, pairs) {

  count <- 0

  # Accessing the large paffile should be minimized!
  # We pull out the relevant columns here instead - much faster to work with them!
  qnames = paffile$qname
  qstarts = paffile$qstart
  qends = paffile$qend
  tstarts = paffile$tstart
  tends = paffile$tend
  strands = paffile$strand
  nmatches = as.numeric(paffile$nmatch)
  alens = as.numeric(paffile$alen)
  
  for (i in dim(pairs)[1]:1) {
    #if (count %% 100 == 0) {
    #  print(paste0("Merging pair ", count, " out of ", dim(pairs)[1]))
    #}

    nl1 <- (pairs[i, 1])
    nl2 <- (pairs[i, 2])
    nl1_bu = nl1
    
    if (strands[nl1] == "+") {
      
      qname1 = qnames[nl1]
      qstart1 = qstarts[nl1]
      qstart2 = qstarts[nl2]
      qend1 = qends[nl1]
      qend2 = qends[nl2]
      tstart1 = tstarts[nl1]
      tstart2 = tstarts[nl2]
      tend1 = tends[nl1]
      tend2 = tends[nl2]
      
      qnames[nl1_bu] <- paste0(
        sub("_[^_]*$", "", qname1),
        "_",
        qstart1,
        "-",
        qend2 - 1
      )

      # query coordinates
      qends[nl1_bu] <- max(qend1, qend2)
      qstarts[nl1_bu] <- min(qstart1, qstart2)

      # target coordinates
      tends[nl1_bu] <- max(tend1, tend2)
      tstarts[nl1_bu] <- min(tstart1, tstart2)
      
      
    } else if (strands[nl1] == "-") {

      # Sort by t, because this is increasing
      if (tstarts[nl1] > tstarts[nl2]){
        # Switch the values of nl1 and nl2
        nl1 = nl2
        nl2 = nl1_bu
      }

      qname1 = qnames[nl1]
      qstart1 = qstarts[nl1]
      qstart2 = qstarts[nl2]
      qend1 = qends[nl1]
      qend2 = qends[nl2]
      tstart1 = tstarts[nl1]
      tstart2 = tstarts[nl2]
      tend1 = tends[nl1]
      tend2 = tends[nl2]

      
      qnames[nl1_bu] <- paste0(
        sub("_[^_]*$", "", qname1),
        "_",
        qend2,
        "-",
        qstart2
      )

      # query coordinates
      qends[nl1_bu] <- min(qend1, qend2)
      qstarts[nl1_bu] <- max(qstart1, qstart2)

      # target coordinates
      tends[nl1_bu] <- max(tend1, tend2)
      tstarts[nl1_bu] <- min(tstart1, tstart2)
    }

    # nmatches and alens
    nmatches[nl1_bu] <- nmatches[nl1] + nmatches[nl2]
    alens[nl1_bu] <- alens[nl1] + alens[nl2]
  
    
    count <- count + 1
  }

  # Re-construct a paffile out of the columns
  paffile = data.frame(qname = qnames,
                       qlen = paffile$qlen,
                       qstart = qstarts,
                       qend = qends,
                       strand = strands,
                       tname = paffile$tname,
                       tlen = paffile$tlen,
                       tstart = tstarts,
                       tend = tends,
                       nmatch = nmatches,
                       alen = alens,
                       mapq = paffile$mapq
                       )
  
  


  paffile <- paffile[!(row.names(paffile) %in% pairs$col), ]

  # Filter
  paffile <- paffile[dplyr::between((abs(as.numeric(paffile$tend) - as.numeric(paffile$tstart)) /
    abs(as.numeric(paffile$qend) - as.numeric(paffile$qstart))), 0.5, 2), ]

  # Print out the longest t length of the paffile
  max_tlen = max(abs(as.numeric(paffile$tend) - as.numeric(paffile$tstart)))
  # Print out the number of alingments in the paffile
  # print(paste0("Number of alignments: ", dim(paffile)[1]))
  # Print this result
  # print(paste0("Max tlen: ", max_tlen))
  # print('done')
  return(paffile)
}



#' Melt Alignment Chunks Back Together
#'
#' @description This function is a crucial part of the process of identifying synteny blocks (SDs)
#' and is the last step of the chunk-minimap2-edit-melt workflow. It takes the chunked PAF file
#' produced earlier in the pipeline and "melts" the alignment pieces back together.
#' Note that improvements may be needed in the future, as the melting process is not
#' entirely straightforward.
#'
#' @param inpaf_link A file path or link to the input chunked PAF file.
#' @param outpaf_link A file path or link to the output "molten" PAF file.
#' @param inpaf_df An optional input PAF data.frame.
#' @param save_outpaf Logical value whether to save the output PAF file (default: TRUE).
#' @param n_quadrants_per_axis The number of quadrants per axis (determines chunking).
#' @param second_run Logical value whether this is the second run of the function (default: FALSE).
#' @param inparam_chunklen The input parameter for chunk length in bp.
#' @param inparam_compression The input parameter for compression in bp.
#' @author Wolfram Hoeps
#' @export
compress_paf_fnct <-
  function(inpaf_link = NULL,
           outpaf_link = NULL,
           inpaf_df = NULL,
           save_outpaf = T,
           n_quadrants_per_axis = NULL,
           second_run = F,
           inparam_chunklen = NULL,
           inparam_compression = NULL) {
    if (is.null(inpaf_df)) {
      # Read and col-annotate the input paf
      inpaf <- read.table(inpaf_link, sep = "\t")
      # Remove redundant rows
      inpaf <- unique(inpaf)

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
      colnames(inpaf)[1:length(colnames_paf)] <- colnames_paf

      inpaf <- transform(
        inpaf,
        qend = ifelse(strand == "-", qstart, qend),
        qstart = ifelse(strand == "-", qend, qstart)
      )

      inpaf <- inpaf[order(inpaf$qstart), ]
      # inpaf = inpaf[(inpaf$tstart > 809000) & (inpaf$qstart > 809000),]
      # # debug input:
      # inpaf = inpaf[abs((inpaf$tstart + inpaf$qstart) -  950000) < 50000,]
      # inpaf = inpaf[inpaf$alen > 5000,]
    } else {
      inpaf <- inpaf_df[order(inpaf_df$qstart), ]
      # Rename rows
      row.names(inpaf) <- 1:dim(inpaf)[1]
    }

    # In paf, start is always smaller than end. For our slope etc calculation, it will be easier
    # to change the order of start and end, if the orientation of the alignment is negative.
    # inpaf = transform(
    #   inpaf,
    #   qend = ifelse(strand == '-', qstart, qend),
    #   qstart = ifelse(strand == '-', qend, qstart)
    # )
    # For safety: sort entries by tstart. Reset row names so they start at 1.
    inpaf <- inpaf[order(inpaf$tstart), ]
    rownames(inpaf) <- NULL

    # We identify rowpairs now.
    # To speed this up, we cut the alignment into chunks, find
    # rowpairs within them and later merge back together.
    range <- c(max(inpaf$tend), max(inpaf$qend))
    if (dim(inpaf)[1] < 15000) {
      minlen_for_merge <- 0
    } else {
      minlen_for_merge <- inparam_chunklen * 0.9
    }

    # Determine number of quadrants

    if (is.null(n_quadrants_per_axis)) {
      if (range[1] < 50000) {
        n_quadrants_per_axis <- 2
      } else if (range[1] >= 50000) {
        n_quadrants_per_axis <- 10
      }
    }

    inpaf_short_alns <- inpaf[inpaf$alen <= minlen_for_merge, ]
    inpaf <- inpaf[inpaf$alen > minlen_for_merge, ]
    row.names(inpaf) <- NULL
    tsteps <- round(seq(1, range[1], length.out = n_quadrants_per_axis))
    tstepsize <- unique(diff(tsteps))[1]
    if (n_quadrants_per_axis == 1) {
      tstepsize <- 1e10
    }
    rowpairs <- data.frame()
    count <- 0

    tstarts = inpaf$tstart
    tends = inpaf$tend
    qstarts = inpaf$qstart
    qends = inpaf$qend

    
    for (tstep in tsteps[1:length(tsteps) - 1]) {
        # Input: we take any alignment that touches our box.
        
        #print(p + geom_rect(aes(xmin=tstep, xmax=tstep+tstepsize, ymin=qstep,ymax=qstep+qstepsize), color='green', alpha=0.1))
        inpaf_q <- inpaf[((tstarts >= tstep) & (tstarts <= tstep + tstepsize) | (tends  >= tstep) & (tends <= tstep + tstepsize)),] 

        if (nrow(inpaf_q) > 1){
          rowpairs <- rbind(rowpairs, merge_paf_entries_intraloop(inpaf_q, second_run, inparam_chunklen, inparam_compression))
          #rowpairs <- merge_paf_entries_intraloop(inpaf_q, second_run, inparam_chunklen, inparam_compression)

        } 
        count <- count + 1
        # print(paste0('Merging step ', count, ' out of ', n_steps))
    }
    
    #print('pairs found!')
    # Sort once again, by row.
    rowpairs <- rowpairs[order(rowpairs$col), ]
    row.names(rowpairs) <- NULL

    # Cleanup
    rowpairs <- unique(rowpairs)

    # Remove rows that want to pair with themselves (should only appear with a fixed tolerance bp,
    # not it the tolerance bp is a fraction of the alignment length)
    rowpairs <- rowpairs[rowpairs$row != rowpairs$col, ]

    # For each pair, show the number of matched bases by the unity.
    # We will use this metric to decide a 'winning' pair if any vectors
    # want to pair with multiple other vectors.
    rowpairs$combined_matchlen <- inpaf[rowpairs$row, ]$nmatch + inpaf[rowpairs$col, ]$nmatch

    # Remove very short stuff
    rowpairs <- rowpairs[rowpairs$combined_matchlen > (inparam_chunklen / 5), ]

    "%>%" <- magrittr::"%>%"

    # Cut down redundant pairs
    rowpairs_singular <- as.data.frame(
      rowpairs %>%
        dplyr::group_by(row) %>% dplyr::top_n(1, combined_matchlen) %>%
        dplyr::group_by(col) %>% dplyr::top_n(1, combined_matchlen)
    )


    # Go through each pair, make the merge. We go through the lines backwards,
    # so that previous merges don't disturb later ones.
    # old_len = sum(inpaf$tend - inpaf$tstart) + sum(inpaf$qend - inpaf$qstart)
    count <- 0
    if (dim(rowpairs_singular)[1] > 0) {
      inpaf <- merge_rows(inpaf, rowpairs_singular)
    }

    # Re-merge with old thing
    inpaf <- rbind(inpaf, inpaf_short_alns)

    # Change orientation BACK. The 'merge_rows' functino is written in that way.
    # inpaf = transform(
    #   inpaf,
    #   qend = ifelse(strand == '-', qstart, qend),
    #   qstart = ifelse(strand == '-', qend, qstart)
    # )

    #print(paste0("PAF compressed to ", dim(inpaf)[1], " alignments."))
    # Save
    if (save_outpaf) {
      inpaf <- transform(
        inpaf,
        qend = ifelse(strand == "-", qstart, qend),
        qstart = ifelse(strand == "-", qend, qstart)
      )
      write.table(
        inpaf,
        file = outpaf_link,
        quote = F,
        col.names = F,
        row.names = F,
        sep = "\t"
      )
    } else {
      return(inpaf)
    }
  }
