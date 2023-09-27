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
  for (i in dim(pairs)[1]:1) {
    if (count %% 100 == 0) {
      print(paste0("Merging pair ", count, " out of ", dim(pairs)[1]))
    }
    nl1 <- (pairs[i, 1])
    nl2 <- (pairs[i, 2])
    paf_work <- paffile[c(nl1, nl2), ]
    if (paf_work[1, "strand"] == "+") {
      paf_work[1, ]$qname <- paste0(
        # HERE is the offending line!!
        # sub("_.*", "", paf_work[1, ]$qname) <- this line removes everything after the first '_'. Actually i want to remove everything after the *last* instance of '_'. I dont know beforehand if there are 1, 2 or more "_". Therefore, the nnew line is:
        sub("_[^_]*$", "", paf_work[1, ]$qname),
        "_",
        paf_work[1, ]$qstart,
        "-",
        paf_work[2, ]$qend - 1
      )

      paf_work_bu <- paf_work
      # query coordinates
      paf_work[1, ]$qend <- max(paf_work_bu[1, ]$qend, paf_work_bu[2, ]$qend)
      paf_work[1, ]$qstart <- min(paf_work_bu[1, ]$qstart, paf_work_bu[2, ]$qstart)

      # target coordinates
      paf_work[1, ]$tend <- max(paf_work_bu[1, ]$tend, paf_work_bu[2, ]$tend)
      paf_work[1, ]$tstart <- min(paf_work_bu[1, ]$tstart, paf_work_bu[2, ]$tstart)
    } else if (paf_work[1, "strand"] == "-") {
      # Sort by t, because this is increasing
      paf_work <- paf_work[order(paf_work$tstart), ]
      paf_work[1, ]$qname <- paste0(
        sub("_[^_]*$", "", paf_work[1, ]$qname),
        "_",
        paf_work[2, ]$qend,
        "-",
        paf_work[2, ]$qstart
      )
      paf_work_bu <- paf_work
      # query coordinates
      paf_work[1, ]$qend <- min(paf_work_bu[1, ]$qend, paf_work_bu[2, ]$qend)
      paf_work[1, ]$qstart <- max(paf_work_bu[1, ]$qstart, paf_work_bu[2, ]$qstart)

      # target coordinates
      paf_work[1, ]$tend <- max(paf_work_bu[1, ]$tend, paf_work_bu[2, ]$tend)
      paf_work[1, ]$tstart <- min(paf_work_bu[1, ]$tstart, paf_work_bu[2, ]$tstart)
    }


    # nmatch
    paf_work[1, ]$nmatch <- paf_work[1, ]$nmatch + paf_work[2, ]$nmatch
    # alen
    paf_work[1, ]$alen <- paf_work[1, ]$alen + paf_work[2, ]$alen

    # Now process the main paf thing
    paffile[nl1, ] <- paf_work[1, ]

    count <- count + 1
  }

  paffile <- paffile[!(row.names(paffile) %in% pairs$col), ]

  # Filter
  paffile <- paffile[dplyr::between((abs(paffile$tend - paffile$tstart) /
    abs(paffile$qend - paffile$qstart)), 0.5, 2), ]

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
    # For safety: sort entries by qstart. Reset row names so they start at 1.
    inpaf <- inpaf[order(inpaf$qstart), ]
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
    qsteps <- round(seq(1, range[2], length.out = n_quadrants_per_axis))
    tstepsize <- unique(diff(tsteps))[1]
    qstepsize <- unique(diff(qsteps))[1]
    if (n_quadrants_per_axis == 1) {
      tstepsize <- 1e10
      qstepsize <- 1e10
    }
    n_steps <- length(tsteps) * length(qsteps)
    rowpairs <- data.frame()
    count <- 0

    for (tstep in tsteps[1:length(tsteps) - 1]) {
      for (qstep in qsteps[1:length(qsteps) - 1]) {
        # Input: we take any alignment that touches our box.



        inpaf_q <- inpaf[(
          (inpaf$tstart >= tstep) & (inpaf$tstart <= tstep + tstepsize) &
            (inpaf$qstart >= qstep) &
            (inpaf$qstart <= qstep + qstepsize)
        ) |
          (
            (inpaf$tend >= tstep) & (inpaf$tend <= tstep + tstepsize) &
              (inpaf$qend >= qstep) &
              (inpaf$qend <= qstep + qstepsize)
          ), ]

        # print(ggplot() + geom_segment(data=inpaf, aes(x=tstart, xend=tend, y=qstart, yend=qend)) +
        #   geom_segment(data=inpaf_q, aes(x=tstart, xend=tend, y=qstart, yend=qend), color='red') +
        #   geom_rect(aes(xmin=tstep, xmax=tstep+tstepsize, ymin=qstep, ymax=qstep+qstepsize), alpha=0.5))


        rowpairs <- rbind(rowpairs, merge_paf_entries_intraloop(inpaf_q, second_run, inparam_chunklen, inparam_compression))
        count <- count + 1
        # print(paste0('Merging step ', count, ' out of ', n_steps))
      }
    }
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

    print(paste0("PAF compressed to ", dim(inpaf)[1], " alignments."))
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
