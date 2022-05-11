#!/usr/local/bin/Rscript

#' Helperfunction (1/1) of compress_paf_fnct
#'
#' @description Melt two alignment entries into one.
#'
#' @param paffile [data frame] loaded paffile with colnames
#' @param nl1 [character/link] line number of first-of-pair
#' @param nl2 [character/link] line number of second-of-pair
#' @return paffile with one row less (bc a pair has been merged)
#'
#' @author Wolfram Höps
#' @export
merge_rows <- function(paffile, nl1, nl2) {
  # paffile = inpaf
  # nl1 = rowpairs[1,1]
  # nl2 = rowpairs[1,2]
  # query name
  paffile[nl1,]$qname = paste0(sub("_.*", "", paffile[nl1,]$qname),
                               "_",
                               paffile[nl1,]$qstart,
                               "-",
                               paffile[nl2,]$qend - 1)
  # query coordinates
  paffile[nl1,]$qend = paffile[nl2,]$qend
  
  # target coordinates
  if (paffile[nl1,]$strand == '+') {
    paffile[nl1,]$tend = paffile[nl2,]$tend
  } else if (paffile[nl1,]$strand == '-') {
    paffile[nl1,]$tstart = paffile[nl2,]$tstart
  }
  # nmatch
  paffile[nl1,]$nmatch = paffile[nl1,]$nmatch + paffile[nl2,]$nmatch
  # alen
  paffile[nl1,]$alen = paffile[nl1,]$alen + paffile[nl2,]$alen
  
  # Remove 2nd line
  paffile = paffile[-nl2,]
  
  return(paffile)
  
}






#' Melt the alignment pieces back together after chunkifying them earlier.
#'
#' @description This is a crucial module for finding the SDs later, and poses
#' the last part of the chunk-minimap2-edit-melt workflow. This might need more
#' work in the future, because melting is not so straight forward.
#'
#' @param inpaf_link [character/link] link to the chunked paf
#' @param outpaf_link [character/link] link to the molten output paf
#' @param quadrantsize [numeric] size of alignments being merged in one piece.
#' Set to larger values (e.g. 250k) for very large alignments (>5Mb).
#'
#' @author Wolfram Höps
#' @export
compress_paf_fnct <-
  function(inpaf_link = NULL,
           outpaf_link = NULL,
           inpaf_df = NULL,
           save_outpaf = T,
           quadrantsize = 100000, 
           second_run = F, 
           inparam_chunklen = NULL,
           inparam_compression = NULL) {
    library(dplyr)
    debug = F
    if (debug) {
      inpaf_link = "/Users/hoeps/PhD/projects/nahrcall/nahrchainer/teppich/res_500k/paf/run_1024_1024_0.90_+_chunked.paf.chunk"
      inpaf_link = "/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/res/outpaf.awked"
      inpaf_link = 'blub.awked'
      inpaf_link = "/Users/hoeps/phd/projects/nahrcall/nahrchainer/seqbuilder/res/outpaf13.awked"
    }
    #inpaf_link = "/Users/hoeps/phd/projects/nahrcall/nahrchainer/seqbuilder/res/outpaf136.awked"
    
    if (is.null(inpaf_df)) {
      # Read and col-annotate the input paf
      inpaf = read.table(inpaf_link, sep = '\t')
      # Remove redundant rows
      inpaf = unique(inpaf)
      
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
      colnames(inpaf)[1:length(colnames_paf)] = colnames_paf
      
      #inpaf = inpaf[(inpaf$tstart > 809000) & (inpaf$qstart > 809000),]
      # # debug input:
      # inpaf = inpaf[abs((inpaf$tstart + inpaf$qstart) -  950000) < 50000,]
      # inpaf = inpaf[inpaf$alen > 5000,]
      # For safety: sort entries by qstart. Reset row names so they start at 1.
      inpaf = inpaf[order(inpaf$qstart),]
      rownames(inpaf) <- NULL
    } else {
      inpaf = inpaf_df[order(inpaf_df$qstart),]
      # Rename rows 
      row.names(inpaf) = 1:dim(inpaf)[1]
    }

    
    # Thin if there are just *too* many alignments
    # browser()
    # ggplot2::ggplot(inpaf) + ggplot2::geom_segment(ggplot2::aes(x=qstart, xend=qend, y=tstart,yend=tend))
    alen_min = 0
    while (dim(inpaf)[1] > 50000){
      alen_min = alen_min + 100
      inpaf = inpaf[inpaf$alen >= alen_min,]
      print(paste0('Trimming length: ', alen_min))
    }
    # Compression here?
    # compression = 1000
    # inpaf[,c('qstart','qend','tstart','tend')] = round(data.frame(inpaf[,c('qstart','qend','tstart','tend')] / compression),0) * compression
    # inpaf = inpaf[(inpaf$qend != inpaf$qstart) & (inpaf$tend != inpaf$tstart) ,]
    
    # We identify rowpairs now.
    # To speed this up, we cut the alignment into chunks, find
    # rowpairs within them and later merge back together.
    range = c(max(inpaf$tend), max(inpaf$qend))
    
    tsteps = c(seq(1, range[1], quadrantsize), max(inpaf$tend))
    qsteps = c(seq(1, range[2], quadrantsize), max(inpaf$qend))
    n_steps = length(tsteps) * length(qsteps)
    rowpairs = data.frame()
    count = 0
    #browser()
    for (tstep in tsteps[1:length(tsteps)-1]) {
      for (qstep in qsteps[1:length(qsteps)-1]) {
        # Input: we take any alignment that touches our box.
        inpaf_q = inpaf[(
          (inpaf$tstart >= tstep) & (inpaf$tstart <= tstep + quadrantsize) &
            (inpaf$qstart >= qstep) &
            (inpaf$qstart <= qstep + quadrantsize)
        ) |
          (
            (inpaf$tend >= tstep) & (inpaf$tend <= tstep + quadrantsize) &
              (inpaf$qend >= qstep) &
              (inpaf$qstart <= qstep + quadrantsize)
          )   ,]
        rowpairs = rbind(rowpairs, merge_paf_entries_intraloop(inpaf_q, second_run, inparam_chunklen, inparam_compression))
        count = count + 1
        #print(paste0('Merging step ', count, ' out of ', n_steps))
      }
    }
    # Sort once again, by row.
    rowpairs = rowpairs[order(rowpairs$col),]
    row.names(rowpairs) <- NULL
    
    # Cleanup
    rowpairs = unique(rowpairs)
    
    # Remove rows that want to pair with themselves (should only appear with a fixed tolerance bp,
    # not it the tolerance bp is a fraction of the alignment length)
    rowpairs = rowpairs[rowpairs$row != rowpairs$col,]
    
    # For each pair, show the number of matched bases by the unity.
    # We will use this metric to decide a 'winning' pair if any vectors
    # want to pair with multiple other vectors.
    rowpairs$combined_matchlen = inpaf$nmatch[rowpairs$row] + inpaf$nmatch[rowpairs$col]
    # rowpairs$tdist1 = inpaf$tend[rowpairs$row] - inpaf$tstart[rowpairs$col]
    # rowpairs$tdist2 = inpaf$tstart[rowpairs$row] - inpaf$tend[rowpairs$col]
    # rowpairs$qdist1 = inpaf$qend[rowpairs$row] - inpaf$qstart[rowpairs$col]
    # rowpairs$qdist2 = inpaf$qstart[rowpairs$row] - inpaf$qend[rowpairs$col]
    
    # Remove very short stuff
    rowpairs = rowpairs[rowpairs$combined_matchlen > (inparam_chunklen/5),]

    # Cut down redundant pairs
    rowpairs_singular = as.data.frame(
      rowpairs %>%
        group_by(row) %>% top_n(1, combined_matchlen) %>%
        group_by(col) %>% top_n(1, combined_matchlen)
    )
    
    # Go through each pair, make the merge. We go through the lines backwards,
    # so that previous merges don't disturb later ones.
    # old_len = sum(inpaf$tend - inpaf$tstart) + sum(inpaf$qend - inpaf$qstart)
    if (dim(rowpairs_singular)[1] > 0) {
      for (nrow in dim(rowpairs_singular)[1]:1) {
        inpaf = merge_rows(inpaf, rowpairs_singular[nrow, 1], rowpairs_singular[nrow, 2])
        # new_len = (sum(inpaf$tend - inpaf$tstart) + sum(inpaf$qend - inpaf$qstart))
        # if (abs(old_len - new_len) > 100000){
        #   browser()
        # }
        #print(old_len - new_len)
      }
    }
    
    print(paste0('PAF compressed to ', dim(inpaf)[1], ' alignments.'))
    # Save
    if (save_outpaf) {
      write.table(
        inpaf,
        file = outpaf_link,
        quote = F,
        col.names = F,
        row.names = F,
        sep = '\t'
      )
    }
    else {
      return(inpaf)
    }
  }


# # runs only when script is run by itself
# if (sys.nframe() == 0) {
#   # Define input
#   inpaf_link = commandArgs(trailingOnly = TRUE)[1]
#   outpaf_link = commandArgs(trailingOnly = TRUE)[2]
#   
#   compress_paf_fnct(inpaf_link, outpaf_link)
# }
