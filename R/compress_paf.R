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
merge_rows <- function(paffile, pairs) {
  
  # paf_mergers = paffile[row.names(paffile) %in% unique(c(rowpairs$row, rowpairs$col)),]
  # paf_nonmergers = paffile[!(row.names(paffile) %in% unique(c(rowpairs$row, rowpairs$col))),]
  # 
  # 
  
  count = 0
  for (i in dim(pairs)[1]:1){
    if (count %% 100 == 0){
      print(paste0('Merging pair ', count, ' out of ', dim(pairs)[1]))
    }
    nl1 = (pairs[i,1])
    nl2 = (pairs[i,2])
    paf_work = paffile[c(nl1,nl2),]
    # paffile = inpaf
    # nl1 = rowpairs[1,1]
    # nl2 = rowpairs[1,2]
    # query name
    
    if (paf_work[1,'strand'] == '+'){
      paf_work[1,]$qname = paste0(sub("_.*", "", paf_work[1,]$qname),
                                   "_",
                                  paf_work[1,]$qstart,
                                   "-",
                                  paf_work[2,]$qend - 1)



    } else if (paf_work[1,'strand'] == '-'){
      # Sort by t, because this is increasing
      paf_work = paf_work[order(paf_work$tstart),]
      paf_work[1,]$qname = paste0(sub("_.*", "", paf_work[1,]$qname),
                                  "_",
                                  paf_work[2,]$qend,
                                  "-",
                                  paf_work[2,]$qstart)

    }
    
    # query coordinates
    paf_work[1,]$qend = paf_work[2,]$qend
    
    # target coordinates
    paf_work[1,]$tend = paf_work[2,]$tend
    
    # nmatch
    paf_work[1,]$nmatch = paf_work[1,]$nmatch + paf_work[2,]$nmatch
    # alen
    paf_work[1,]$alen = paf_work[1,]$alen + paf_work[2,]$alen
    
    # Now process the main paf thing
    paffile[nl1,] = paf_work[1,]
    
    count = count+1
  }
  
  paffile = paffile[!(row.names(paffile) %in% pairs$col),]
  
  #Filter
  paffile = paffile[  between((abs(paffile$tend - paffile$tstart) / 
                          abs(paffile$qend - paffile$qstart)), 0.5, 2),  ]
  
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
#' Set to larger values (e.g. 250k) for very large alignments (>5Mb).
#'
#' @author Wolfram Höps
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

    } else {
      inpaf = inpaf_df[order(inpaf_df$qstart),]
      # Rename rows 
      row.names(inpaf) = 1:dim(inpaf)[1]
    }

    # In paf, start is always smaller than end. For our slope etc calculation, it will be easier
    # to change the order of start and end, if the orientation of the alignment is negative.
    inpaf = transform(
      inpaf,
      qend = ifelse(strand == '-', qstart, qend),
      qstart = ifelse(strand == '-', qend, qstart)
    )
    # For safety: sort entries by qstart. Reset row names so they start at 1.
    inpaf = inpaf[order(inpaf$qstart),]
    rownames(inpaf) <- NULL
    # Thin if there are just *too* many alignments
    # browser()
    # ggplot2::ggplot(inpaf) + ggplot2::geom_segment(ggplot2::aes(x=qstart, xend=qend, y=tstart,yend=tend))
    # alen_min = 0
    # while (dim(inpaf)[1] > 10000){
    #   alen_min = alen_min + 500
    #   inpaf = inpaf[inpaf$alen >= alen_min,]
    #   print(paste0('Trimming length: ', alen_min))
    # }
    # Compression here?
    # compression = 1000
    # inpaf[,c('qstart','qend','tstart','tend')] = round(data.frame(inpaf[,c('qstart','qend','tstart','tend')] / compression),0) * compression
    # inpaf = inpaf[(inpaf$qend != inpaf$qstart) & (inpaf$tend != inpaf$tstart) ,]
    
    # We identify rowpairs now.
    # To speed this up, we cut the alignment into chunks, find
    # rowpairs within them and later merge back together.
    range = c(max(inpaf$tend), max(inpaf$qend))
    if (dim(inpaf)[1] < 15000){
      minlen_for_merge = 0
    } else {
      minlen_for_merge = inparam_chunklen*0.9
    }
    
    # Determine number of quadrants
    
    if (is.null(n_quadrants_per_axis)){
      if (range[1] < 50000){
        n_quadrants_per_axis = 2
      } else if (range[1] >= 50000) {
        n_quadrants_per_axis = 10
      }
    }
    
    inpaf_short_alns = inpaf[inpaf$alen <= minlen_for_merge,] 
    inpaf = inpaf[inpaf$alen > minlen_for_merge,] 
    row.names(inpaf) = NULL
    tsteps = round(seq(1, range[1], length.out=n_quadrants_per_axis))
    qsteps = round(seq(1, range[2], length.out=n_quadrants_per_axis))
    tstepsize = unique(diff(tsteps))[1]
    qstepsize = unique(diff(qsteps))[1]
    if (n_quadrants_per_axis == 1){
      tstepsize = 1e10
      qstepsize = 1e10
    }
    n_steps = length(tsteps) * length(qsteps)
    rowpairs = data.frame()
    count = 0
    for (tstep in tsteps[1:length(tsteps)-1]) {
      for (qstep in qsteps[1:length(qsteps)-1]) {
        # Input: we take any alignment that touches our box.


         inpaf_q = inpaf[(
          (inpaf$tstart >= tstep) & (inpaf$tstart <= tstep + tstepsize) &
            (inpaf$qstart >= qstep) &
            (inpaf$qstart <= qstep + qstepsize)
        ) |
          (
            (inpaf$tend >= tstep) & (inpaf$tend <= tstep + tstepsize) &
              (inpaf$qend >= qstep) &
              (inpaf$qend <= qstep + qstepsize)
          )   ,] 
         

        rowpairs = rbind(rowpairs,merge_paf_entries_intraloop(inpaf_q, second_run, inparam_chunklen, inparam_compression))
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
    rowpairs$combined_matchlen = inpaf[rowpairs$row,]$nmatch + inpaf[rowpairs$col,]$nmatch
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
    count = 0
    if (dim(rowpairs_singular)[1] > 0) {
      
        inpaf = merge_rows(inpaf, rowpairs_singular)
    }
    
    # Re-merge with old thing
    inpaf = rbind(inpaf, inpaf_short_alns)
    
    # Change orientation BACK. The 'merge_rows' functino is written in that way. 
    inpaf = transform(
      inpaf,
      qend = ifelse(strand == '-', qstart, qend),
      qstart = ifelse(strand == '-', qend, qstart)
    )
    
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
