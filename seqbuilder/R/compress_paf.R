#!/usr/local/bin/Rscript

merge_rows <- function(paffile, nl1, nl2){
  
  # query name
  paffile[nl1,]$qname = paste0(
    sub("_.*", "", paffile[nl1,]$qname), 
    "_", 
    paffile[nl1,]$qstart, 
    "-",
    paffile[nl2,]$qend - 1
  )
  # query coordinates
  paffile[nl1,]$qend = paffile[nl2,]$qend
  
  # target coordinates
  if (paffile[nl1,]$strand == '+'){
    paffile[nl1,]$tend = paffile[nl2,]$tend
  } else if (paffile[nl1,]$strand == '-'){
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

compress_paf_fnct <- function(inpaf_link, outpaf_link){
  inpaf = read.table(inpaf_link, sep='\t')
  colnames_paf = c('qname','qlen','qstart','qend',
                   'strand','tname','tlen','tstart',
                   'tend','nmatch','alen','mapq')
  
  colnames(inpaf)[1:length(colnames_paf)] = colnames_paf
  
  # For safety: sort entries by qstart
  inpaf = inpaf[order(inpaf$qstart),]
  
  # Identify alignments that border each other: same strand, and end of of is the start
  # of the other.

  tolerance_bp = 10
  rowpairs = which( ( abs(outer(inpaf$qend, inpaf$qstart, '-')) < tolerance_bp) &
                      (
                        (abs(outer(inpaf$tend, inpaf$tstart, '-')) < tolerance_bp) |
                        (abs(outer(inpaf$tstart, inpaf$tend, '-')) < tolerance_bp)
                      ) &
                      (outer(inpaf$strand, inpaf$strand, '==')),
                    arr.ind = T)
  # Go through each pair, make the merge. We go through the lines backwards,
  # so that previous merges don't disturb later ones.
  for (nrow in dim(rowpairs)[1]:1){
    inpaf = merge_rows(inpaf, rowpairs[nrow, 1], rowpairs[nrow, 2])
  }  
  # Save
  write.table(inpaf, file=outpaf_link, quote = F, col.names = F, row.names = F, sep='\t')
}


# runs only when script is run by itself
if (sys.nframe() == 0){
  
  # Define input
  inpaf_link = commandArgs(trailingOnly=TRUE)[1]
  outpaf_link = commandArgs(trailingOnly=TRUE)[2]

  compress_paf_fnct(inpaf_link, outpaf_link)
}




