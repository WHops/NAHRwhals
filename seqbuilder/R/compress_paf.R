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
#' @rdname alignment
#' @export
merge_rows <- function(paffile, nl1, nl2){
  
  # paffile = inpaf
  # nl1 = rowpairs[1,1]
  # nl2 = rowpairs[1,2]
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

#' Melt the alignment pieces back together after chunkifying them earlier.   
#' 
#' @description This is a crucial module for finding the SDs later, and poses
#' the last part of the chunk-minimap2-edit-melt workflow. This might need more 
#' work in the future, because melting is not so straight forward. 
#' 
#' @param inpaf_link [character/link] link to the chunked paf
#' @param outpaf_link [character/link] link to the molten output paf
#' 
#' @author Wolfram Höps
#' @rdname alignment
#' @export
#' 
#' 
compress_paf_fnct <- function(inpaf_link, outpaf_link){
  
  

  
  debug = F
  if (debug){
  inpaf_link = "/Users/hoeps/PhD/projects/nahrcall/nahrchainer/teppich/res_500k/paf/run_1024_1024_0.90_+_chunked.paf.chunk"
  inpaf_link = "/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/res/outpaf.awked"
  inpaf_link = 'blub.awked'
  } 
  
  # Read and col-annotate the input paf
  inpaf = read.table(inpaf_link, sep='\t')
  
  colnames_paf = c('qname','qlen','qstart','qend',
                   'strand','tname','tlen','tstart',
                   'tend','nmatch','alen','mapq')
  colnames(inpaf)[1:length(colnames_paf)] = colnames_paf
  
  # For safety: sort entries by qstart
  inpaf = inpaf[order(inpaf$qstart),]
  
  # Identify alignments that border each other: same strand, and end of of is the start
  # of the other. With some tolerance
  tolerance_bp = 10
  rowpairs = data.frame(which( ( abs(outer(inpaf$qend, inpaf$qstart, '-')) < tolerance_bp) & # Take only one half of the minus matrix so pairs dont appear twice.
                     (
                       (abs(outer(inpaf$tend, inpaf$tstart, '-')) < tolerance_bp) |
                       (abs(outer(inpaf$tstart, inpaf$tend, '-')) < tolerance_bp)
                     ) &
                     (outer(inpaf$strand, inpaf$strand, '==')),
                   arr.ind = T))

  
  # Plot those pairs?
  #ggplot2::ggplot(rowpairs) + ggplot2::geom_point(ggplot2::aes(x=row, y=col))
  
  # Resolve multi-pairs 
  # By taking blobs of pairs and keeping only the top two
  # TODO. 
  #n_occur <- data.frame(table(rowpairs$row))
  #n_occur[n_occur$Freq > 1,]
  #rowpairs[rowpairs$row %in% n_occur$Var1[n_occur$Freq > 1],]
  
  
  # Go through each pair, make the merge. We go through the lines backwards,
  # so that previous merges don't disturb later ones.
  # TODO: what if one alignment borders multiple others? 
  if (dim(rowpairs)[1] > 0){
    for (nrow in dim(rowpairs)[1]:1){
      inpaf = merge_rows(inpaf, rowpairs[nrow, 1], rowpairs[nrow, 2])
    }  
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




