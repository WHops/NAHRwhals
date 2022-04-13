#' Tiny undocumented helperfunction.
#' @author Wolfram Höps
#' @export
merge_paf_entries_intraloop <- function(inpaf, second_run=F, inparam_chunklen=NULL, inparam_compression=NULL) {
  inpaf_rownames = row.names(inpaf)
  
  # For safety: sort entries by qstart. Reset row names so they start at 1.
  inpaf = inpaf[order(inpaf$qstart),]
  rownames(inpaf) <- NULL
  #browser()
  # We consider alignments as 'potential neighbours' if their distance in any direction
  # (+-x, +-y) is less than 5% of their alignment length.
  #tolerance_bp = 10 #
  if (second_run){
    # Not setting to 1 because i'm scared of edge cases.. But should be one. 
    #tolerance_bp = 2
    
    # Fixing a bug here which caused problems with wrong stuff being merged,
    # If their start AND end points are both in reach of one point.
    chunklen = min(inpaf$qlen)
    tolerance_bp_unlimited = 0.5 * inparam_compression
    tolerance_bp = min(tolerance_bp_unlimited, chunklen * 0.9)
    print(tolerance_bp)
    
  } else {
    # On second thought, let's have this cap on the length of the input chunklen. 
    tolerance_bp = 0.05 * inparam_chunklen
    
    # # I want to cap this, so e.g. a 1MB stretch can not just acquire everything in 50kbp distance. 
    # tolerance_bp_unlimited = 0.05 * (outer(inpaf$alen, inpaf$alen, '+') / 2)
    # tolerance_bp = min(tolerance_bp_unlimited, )
  }
  # Identify alignments that border each other: same strand, and end of of is the start
  # of the other. With some tolerance
  # Take only one half of the minus matrix so pairs dont appear twice.
  inpaf_t2 = inpaf
  inpaf_t =  transform(
    inpaf_t2,
    qend = ifelse(strand == '-', qstart, qend),
    qstart = ifelse(strand == '-', qend, qstart)
  )
  
  # Quick and dirty fix for problem of matches if both 
  # start and end are inside of the threshold. 
  if (tolerance_bp > 100){
    inpaf_t = inpaf_t[(inpaf_t$alen > tolerance_bp),]
  }
  rowpairs = data.frame(which(
    
    # only first half of matrix
    #(outer(inpaf_t$tend, inpaf_t$tstart, '-') <= tolerance_bp) &
      
    (((abs(outer(inpaf_t$tend, inpaf_t$tstart, '-')) < tolerance_bp) &
      (abs(outer(inpaf_t$qend, inpaf_t$qstart, '-')) < tolerance_bp)) | 
         
    ((abs(outer(inpaf_t$tstart, inpaf_t$tend, '-')) < tolerance_bp) &
     (abs(outer(inpaf_t$qstart, inpaf_t$qend, '-')) < tolerance_bp))) &  

    # Same directionality
    (outer(inpaf_t$strand, inpaf_t$strand, '==')),
    arr.ind = T))
  
  rowpairs = rowpairs[rowpairs$col > rowpairs$row,]

  # rowpairs = data.frame(which(
  #   
  #   
  #   (((abs(outer(inpaf$tend, inpaf$tstart, '-')) < tolerance_bp) &
  #      (abs(outer(inpaf$qend, inpaf$qstart, '-')) < tolerance_bp)) | 
  #     
  #      
  #        ((abs(outer(inpaf$tstart, inpaf$tend, '-')) < tolerance_bp) &
  #      
  #        (abs(outer(inpaf$qstart, inpaf$qend, '-')) < tolerance_bp))) &
  #     
  #     (outer(inpaf$strand, inpaf$strand, '==')),
  #   arr.ind = T))
  
  rowpairs$row = as.numeric(inpaf_rownames[rowpairs$row])
  rowpairs$col = as.numeric(inpaf_rownames[rowpairs$col])
  
  
  return(rowpairs)
}
