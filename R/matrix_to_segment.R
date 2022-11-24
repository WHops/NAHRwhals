#' helperfunction
#' @export
absmin <- function(x) { x[which.min( abs(x) )]}


#' helperfunction
#' @export
gm_to_segvec <- function(gridmatrix, groundlen = 10000){
  

  
  # IDS
  ids = cumsum(as.numeric(colnames(gridmatrix)) / groundlen)
  
  # Matchsizes
  sizes = data.frame(rname = ids, size = as.numeric(colnames(gridmatrix)))
  row.names(sizes) = sizes$rname
  
  # Turn into increasing numbers
  m2 = t(t(sign(gridmatrix)) * ids)
  
  
  
  # Make empty output vector
  outvec = rep(0,dim(m2)[1])
  
  # For each row, determine the index with the lowest absolute number. This is 
  # The index we go for. 
  for (matrow in 1:dim(m2)[1]){
    
    # For the segment ID, we use its lowest_index,
    # I.e. the index of the first time this was observed 
    # on the regerence. 
    segment_ID =  absmin(m2[matrow, ][m2[matrow, ] != 0])
    m2[matrow, ] = abs(sign(m2[matrow, ])) * segment_ID
    outvec[matrow] = segment_ID
  }
  

  ns = (sizes[as.character(abs(outvec)),]$size * sign(outvec)) / groundlen
  for (n_entry in length(outvec):1){
    if (ns[n_entry] > 1){
      replacement = list((outvec[n_entry] - (ns[n_entry] -1)) : outvec[n_entry])
      outvec = unlist(replace(outvec, n_entry, replacement))
    } else if (ns[n_entry] < -1){
      replacement = list( outvec[n_entry] : (outvec[n_entry] - (ns[n_entry]+ 1) ) )
      outvec = unlist(replace(outvec, n_entry, replacement))
    }
  }
    
  print(m2)
  
  return(outvec)

  
}



# # 
# 
# m1 = matrix(c(0,0,0,0,1,
#               0,-3,0,1,0,
#               0,0,1,0,0,
#               0,1,0,-1,0,
#               1,0,0,0,0), nrow=5)
# m1 = m1 * 10000
# colnames(m1) = c(10000,20000,50000,20000,10000)
# row.names(m1) = c(10000,20000,50000,20000,10000)
# m1 <- m1[seq(dim(m1)[1],1),]
# gm_to_segvec(m1)
