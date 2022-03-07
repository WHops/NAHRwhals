#' flip_bitl_y_if_needed
#' 
#' Early/quick and dirty version. Might need to be revised in the future. 
#'
#' @author  Wolfram Hoeps
#' @export
flip_bitl_y_if_needed <- function(bitl){
  
  cornerstones = c(0,0,0,0)
  i = 1
  while(all(cornerstones == 0)){
    cornerstones = c(
      bitl[i,i],
      bitl[nrow(bitl),ncol(bitl)],
      bitl[nrow(bitl),i],
      bitl[i,ncol(bitl)]
    )
    i = i+1
  }
  
  pos_corners = cornerstones[1:2]
  neg_corners = cornerstones[3:4]
  
  if (sum(pos_corners) == 0 & sum(neg_corners) <= 0){
    print('Inverse y axis detected. Flipping ... ')
    bitl = -bitl[nrow(bitl):1,]
    #flip
  }
  
  return(bitl)
}