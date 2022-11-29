#' Find SV opportunities
#'
#' @description Building block: Take in a 'bitlocus',
#' and return a list of NAHR-driven rearrangements that
#' are possible on the reference sequence.
#'
#' @param sample  an nxm matrix (n=gridpoints_x, m=gridpoints_y)
#' @return A dataframe with columns 'p1', 'p2', 'sv'
#'
#' @author Wolfram Höps
#' @export
find_sv_opportunities <- function(sample) {
  double_SD_rows = which(rowSums(sample != 0) >= 2)
  all_opportunities = replicate(3, 
                                rep(NA, ((length(double_SD_rows)[1] ** 2) + 2)), 
                                simplify=T)
  colnames(all_opportunities) = NULL
  all_opp_len = dim(all_opportunities)[1]
  n_opportunity = 1
  
  # Rather disgusting block of if, else, ....
  # Probably can be implemented way more efficient.
  
  # For every row of the dot matrix, check first if this row
  # has more than one non-zero value in the matrix.
  for (nrow in 1:dim(sample)[1]) {
    if (!(nrow %in% double_SD_rows)) {
      next()
    }
    
    # ... If if does, now find the colnames of all those.
    identical_columns = as.numeric(which(sample[nrow, ] != 0))
    
    # make a matrix containing every combination of identical columns.
    # Imagine we have one segment repeated in columns 4, 5 and 9.
    # Then we get NAHR possibilities in 4-5, 4-9 and 5-9.
    # If this is only one pair, then the one pair stays the same.
    combs = combn(identical_columns, 2)
    
    for (ncomb in (1:dim(combs)[2])) {
      
      if (n_opportunity == all_opp_len){
        
        print('[Debug]: Extending size of all-opp matrix.')
        all_opportunities_addendum = replicate(3, 
                                      rep(NA, dim(all_opportunities)[1]), 
                                      simplify=T)
                                      
        all_opportunities = rbind(all_opportunities, 
                                  all_opportunities_addendum)     
        
        all_opp_len = dim(all_opportunities)[1]

      }
      
      identical_column_pair = combs[, ncomb]
      row_of_interest = sample[nrow,]
      
      # Determine relative orientation
      orientation = sign(row_of_interest[identical_column_pair[1]] * row_of_interest[identical_column_pair[2]])
      
      result_pair = c(identical_column_pair, orientation)
      all_opportunities[n_opportunity, ] = result_pair
      
      n_opportunity = n_opportunity + 1
      
    }
  }
  

  #all_opportunities = all_opportunities[rowSums(is.na(all_opportunities)) != ncol(all_opportunities),]
  all_opportunities = all_opportunities[1:(n_opportunity-1),, drop=F] # Replaces the previous line; is faster.
  colnames(all_opportunities) = c('p1', 'p2', 'sv')
  all_opportunities = na.omit(as.data.frame(unique(all_opportunities)))
  
  
  ## INV or DEL/DUP pair? ##
  # If at least one of the pair is negative, it's an 'inv'.
  if (dim(all_opportunities[all_opportunities$sv == -1, ])[1] > 0) {
    all_opportunities[all_opportunities$sv == -1 , ]$sv = 'inv'
  }
  # If at least one of the pair is positive, it's a 'del+dup' (???)
  if (dim(all_opportunities[all_opportunities$sv == 1, ])[1] > 0) {
    all_opportunities[all_opportunities$sv == 1 , ]$sv = 'del'
    # Duplicate plus pairs for dup and del.
    pluspairs = all_opportunities[all_opportunities$sv == 'del', ]
    pluspairs$sv = 'dup'
    all_opportunities = rbind(all_opportunities, pluspairs)
  }
  
  # No length-1-inversions
  all_opportunities = all_opportunities[!(
      ((all_opportunities$p2 - all_opportunities$p1) == 1) & 
      (all_opportunities$sv == 'inv')
    ),]
  
  return(all_opportunities)
  
}
