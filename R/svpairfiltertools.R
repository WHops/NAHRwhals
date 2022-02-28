
#' find_nb_alt
#'
#' @description Cut down a list of pairs to representatives of a list. 
#' @param pairs dataframe with columns p1 and p2
#' @param mode inv or deldup.
#' @return pairs dataframe, but with fewer entries
#'
#' @author Wolfram Höps
#' @export
find_nb_alt <- function(pairs, mode){
  
  # Make sure mode is inv or deldup.
  stopifnot("mode must be inv or deldup." = mode %in% c('inv','deldup'))
  
  # Make sure there is at least one pair in the dataframe.
  if (dim(pairs)[1] == 0){
    return(NULL)
  }
  
  pairs_min1 = pairs
  
  if (mode == 'inv'){
    pairs_min1$p1 = pairs_min1$p1-1
    pairs_min1$p2 = pairs_min1$p2+1
  } else if (mode == 'deldup'){
    pairs_min1$p1 = pairs_min1$p1+1
    pairs_min1$p2 = pairs_min1$p2+1
  }
  remainder = suppressMessages(dplyr::anti_join(pairs, pairs_min1))
  
  return(remainder)
  
}



#' filter_length_sense
#'
#' @author Wolfram Höps
#' @export
filter_length_sense <-
  function(oldpairs, newpairs, del_dup_direction) {
    # -1: we need at least 1 deletion
    if (del_dup_direction == -1) {
      if (!('del' %in% oldpairs$sv)) {
        newpairs = newpairs[newpairs$sv == 'del', ]
      }
    }
    
    if (del_dup_direction == 1) {
      if (!('dup' %in% oldpairs$sv)) {
        newpairs = newpairs[newpairs$sv == 'dup', ]
      }
    }
    
    return(newpairs)
  }


#' remove_circular_chains
#'
#' @author Wolfram Höps
#' @export
remove_circular_chains <- function(oldpair, newpairs) {
  # If the last one was an inversion, don't invert the same back.
  # If the last one was a duplication, don't delete away the changed part.
  if (oldpair$sv == 'inv') {
    newpairs = newpairs[!((newpairs$p1 == oldpair$p1) &
                            (newpairs$p2 == oldpair$p2)),]
  } else if (oldpair$sv == 'dup') {
    # In order of execution:
    # deleting the part that was just duplicated
    # deleting the whole thing
    # deleting the duplicated part
    # deleting any newly duplicated part.
    newpairs = newpairs[!(((newpairs$p1 == oldpair$p1) &
                             (newpairs$p2 == oldpair$p2) & (newpairs$sv == 'del')
    ) |
      
      ((newpairs$p1 == (
        oldpair$p1 + (oldpair$p2 - oldpair$p1)
      )) &
        (newpairs$p2 == (
          oldpair$p2 + (oldpair$p2 - oldpair$p1)
        )) & (newpairs$sv == 'del')) |
      ((newpairs$p1 == oldpair$p1) &
         (newpairs$p2 == (oldpair$p2 + (
           oldpair$p2 - oldpair$p1
         ))) & (newpairs$sv == 'del')
      ) |
      ((newpairs$p1 > oldpair$p1) &
         (newpairs$p1 < oldpair$p2) &
         ((newpairs$p2 - newpairs$p1) == (oldpair$p2 - oldpair$p1))
      )),]
  }
  
  return(newpairs)
  
}




#' remove_equivalent_pairs
#'
#' @description Wrapper for the neighbour-removal.
#' @param pairs data frame, p1, p2
#' @return the same data frame, but fewer entries. 
#'
#' @author Wolfram Höps
#' @export
remove_equivalent_pairs <- function(pairs){
  
  # First deal with 'inversion' pairs.
  pairs_inv = pairs[pairs$sv == 'inv',]
  pairs_inv = pairs_inv[order(pairs_inv$p1),]
  svs_inv = find_nb_alt(pairs_inv, mode='inv')
  
  # Also remove inversion of palindrome.
  svs_inv = svs_inv[svs_inv$p1+1 != svs_inv$p2,]
  
  # Next, deal with del / dup pairs. 
  pairs_del = pairs[pairs$sv == 'del',]
  pairs_del = pairs_del[order(pairs_del$p1),]
  svs_del = find_nb_alt(pairs_del, mode='deldup')
  
  # Having filered dup pairs, we make a new 'del' pairlist.
  svs_dup = svs_del
  if (!is.null(svs_del)){
    svs_dup$sv = 'dup'
  }
  
  # Stitch everything back together
  pairs_all = rbind(svs_inv, rbind(svs_del, svs_dup))
  
  # Return. 
  if (is.null(pairs_all)){
    return(data.frame())
  } else {
  return(pairs_all)
  }
}
