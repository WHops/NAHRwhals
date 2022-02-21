#' ColMax function
#'
#' From Stackoverflow. Does what you expect it to do.
#' @export
colMax <- function(data)
  sapply(data, max, na.rm = TRUE)


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
  all_opportunities = do.call(data.frame, 
                              replicate(3, 
                                        rep(NA, ((length(double_SD_rows)[1] ** 2) / 2)), 
                                        simplify=FALSE))
  colnames(all_opportunities) = NULL
  
  n_opportunity = 1
  # Rather disgusting block of if, else, ....
  # Probably can be implemented way more efficient.
  
  # For every row of the dot matrix, check first if this row
  # has more than one non-zero value in the matrix.
  for (nrow in 1:dim(sample)[1]) {
    if (nrow %in% double_SD_rows) {
      # ... If if does, now find the colnames of all those.
      identical_columns = as.numeric(which(sample[nrow, ] != 0))
      
      # make a matrix containing every combination of identical columns.
      # Imagine we have one segment repeated in columns 4, 5 and 9.
      # Then we get NAHR possibilities in 4-5, 4-9 and 5-9.
      # If this is only one pair, then the one pair stays the same.
      combs = combn(identical_columns, 2)
      
      for (ncomb in (1:dim(combs)[2])) {
        identical_column_pair = combs[, ncomb]
        
        # Determine relative orientation
        orientation = sign(sample[nrow, identical_column_pair[1]] * sample[nrow, identical_column_pair[2]])
        
        result_pair = c(identical_column_pair, orientation)
        all_opportunities[n_opportunity, ] = result_pair
        
        n_opportunity = n_opportunity + 1
        
      }
    }
  }
  
  all_opportunities = all_opportunities[rowSums(is.na(all_opportunities)) != ncol(all_opportunities),]
  
  colnames(all_opportunities) = c('p1', 'p2', 'sv')
  all_opportunities = unique(all_opportunities)
  
  
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
  
  
  return(all_opportunities)
  
}





#' Find SV carry_out_compressed_sv
#'
#' @description TBD
#'
#' @param bitl  an nxm matrix (n=gridpoints_x, m=gridpoints_y)
#' @param input_ins a vector with instruction: p1, p2, sv.
#'
#' @author Wolfram Höps
#' @export
carry_out_compressed_sv <- function(bitl, input_ins) {
  pair = input_ins[1:2]
  action = input_ins[3]
  
  
  # Modify the matrix accordingly
  if (action == 'del') {
    bitl_mut = bitl[, -c(as.numeric(pair[1]):(as.numeric(pair[2]) - 1))]
  } else if (action == 'dup') {
    bitl_mut = cbind(bitl[, 1:as.numeric(pair[2])],
                     bitl[, as.numeric(pair[1] + 1):dim(bitl)[2]])
  } else if (action == 'inv') {
    bitl_mut = cbind(cbind(bitl[, 1:as.numeric(pair[1])],-bitl[, (as.numeric(pair[2]) -
                                                                    1):(as.numeric(pair[1]) + 1)]),
                     bitl[, as.numeric(pair[2]):dim(bitl)[2]])
  }
  
  colnames(bitl_mut) = 1:dim(bitl_mut)[2]
  
  return(bitl_mut)
  
}


#' filter_length_sense
#'
#' @author Wolfram Höps
#' @export
filter_length_sense <- function(oldpairs, newpairs, del_dup_direction){
  
  # -1: we need at least 1 deletion
  if (del_dup_direction == -1){
    if (! ('del' %in% oldpairs$sv)){
      newpairs = newpairs[newpairs$sv == 'del',]
    }
  }
  
  if (del_dup_direction == 1){
    if (! ('dup' %in% oldpairs$sv)){
      newpairs = newpairs[newpairs$sv == 'dup',]
    }
  }
  
  return(newpairs)
}


#' filter_newpairs_for_sense
#'
#' @author Wolfram Höps
#' @export
remove_circular_chains <- function(oldpair, newpairs){
  

  # If the last one was an inversion, don't invert the same back.
  # If the last one was a duplication, don't delete away the changed part. 
  if (oldpair$sv == 'inv'){
    newpairs = newpairs[!((newpairs$p1 == oldpair$p1) & (newpairs$p2 == oldpair$p2)), ]
  } else if (oldpair$sv == 'dup'){
    
    
    newpairs = newpairs[!( ((newpairs$p1 == oldpair$p1) & (newpairs$p2 == oldpair$p2) & (newpairs$sv == 'del')) |   # deleting the part that was just duplicated
                             
                           ((newpairs$p1 == (oldpair$p1 + (oldpair$p2 - oldpair$p1))) & (newpairs$p2 == (oldpair$p2 + (oldpair$p2 - oldpair$p1))) & (newpairs$sv == 'del')) | # deleting the whole thing
                             
                           ((newpairs$p1 == oldpair$p1) & (newpairs$p2 == (oldpair$p2 + (oldpair$p2 - oldpair$p1))) & (newpairs$sv == 'del')) | # deleting the duplicated part
                             
                          ((newpairs$p1 > oldpair$p1) & (newpairs$p1 < oldpair$p2) & ((newpairs$p2 - newpairs$p1) == (oldpair$p2 - oldpair$p1)))  # deleting any newly duplicated part. 
                             ), ]
  }
  
  return(newpairs)
  
}


#' explore_mutation_space
#'
#' @description Main workhorse, tying together the pieces.
#' @param bitlocus matrix, nxm
#' @param depth How many SVs in sequence should be simulated?
#' @return evaluation matrix
#'
#' @author Wolfram Höps
#' @export
explore_mutation_space <- function(bitlocus, depth) {
  
  sample = bitlocus
  pairs = find_sv_opportunities(sample)
  
  # Del-dup-direction: do we need to delete or duplicate overall? 
  del_dup_direction = ((dim(bitlocus)[1] / dim(bitlocus)[2]) - 1)
  if (abs(del_dup_direction) < 0.1){ del_dup_direction == 0}
  del_dup_direction = sign(del_dup_direction)
  
  
  res = (matrix(ncol = depth + 1, nrow = (dim(pairs)[1] ** depth) * 5))
  
  colnames(res) = c('eval', paste0('mut', 1:depth))
  res[1, ] = unlist(c(eval_mutated_seq(bitlocus), 'ref', rep('NA', depth - 1)))
  # Init res with the reference result.
  
  rescount = 2
  
  
  for (npair_level1 in 1:dim(pairs)[1]) {
    print(paste0('Entering front layer ', npair_level1, ' of ', dim(pairs)[1]))
    
    pair_level1 = pairs[npair_level1, ]
    
    # Trio of function: Mutate, Evaluate, Observate
    bitl_mut = carry_out_compressed_sv(sample, pair_level1)
    
    res[rescount, 1:2] = add_eval(res,
                                  m = bitl_mut,
                                  layer = 1,
                                  pair_level1,
                                  NULL,
                                  NULL)

    rescount = rescount + 1
    
    if (depth > 1) {
      
      newpairs = find_sv_opportunities(bitl_mut)
      newpairs = remove_circular_chains(pair_level1, newpairs)
      
      if (depth == 2){
        newpairs = filter_length_sense(pair_level1, newpairs, del_dup_direction)
      }
      
      if (dim(newpairs)[1] > 0) {
        for (npair_level2 in 1:dim(newpairs)[1]) {
          #print(paste0('Entering second layer ', npair_level2, ' of ', dim(newpairs)[1]))
          
          pair_level2 = newpairs[npair_level2, ]
          
          # Trio of function: Mutate, Evaluate, Observate
          bitl_mut2 = carry_out_compressed_sv(bitl_mut, pair_level2)
          res[rescount, 1:3] = add_eval(res,
                                        m = bitl_mut2,
                                        layer = 2,
                                        pair_level1,
                                        pair_level2,
                                        NULL)
          rescount = rescount + 1
          
          if (depth > 2) {
            newpairs3 = find_sv_opportunities(bitl_mut2)
            newpairs3 = remove_circular_chains(pair_level2, newpairs3)
            
            if (depth == 3){
              newpairs3 = filter_length_sense(rbind(pair_level1, pair_level2), newpairs3, del_dup_direction)
            }
            
            if (dim(newpairs3)[1] > 0) {
              for (npair_level3 in 1:dim(newpairs3)[1]) {
                pair_level3 = newpairs3[npair_level3, ]
                
                # Duo of function: Mutate, Evaluate
                bitl_mut3 = carry_out_compressed_sv(bitl_mut2, pair_level3)
                res[rescount, 1:4] = add_eval(
                  res,
                  m = bitl_mut3,
                  layer = 3,
                  pair_level1,
                  pair_level2,
                  pair_level3
                )
                
                rescount = rescount + 1
                
              } #loop 3 end
            } #if depth > 1 end
          }
        }
      } # loop 2 end
    } # if depth > 1 end
  } # loop 1 end
  
  print('hi')
  res_df = as.data.frame(res)
  res_df$eval = as.numeric(res_df$eval)
  
  # Remove NA rows
  res_df = res_df[rowSums(is.na(res_df)) != ncol(res_df),]
  
  print(paste0(
    'Finished ',
    dim(res_df)[1],
    ' mutation simulations (',
    dim(res_df[res_df$eval != 0, ])[1],
    ' with eval calculated)'
  ))
  
  return(res_df)
}



# Sample case
# #
# C = gimme_sample_matrix(mode = 'diff')
# a = explore_mutation_space_(C, depth = 3)
