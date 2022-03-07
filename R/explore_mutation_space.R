#' explore_mutation_space
#'
#' @description Main workhorse, tying together the pieces.
#' Hopefully we can replace this whole function with something nicer
#' one day. It is very s***y :P
#' @param bitlocus matrix, nxm
#' @param depth How many consecutive SVs should be simulated?
#' @return evaluation matrix
#'
#' @author Wolfram Höps
#' @export
explore_mutation_space <- function(bitlocus, depth) {
  

  # Consider to flip y axis. 
  
  bitlocus = flip_bitl_y_if_needed(bitlocus)

  pairs = find_sv_opportunities(bitlocus)
  
  if (dim(pairs)[1] == 0){
    res = (matrix(ncol = depth + 1, nrow = 1))
    colnames(res) = c('eval', paste0('mut', 1:depth))
    res[1, ] = unlist(c(calc_coarse_grained_aln_score(bitlocus, forcecalc=T), 'ref', rep('NA', depth - 1)))
    return(as.data.frame(res))
  }
  

  # Del-dup-direction: do we need to delete or duplicate overall? 
  del_dup_direction = ((dim(bitlocus)[1] / dim(bitlocus)[2]) - 1)
  if (abs(del_dup_direction) < 0.1){ del_dup_direction == 0}
  del_dup_direction = sign(del_dup_direction)
  
  pairs = remove_equivalent_pairs(pairs)
  
  res = (matrix(ncol = depth + 1, nrow = (dim(pairs)[1] ** depth) * 5))
  colnames(res) = c('eval', paste0('mut', 1:depth))
  res[1, ] = unlist(c(calc_coarse_grained_aln_score(bitlocus, forcecalc=T), 'ref', rep('NA', depth - 1)))

  if (as.numeric(res[1,'eval'] > 95)){
    print('Initial alignment is better than 95%. Not attempting to find SV.')
    res[2,] = res[1,]
    return(as.data.frame(res[1:2,]))
  }
  
  rescount = 2
  
  for (npair_level1 in 1:dim(pairs)[1]) {
    print(paste0('Entering front layer ', npair_level1, ' of ', dim(pairs)[1]))
    
    pair_level1 = pairs[npair_level1, ]
    
    # Trio of function: Mutate, Evaluate, Observate
    bitl_mut = carry_out_compressed_sv(bitlocus, pair_level1)
    
    res[rescount, 1:2] = add_eval(res,
                                  m = bitl_mut,
                                  layer = 1,
                                  pair_level1,
                                  NULL,
                                  NULL)

    rescount = rescount + 1
    
    if (depth > 1) {
      
      newpairs = find_sv_opportunities(bitl_mut)
      newpairs = remove_equivalent_pairs(newpairs)
      
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
            
            newpairs3 = remove_equivalent_pairs(newpairs3)
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
  
  #print('hi')
  res_df = as.data.frame(res)
  res_df$eval = as.numeric(res_df$eval)
  
  # Remove NA rows
  res_df = res_df[rowSums(is.na(res_df)) != ncol(res_df),]
  
  # Remove entries without evals
  #res_df = res_df[!is.na(res_df$eval),]
  
  
  print(paste0(
    'Finished ',
    dim(res_df)[1],
    ' mutation simulations (',
    dim(res_df[!is.na(res_df$eval), ])[1],
    ' with eval calculated)'
  ))
  
  return(res_df)
}


