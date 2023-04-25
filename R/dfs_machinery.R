#' DFS caller function
#' 
#' This function serves as the gateway to entering the DFS process. It calls the dfsutil function recursively
#' to explore possible paths within the search space.
#'
#' @param bitlocus A matrix representing the search space.
#' @param maxdepth The maximum depth of the search tree.
#' @param increase_only A boolean indicating whether only increasing mutations are allowed.
#' @param earlystop The maximum number of evaluations before the search is terminated.
#' @param history An optional argument representing previous search history.
#' 
#' @return A list containing the search history and the evaluation results.
#' 
#' @export
#' @author Wolfram Höps
dfs <- function(bitlocus, maxdepth = 3, increase_only=F, earlystop = Inf, history = NULL){
  

  
  # Prepare initial bitlocus bit
  bitl = flip_bitl_y_if_needed(bitlocus)
  bitl = matrix_remove_zero_pads(bitl)
  # Get the orig_symmm value
  orig_symm = calc_symm(bitl)
  # Initialize a hash. 
  visited = hash::hash()
  ref_mut_pair = data.frame(p1='1', p2='1', sv='ref')
  ref_mut_pair = annotate_pairs_with_hash(bitl, ref_mut_pair)
  bitl_ref_mut = carry_out_compressed_sv(bitl, ref_mut_pair)
  
  # This is if it's called from a history mode
  if (!is.null(history)){
    print('Continuing a previous run.')
    pairhistory = history$mut_path
    depth = as.numeric(history$depth)
  } else {
    pairhistory = paste(ref_mut_pair[,1:3], collapse='_')
    depth = 0
  }
  
  # And enter the dfs rabbithole! (you'll stay there for a while,
  # if calls itself recursively.)
  visited_outputdf_list = dfsutil(visited=visited, 
                                  pair=ref_mut_pair,
                                  pairhash = ref_mut_pair$hash,
                                  mutator=bitl_ref_mut, 
                                  depth=depth, 
                                  maxdepth=maxdepth,
                                  pairhistory=pairhistory, 
                                  increase_only=increase_only,
                                  orig_symm = orig_symm,
                                  earlystop = earlystop)
  
  # Make sure output eval is numeric (perhaps not needed but unsure)
  visited_outputdf_list[[2]] = na.omit(as.data.frame(visited_outputdf_list[[2]]))
  visited_outputdf_list[[2]]$eval = as.numeric(visited_outputdf_list[[2]]$eval)
  
  return(visited_outputdf_list)

}




#' DFS workhorse function
#' 
#' This function recursively explores possible paths within the search space and returns a list containing
#' the search history and the evaluation results.
#'
#' @param visited A hash table containing the visited nodes.
#' @param pair A data frame representing a pair of points to be mutated.
#' @param pairhash A hash value representing the mutated pair.
#' @param mutator A matrix representing the mutated search space.
#' @param depth The current depth of the search tree.
#' @param maxdepth The maximum depth of the search tree.
#' @param pairhistory The search history up to the current pair.
#' @param df_output A data frame representing the evaluation results.
#' @param increase_only A boolean indicating whether only increasing mutations are allowed.
#' @param orig_symm The original symmetry value of the search space.
#' @param last_eval The score of the previous evaluation.
#' @param earlystop The maximum number of evaluations before the search is terminated.
#' 
#' @return A list containing the search history and the evaluation results.
#' 
#' @export
dfsutil <- function(visited, pair, pairhash, mutator, depth, maxdepth = 3, pairhistory=NULL, df_output=NULL, increase_only=NULL, orig_symm=1, last_eval=0, earlystop=Inf){
  

  # Calc score of a node. Force the calculation if we have ref (there we definitely want to know the value)
  aln_score = calc_coarse_grained_aln_score(mutator, forcecalc = (pair$sv == 'ref'), orig_symm = orig_symm)
  # if ((depth==2)){#} & (pair$p1 == 25) & (pair$p2 == 41)){#'25_41_inv'){
  #   browser()
  # }
  # Initiate an output_df if there is none
  if (is.null(df_output)){
    df_output = replicate(4, 
                          rep(NA, 1000), 
                          simplify=T)
    colnames(df_output) = c('hash','depth','mut_path','eval')
    len_df_output <<- dim(df_output)[1]
    n_res_counter <<- 1
    
  }
    

  # Make an entry to the output IF there is no NA and if there is
  # a calculated output
  if (!any(is.na(c(pairhash, depth, pairhistory, aln_score)))){
    if ((aln_score > 90) | (all(is.na(df_output))) |  (pairhistory == '1_1_ref')){
      df_output[n_res_counter,] = c(pairhash, depth, pairhistory, aln_score)
      n_res_counter <<- n_res_counter + 1
    }
  }
  # Update the visited hash to reflect that we have seen and processed that node. 
  visited[[pairhash]] = c(depth, pairhistory, aln_score)
  
  
  
  ########################### THIS PART NEEDS CONCEPTUAL DEVELOPMENT ########################
  # Criteria to continue. 
  #   2) If our node is perfect, no need to continue on that branch.
  #   3) If our node returns NA, we are somewhere weird
  #   4) on increasing depth, we have a cutoff crit for aln score. 
  

  # Make criterion 4:
  score_ref = as.numeric(df_output[1,'eval'])
  if (increase_only){
    min_score = rep(
        max(score_ref, max( as.numeric(df_output[which(df_output[,'depth'] == '0'),'eval'] ))),10
        )
  } else {
    min_score = c(0,0,0,0) # The second iteration has to come back to 75% of score_ref, otherwise we discontinue.
        # The rest here is usually not used...
        # (100-((100-score_ref)/0.6)), 
        # 0,0,0,(100-((100-score_ref)/0.8)),
        # (100-(100-score_ref)))
        # min_score = c(0, score_ref * 0.6, score_ref * 0.8, score_ref)
  }
  continue = ((
    #(aln_score < 100) &
      !is.na(aln_score) &
      (aln_score >= min_score[depth+1]) &
      (depth < maxdepth)))

  # continue = ((
  #     !is.na(aln_score) & 
  #     (depth < maxdepth)))
  ########################### END OF THIS PART  #######################################

  # We calculate the next batch of children
  if (!continue){
    n_discontinued <<- n_discontinued + 1
    return(list(visited, df_output))
  }
  
  pairs = find_sv_opportunities(mutator)
  
  if (dim(pairs)[1] > 0){
    # Give the pairs an initial hash
    pairs$hash = NA
  }
  #pairs = annotate_pairs_with_hash(mutator, pairs)
  
  for (npair in seq_along(row.names(pairs))){
    
    # Extend df_output when needed
    if (n_res_counter >= (len_df_output * 0.9)){
      print('Extending output matrix')
      df_output_addendum = replicate(4, 
                                     rep(NA, len_df_output), 
                                     simplify=T)
      df_output = rbind(df_output, 
                        df_output_addendum)
      len_df_output <<- dim(df_output)[1]
      
    } 
    
    # Be verbose 
    if (depth == 0){
      print(paste0('Processing branch ', npair, ' of ', dim(pairs)[1]))
    }
    
    if ((depth == 0) & (npair == earlystop)){
      return(list(visited, df_output))
    }
    
    # Carry out mutation
    bitl_mut = carry_out_compressed_sv(mutator, pairs[npair,1:3])
    
    # Add the hash of that bitl_mut to the visited
    hash = rlang::hash(return_diag_values_new(bitl_mut, fraction = 0.05))
    
    # Here comes the hash filter
    # Leave branch if we know the node already (?)
    node_is_novel_bool = is.null(visited[[hash]])
    if (!node_is_novel_bool){
      n_hash_excluded <<- n_hash_excluded + 1
      
      # The value of this hash has been computed
    if ((as.numeric(visited[[hash]][3])) > 1){
    vector_to_add_to_df_out = c(hash, 
                      as.numeric(depth+1), 
                      paste(pairhistory, paste0(pairs[npair,1:3], collapse = '_'), sep='+'),
                      as.numeric(visited[[hash]][3]))

    
      if  (!any(is.na(vector_to_add_to_df_out))){
        df_output[n_res_counter,] = vector_to_add_to_df_out
        n_res_counter <<- n_res_counter + 1
      }
    }
      next()
    }
    

    # Leave if the mutation is too unsymmetrical (unlikely path)
    node_passes_symmetry_crit = T#decide_loop_continue_symmetry(bitl_mut, orig_symm = orig_symm)
    if (!node_passes_symmetry_crit){
      visited[[hash]] = 0 #huh?
      next()
    }
  
    list_visit = dfsutil(visited=visited, 
                         pair=pairs[npair,], 
                         pairhash = hash,
                         mutator=bitl_mut, 
                         depth=depth+1, 
                         maxdepth=maxdepth,
                         pairhistory = paste(pairhistory, paste0(pairs[npair,1:3], collapse = '_'), sep='+'),
                         df_output = df_output,
                         increase_only = increase_only,
                         last_eval = aln_score,
                         orig_symm = orig_symm)
    
    visited = list_visit[[1]]
    df_output = list_visit[[2]]
  

    
  }
  

  return(list(visited, df_output))
}


