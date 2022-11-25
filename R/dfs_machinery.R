#' DFS workhorse function
#' TODO: description
#' @export
dfsutil <- function(visited, pair, mutator, depth, maxdepth = 3, pairhistory=NULL, df_output=NULL, increase_only=NULL, orig_symm=1, last_eval=0){
  
  # Calc score of a node
  if (pair$sv == 'ref'){
    aln_score = calc_coarse_grained_aln_score(mutator, forcecalc = T, orig_symm = orig_symm)
  } else {
    aln_score = calc_coarse_grained_aln_score(mutator, forcecalc = F, orig_symm = orig_symm)
    
  }
  # Make an entry to the output df. 
  if (is.null(df_output)){
    df_output = data.frame(matrix(nrow = 1, ncol=4)) 
    colnames(df_output) = c('hash','depth','mut_path','eval')
  }
  
  
  df_output = na.omit(rbind(df_output, c(pair[1,'hash'], depth, pairhistory, aln_score)))
  #df_output = (rbind(df_output, c(pair[1,'hash'], depth, pairhistory, aln_score)))
  
  
  # Update the visited hash to reflect that we have seen and processed that node. 
  visited[[pair$hash]] = c(depth, pairhistory, aln_score)
  
  # Criteria to continue. 
  #   2) If our node is perfect, no need to continue on that branch.
  #   3) If our node returns NA, we are somewhere weird
  #   4) on increasing depth, we have a cutoff crit for aln score. 
  
  # Make criterion 4:
  score_ref = as.numeric(df_output[df_output$mut_path == '1_1_ref',]$eval)
  if (increase_only){
    min_score = rep(c(max(score_ref, max(as.numeric(df_output[df_output$depth==depth,]$eval)))),10)
  } else {
    min_score = c(0, 
                  0, 
                  (100-((100-score_ref)/0.6)), 
                  (100-((100-score_ref)/0.8)),
                  (100-(100-score_ref)))
    #min_score = c(0, score_ref * 0.6, score_ref * 0.8, score_ref)
  }
  continue = ((
    (aln_score < 100) & 
      !is.na(aln_score) & 
      (aln_score >= min_score[depth+1]) &
      (depth < maxdepth)))
  #print(continue)
  # We caclualte the next batch of children
  if (continue){
    pairs = find_sv_opportunities(mutator)
    pairs = annotate_pairs_with_hash(mutator, pairs)
    
    
    for (npair in seq_along(pairs$hash)){
      if (depth == 0){
        print(paste0('Processing branch ', npair, ' of ', dim(pairs)[1]))
      }
      # if(paste(pairhistory, paste0(pairs[npair,1:3], collapse = '_'), sep='+') == '16_20_dup+16_20_del'){
      #   browser()
      # }
      node_is_novel_bool = is.null(visited[[pairs[npair, 'hash']]])
      # if (!node_is_novel_bool){
      #   browser()
      #   print('skipping a visited node.')
      # }
      if (node_is_novel_bool){
        bitl_mut = carry_out_compressed_sv(mutator, pairs[npair,1:3])
        
        node_passes_symmetry_crit = T#decide_loop_continue_symmetry(bitl_mut, orig_symm = orig_symm)
        
        if (node_passes_symmetry_crit){
          list_visit = dfsutil(visited=visited, 
                               pair=pairs[npair,], 
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
        } else {
          visited[[pairs[npair,'hash']]] = 0
        }
        # If node is NOT novel, BUT has only been seen at higher depth,
        # update the path and depth. 
      # } else if (visited[[pairs[npair, 'hash']]][1] > (depth + 1)){
      #   #browser()
      #   visited[[pairs[npair, 'hash']]][c(1,2)] = c(
      #     depth, 
      #     paste(pairhistory, paste0(pairs[npair,1:3], collapse = '_'), sep='+'))
      #   df_output[df_output$hash == pairs[npair, 'hash'],c(2,3,4)] =  c(depth+1, 
      #                                                                   paste(pairhistory, paste0(pairs[npair,1:3], collapse = '_'), sep='+'), 
      #                                                                   visited[[pairs[npair, 'hash']]][3])
       
      } else {
        n_hash_excluded <<- n_hash_excluded + 1
        df_output = na.omit(rbind(df_output, c(pairs[npair,'hash'], 
                                               as.numeric(depth+1), 
                                               paste(pairhistory, paste0(pairs[npair,1:3], collapse = '_'), sep='+'),
                                               as.numeric(visited[[pairs[npair, 'hash']]][3]))))

      }
      
    }
  }
  return(list(visited, df_output))
}


#' DFS caller function
#' TODO: description
#' @export
dfs <- function(bitlocus, maxdepth = 3, increase_only=F){
  
  # Prepare bitlocus that we will be working on. 
  bitl = flip_bitl_y_if_needed(bitlocus)
  
  n_pairs = dim(find_sv_opportunities(bitl))[1]
  if ((n_pairs > 200) & (maxdepth == 3) & (increase_only==F)){
    print('Uh oh that is a bit large. Reducing depth to 2. Try to avoid producing such large alignments.')
    maxdepth = 2
  }
  
  # if ((n_pairs > 500) & (increase_only==F)){
  #   print('Huge Alignment! Going for depth 1. ')
  #   maxdepth = 1
  # }
  
  # What is the initial symmetry of the bitlocus? 
  climb_up_cost = as.numeric(row.names(bitl))
  walk_right_cost = as.numeric(colnames(bitl))
  orig_symm = min(sum(climb_up_cost), sum(walk_right_cost)) / max(sum(climb_up_cost), sum(walk_right_cost))
  
  # Initialize a hash. 
  visited = hash::hash()
  ref_mut_pair = data.frame(p1='1', p2='1', sv='ref')
  ref_mut_pair = annotate_pairs_with_hash(bitl, ref_mut_pair)
  bitl_ref_mut = carry_out_compressed_sv(bitl, ref_mut_pair)
  
  visited_outputdf_list = dfsutil(visited=visited, 
                                  pair=ref_mut_pair,
                                  mutator=bitl_ref_mut, 
                                  depth=0, 
                                  maxdepth=maxdepth,
                                  pairhistory=paste(ref_mut_pair[,1:3], collapse='_'), 
                                  increase_only=increase_only,
                                  orig_symm = orig_symm)
  
  
  visited_outputdf_list[[2]]$eval = as.numeric(visited_outputdf_list[[2]]$eval)
  
  return(visited_outputdf_list)
  
  
}