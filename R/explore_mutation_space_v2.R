


# 
# library(devtools)
# devtools::load_all()
# 
# # # bitlocusfull = readRDS('~/Desktop/gridmatrix_sample')
# # #bitlocusfull = readRDS('~/Desktop/blub.tsv')
#bitlocusfull = readRDS('~/Desktop/n2')
# # 
#bitlocus = bitlocusfull#[15:35, 15:35]

#' Turn a hash into a data.frame
#'
#' Transform a hash into a nrowx2 data.frame
#'  
#' @param x hash
#' @param key string; name for the column of keys (Default: "key").
#' @param value string; name for the column of values (Default: "value")
#' @param ... additional arguments to values
#'  
#' @return 
#' a data.frame with as many row as elements in `x` and two columns one for the
#' keys (always character) and value a list column
#'      
#'                
#' @examples 
#' 
#'  h <- hash( a=1, b=1:2, c=1:3)
#'  as.data.frame(h)
#'  
#' @export
as.data.frame.hash <- function(x, ..., key="key", value="value" ){
  df <- as.data.frame( list( keys(x) ), col.names=key )
  df[[value]] <- values(x)
  df 
}

#' annotate_pairs_with_hash
#' TODO: description
#' @export
annotate_pairs_with_hash <- function(bitlocus, pairs){
  
  # border cases
  # May 11th, Introducing this again after finding error $4.
  # I'm surprised why this was once here but apparently commented out. 
  if (dim(pairs)[1] == 0){
    return(pairs)
  }

  pairs$hash = 'NA'
  for (npair in seq_along(pairs$hash)){
    pairs[npair, 'hash'] = rlang::hash(as.numeric(carry_out_compressed_sv(bitlocus, pairs[npair,1:3])))
  }
  
  return(pairs)
}


decide_loop_continue <- function(bitl_f, symm_cutoff = 0.80, orig_symm = 1){

  climb_up_cost = as.numeric(row.names(bitl_f))
  walk_right_cost = as.numeric(colnames(bitl_f))
  symmetry = min(sum(climb_up_cost), sum(walk_right_cost)) / max(sum(climb_up_cost), sum(walk_right_cost))

  # Run away if there are at least 5 columns, and we have less than 75% symmetry
  if ((symmetry < (orig_symm * symm_cutoff)) & (dim(bitl_f)[1] > 5)){
    return(F)
  } 
  
  return(T)
}


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
      if (node_is_novel_bool){
        bitl_mut = carry_out_compressed_sv(mutator, pairs[npair,1:3])

        node_passes_symmetry_crit = T#decide_loop_continue(bitl_mut, orig_symm = orig_symm)
        
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
      } else if (visited[[pairs[npair, 'hash']]][1] > (depth + 1)){
        #browser()
        visited[[pairs[npair, 'hash']]][c(1,2)] = c(
          depth, 
          paste(pairhistory, paste0(pairs[npair,1:3], collapse = '_'), sep='+'))
        df_output[df_output$hash == pairs[npair, 'hash'],c(2,3,4)] =  c(depth+1, 
                                                                          paste(pairhistory, paste0(pairs[npair,1:3], collapse = '_'), sep='+'), 
                                                                          visited[[pairs[npair, 'hash']]][3])
        
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
  
  if ((n_pairs > 500) & (increase_only==F)){
    print('Huge Alignment! Going for depth 1. ')
    maxdepth = 1
  }
  
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


#' DFS caller function
#' TODO: description
#' @export
solve_mutation_slimmer <- function(bitlocus, depth){

  # Run a shallow run. 
  vis_list = (dfs(bitlocus, maxdepth = 1, increase_only = T))
  res_df = vis_list[[2]]
  res_df_sort = sort_new_by_penalised(bitlocus, res_df)

  # Did the shallow run give a result?
  if (res_df_sort1[1,'eval'] > find_threshold(bitlocus, res_df_sort1[1,'depth'])){
    print('Good news! Easy solution found :) No need to continue calculating.')
    res_out = transform_res_new_to_old(res_df_sort)
  } else { # If not, then we have to do the full calculus.
    print('No easy solution found. Continuing with complicated ones.')
    vis_list = (dfs(bitlocus, maxdepth = depth, increase_only = F))
    res_df = vis_list[[2]]
    res_df_sort = sort_new_by_penalised(bitlocus, res_df)
    res_out = transform_res_new_to_old(res_df_sort)
    print('Got so far.')
  }

  return(res_out)
}



#' DFS caller function
#' TODO: description
#' @export
solve_mutation <- function(bitlocus, depth){
  # reject weird stuff
  if (is_cluttered(bitlocus)){
    print("ÜÄH disgusting locus, make frame larger.")
    res_out = data.frame(eval=0, mut1='ref')

    return(res_out)
  }
  
  conclusion_found = F
  attempt = 1
  while (!conclusion_found){
    if (attempt == 1){
      # First, try if there is an easy solution on depth 1.
      vis_list = (dfs(bitlocus, maxdepth = 1, increase_only = T))
      res_df = vis_list[[2]]
      res_df$eval = as.numeric(res_df$eval)
      conclusion_found = (max(as.numeric(res_df$eval)) == 100)
      print(conclusion_found)
      if (conclusion_found){
        print('Good news! Very easy solution found :) No need to continue calculating.')
        res_df_sort = sort_new_by_penalised(bitlocus, res_df)
        res_out = transform_res_new_to_old(res_df_sort)
        res_out$eval = as.numeric(res_out$eval)
      }
    } else if (attempt == 2){
      print('No very easy solution found. Continuing steepest-descent solution search.')
      # Next, try if there is a combination of individual adaptations.
        n = 1
        #vis_list = dfs(bitlocus, maxdepth = 1, increase_only = T)
        res_df = vis_list[[2]]
        res_df_sort = sort_new_by_penalised(bitlocus, res_df)
        res_preferred = transform_res_new_to_old(res_df_sort)[1,]
        
        if (res_preferred[2] != 'ref'){
          bitlocus_new = carry_out_compressed_sv(bitlocus, data.frame(p1=strsplit(res_preferred[1,2],'_')[[1]][1],
                                                                      p2=strsplit(res_preferred[1,2],'_')[[1]][2],
                                                                      sv=strsplit(res_preferred[1,2],'_')[[1]][3]))
        } else {
          bitlocus_new = bitlocus
        }
        res_memory = res_preferred
        res_ref = transform_res_new_to_old(res_df_sort[res_df_sort$mut_path=='1_1_ref',])
        
        #browser()
        while ((res_preferred[2] != 'ref') & (n < 10)){
          n = n+1
          vis_list = (dfs(bitlocus_new, maxdepth = 1, increase_only = T))
          res_df = vis_list[[2]]
          res_df_sort = sort_new_by_penalised(bitlocus_new, res_df)
          res_preferred = transform_res_new_to_old(res_df_sort)[1,]
          if (res_preferred[2] != 'ref'){
            bitlocus_new = carry_out_compressed_sv(bitlocus_new, data.frame(p1=strsplit(res_preferred[1,2],'_')[[1]][1],
                                                                        p2=strsplit(res_preferred[1,2],'_')[[1]][2],
                                                                        sv=strsplit(res_preferred[1,2],'_')[[1]][3]))
            res_memory[['eval']] = res_preferred$eval
            res_memory[[paste0('mut',n)]] = res_preferred$mut1
          }
        }
        
        res_memory$eval = as.numeric(res_memory$eval)
        conclusion_found = res_memory$eval >= find_threshold(bitlocus_new, dim(res_memory)[2]-1)
        res_out = dplyr::bind_rows(res_memory, res_ref)
        
        
        print(conclusion_found)
    } else if (attempt == 3){
      print('No easy solution found. Continuing with complicated ones.')
      vis_list = (dfs(bitlocus, maxdepth = depth))
      res_df = vis_list[[2]]
      res_df_sort = sort_new_by_penalised(bitlocus, res_df)
      res_out = transform_res_new_to_old(res_df_sort)
      conclusion_found = T
    }
    attempt = attempt + 1  
  }
    
  # # Decide if a full run is necessary
  # if (max(as.numeric(res_df$eval)) == 100){
  #   print('Found an easy solution on depth=1. Not continuing the search.')
  # } else {
  # 
  # }
  # Make an entry to the output logfile #
  if (exists('log_collection')){
    log_collection$depth <<- depth
    log_collection$mut_simulated <<- dim(res_df)[1]
    log_collection$mut_tested <<- dim(res_df)[1]
  }

  res_out$eval = as.numeric(res_out$eval)
  return(res_out)
}

#' Transformation block. Ugly and undocumented i don't care sue me.
#' @export
transform_res_new_to_old <- function(res_df_f){

  #print(res_df_f) 
  if ((dim(res_df_f)[1] == 1) & (res_df_f$mut_path[1] == '1_1_ref')){
    res_out = data.frame(eval = res_df_f$eval, 
                         mut1 = 'ref')
    return(res_out)
  }
  max_mut = max(lengths(strsplit(res_df_f$mut_path, '\\+')))
  lits = (strsplit(res_df_f$mut_path, '\\+'))
  l <- lapply(lits, function(v) { c(v, rep(NA, max_mut-length(v)))})
  res_out = as.numeric(res_df_f$eval)
  for (i in 2:max_mut){
    res_out = cbind(res_out, sapply(l, "[[", i))
  }
  res_out = as.data.frame(res_out)
  colnames(res_out) = c('eval', paste0('mut', 1:(max_mut-1)))
  res_out[is.na(res_out$mut1),'mut1'] = 'ref'
  return(res_out)
}

#' Blub
#' @export
sort_new_by_penalised <- function(bitlocus_f, res_df_f){
  
  res_df_f$eval = as.numeric(res_df_f$eval)
  
  min_dotsize = min(as.numeric(row.names(bitlocus_f)), as.numeric(colnames(bitlocus_f)))
  step_penalty = (min_dotsize / sum(as.numeric(colnames(bitlocus_f))) * 100) * 1.5
  res_df_f$eval_penalised = as.numeric(res_df_f$eval) - (as.numeric(res_df_f$depth) * step_penalty)
  
  res_df_f = res_df_f[order(as.numeric(res_df_f$eval_penalised), decreasing=T ),]
  
  return(res_df_f)
}

#' blub2
#' @export
find_threshold <- function(bitlocus_f, nmut){
  min_dotsize = min(as.numeric(row.names(bitlocus_f)), as.numeric(colnames(bitlocus_f)))
  step_penalty = (min_dotsize / sum(as.numeric(colnames(bitlocus_f))) * 100) 
  
  threshold = 100 - (step_penalty * nmut)
  
  return(threshold)
}

#' to be documented
#' @export
is_cluttered <- function(bitlocus_f){
  is_large = (dim(bitlocus_f)[1] > 5 & dim(bitlocus_f)[2] > 5)
  
  if (!is_large){
    return(F)
  }
  clutter1 = sum(bitlocus_f[1                 ,2:(dim(bitlocus_f)[2]-1)] != 0)
  clutter2 = sum(bitlocus_f[dim(bitlocus_f)[1],2:(dim(bitlocus_f)[2]-1)] != 0)
  clutter3 = sum(bitlocus_f[2:(dim(bitlocus_f)[1]-1),1] != 0)
  clutter4 = sum(bitlocus_f[2:(dim(bitlocus_f)[1]-1),dim(bitlocus_f)[2]] != 0)
  clutter_all = c(clutter1, clutter2, clutter3, clutter4)
  
  if (any(clutter_all >= 5 )){
    return(T)
  }
  
  return(F)
}
# 
# bitlocusfull = readRDS('~/Desktop/latest')
# bitlocus = bitlocusfull#[15:35, 15:35]

