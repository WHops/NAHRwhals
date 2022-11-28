# Whoeps, 24th Nov 2022
# This file harbours second-level wrapper functions. The 
# solving 'workhorse' fonctions dfs and dfs_utils are found 
# in a seaprate script, dfs_utils.R




#' Wrapper for checking depth 1 and continue to full depth if this failed.
#' TODO: description
#' @export
solve_mutation_slimmer <- function(bitlocus, depth){

  # Run a shallow run. 
  vis_list = (dfs(bitlocus, maxdepth = 1, increase_only = T))
  res_df = vis_list[[2]]
  res_df_sort = sort_new_by_penalised(bitlocus, res_df)

  # Did the shallow run give a result? If yes, we can leave now. 
  if (res_df_sort1[1,'eval'] > find_threshold(bitlocus, res_df_sort1[1,'depth'])){
    print('Good news! Easy solution found :) No need to continue calculating.')
    res_out = transform_res_new_to_old(res_df_sort)
    
    return(res_out)
  }
  
  # Not sure if this makes a lot of sense. 
  print('No easy solution found. Continuing with complicated ones.')
  vis_list = (dfs(bitlocus, maxdepth = depth, increase_only = F))
  res_df = vis_list[[2]]
  res_df_sort = sort_new_by_penalised(bitlocus, res_df)
  res_out = transform_res_new_to_old(res_df_sort)

  return(res_out)
}


#' Alternative wrapper for calling dfs several times.  
#' THis function is called as a standard.
#' TODO: description
#' @export
solve_mutation <- function(bitlocus, maxdepth, earlystop = Inf){
  
  n_eval_total <<- 0
  n_eval_calc <<- 0
  n_hash_excluded <<- 0
  
  
  # reject bitloci with cut-off alignments at the borders ('clutter')
  if (is_cluttered(bitlocus)){
    print("Alignments overlap with borders. Make frame larger!")
    res_out = data.frame(eval=0, mut1='ref')
    
    return(res_out)
  }
  
  # Optional: adjust depth
  #maxdepth = reduce_depth_if_needed(bitlocus, increase_only = F, maxdepth)
  
  # Run iterations of increasing depth. 
  # Always run the whole thing, but also don't go 
  # deeper once a good solution is found. 
  conclusion_found = F
  current_depth = 1
  
  while ((!conclusion_found) & (current_depth <= maxdepth)){
    print(paste0('Running depth layer: ', current_depth))
    
    # Run the dfs machinery (main work horse)
    vis_list = (dfs(bitlocus, maxdepth = current_depth, increase_only = F, earlystop = earlystop))
    
    # Extract res_df
    res_df = vis_list[[2]]

    # Are we happy with the best result? 
    solve_th = find_threshold(bitlocus, current_depth)
    conclusion_found = (max(as.numeric(res_df$eval)) >= solve_th)
    
    if (conclusion_found){
      print("Conclusion found!!")
    }
    
    # Higher depth for next iteration!
    current_depth = current_depth + 1  
  }
  
  # Prepare output for returns
  res_df_sort = sort_new_by_penalised(bitlocus, res_df)
  res_out = transform_res_new_to_old(res_df_sort)
  res_out$eval = as.numeric(res_out$eval)
  
  
  # Make an entry to the output logfile #
  if (exists('log_collection')){
    log_collection$depth <<- maxdepth
    log_collection$mut_simulated <<- dim(res_df)[1]
    log_collection$mut_tested <<- dim(res_df)[1]
  }
  
  print(paste0('Eval total: ', n_eval_total))
  print(paste0('Eval calced: ', n_eval_calc))
  print(paste0('Hash excluded: ', n_hash_excluded))
  return(res_out)

}
  

















#' Alternative wrapper for calling dfs several times.  
#' TODO: description
#' @export
solve_mutation_old <- function(bitlocus, depth, discovery_exact){
  
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
      vis_list = (dfs(bitlocus, maxdepth = 1, increase_only = T)) # Increase only: we go uphill only. no local minima explored
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
      
      if (discovery_exact == T){
        print('No very easy solution found. Skipping to exact solution search')
        conclusion_found = F
      } else {
        
      
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
      }
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





# 
# bitlocusfull = readRDS('~/Desktop/latest')
# bitlocus = bitlocusfull#[15:35, 15:35]


# 

# 
bitlocusfull = readRDS('~/Desktop/geruempel/desktop_23rd_may_2022/transfer2trash/gridmatrix_sample')
# # #bitlocusfull = readRDS('~/Desktop/geruempel/desktop_23rd_may_2022/transfer2trash/blub.tsv')
#bitlocusfull = readRDS('~~/Desktop/geruempel/desktop_23rd_may_2022/transfer2trash/n2')
