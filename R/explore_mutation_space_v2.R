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
solve_mutation <- function(bitlocus, maxdepth_input, earlystop = Inf){
  # 
  # browser()
  # bitlocus = readRDS('bitme_large')
  # maxdepth_input = 3
  # earlystop = Inf
  #ref_mut_pair = data.frame(p1='12', p2='25', sv='inv')
  #bitlocus = carry_out_compressed_sv(bitlocus, ref_mut_pair)
  est_highest <<- 0
  n_eval_total <<- 0
  n_eval_calc <<- 0
  n_hash_excluded <<- 0
  n_discontinued <<- 0
  
  # reject bitloci with cut-off alignments at the borders ('clutter')
  if (is_cluttered(bitlocus)){
    print("Alignments overlap with borders. Make frame larger!")
    res_out = data.frame(eval=0, mut1='ref')
    res_out = annotate_res_out_with_positions_lens(res_out, bitlocus)
    
    return(res_out)
  }
  
  # Optional: adjust depth
  maxdepth = 3#reduce_depth_if_needed(bitlocus, increase_only = F, maxdepth_input)
  print(maxdepth)
  # Run iterations of increasing depth. 
  # Always run the whole thing, but also don't go 
  # deeper once a good solution is found. 
  conclusion_found = F
  current_depth = 1
  
  # The tree is explored in an iterative manner, in which all mutations from the one depth are
  # explored. If one depth does not yield a solution, the next depth is explored, until
  # the largest depth (3) is reached. 
  # In cases where >1M mutations are expected in depth 3, (corresponding to 100 mutations on depth 1),
  # the tree is explored until depth 2, and only the 100 highest-scoring depth-2 mutations are followed up
  # in layer three. This approach cuts computation time drastically, but may miss mutation chains behind 
  # local minima in rare instances. Additionally, it is not guaranteed that all alternative paths to an 
  # optimal result are identified with this. 
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
  
  # OPTIONAL: THIS IS THE STEEPEST DESCENT OF LAYER 3
  # IF max depth was reduced due to size, and if no conclusion was found, 
  # then we go into a steepest-descent bonus round. 
  if ((current_depth == 2) & (maxdepth < maxdepth_input) & (conclusion_found == F)){

    res_df = run_steepest_descent_one_layer_down(bitlocus, res_df, starting_depth = current_depth - 1)
    
    if (conclusion_found == F){
      res_df = run_steepest_descent_one_layer_down(bitlocus, res_df, starting_depth = current_depth)
    }
    
  }  
  
  if ((current_depth == 3) & (maxdepth < maxdepth_input) & (conclusion_found == F)){

    res_df = run_steepest_descent_one_layer_down(bitlocus, res_df, starting_depth = current_depth - 1)
  } 
  

  print('Sorting results')

  # Prepare output for returns
  res_df_sort = sort_new_by_penalised(bitlocus, res_df)
  res_out = transform_res_new_to_old(res_df_sort)
  res_out$eval = as.numeric(res_out$eval)
  
  res_out = annotate_res_out_with_positions_lens(res_out, bitlocus)
  

  # Make an entry to the output logfile #
  if (exists('log_collection')){
    log_collection$depth <<- maxdepth
    log_collection$mut_simulated <<- dim(res_df)[1]
    log_collection$mut_tested <<- dim(res_df)[1]
  }

  print(paste0('Nodes considered: ', n_eval_total + n_hash_excluded))  
  print(paste0('Eval attempted: ', n_eval_total))
  print(paste0('Eval calced: ', n_eval_calc))
  print(paste0('Hash excluded: ', n_hash_excluded))

  return(res_out)

}
  


#' Description...
#' @export
annotate_res_out_with_positions_lens <- function(res_out, bitlocus){
  

  n_mut = dim(res_out)[2] - 1
  res_out$mut1_start = NA
  res_out$mut1_end = NA
  res_out$mut1_pos_pm = NA
  res_out$mut1_len = NA
  res_out$mut1_len_pm = NA
  
  res_out$mut2_len = NA
  res_out$mut2_len_pm = NA
  res_out$mut3_len = NA
  res_out$mut3_len_pm = NA
  
  if (is.null(bitlocus)){
    return(res_out)
  }
  
  colnumbers = as.numeric(colnames(bitlocus))
  
  count = 1
  # Determine location of mut #1
  for (nrow in row.names(res_out)){
    count = count + 1
    
    if (count > 2000){
      print('Stopping to annotate results to save some time.')
      return(res_out)
    }
    res_vec = res_out[nrow, ]
    nmut = 
      if (res_vec$mut1 == 'ref'){
        next()
      }
    
    # mut 1
    m1_places = as.numeric(strsplit(res_vec$mut1, '_')[[1]][1:2]) 
    start_mean = round(sum(colnumbers[1:(m1_places[1] -1)]) + (colnumbers[m1_places[1]] / 2))
    end_mean = round(sum(colnumbers[1:(m1_places[2] -1)]) + (colnumbers[m1_places[2]] / 2))
    
    # PM: plus/minus. 
    start_pm =  round((colnumbers[m1_places[1]] / 2))
    end_pm =  round((colnumbers[m1_places[2]] / 2))
    
    len_range = c(
      (end_mean-end_pm)-(start_mean+start_pm), # Shortest possible
      (end_mean+end_pm)-(start_mean-start_pm) # Largest possible
    )
    
    res_out[nrow, 'mut1_start'] = start_mean
    res_out[nrow, 'mut1_end'] = end_mean
    res_out[nrow, 'mut1_pos_pm'] = start_pm
    res_out[nrow, 'mut1_len'] = mean(len_range)
    res_out[nrow, 'mut1_len_pm'] = len_range[2] - mean(len_range)
    
    if ('mut2' %in% colnames(res_vec[,colSums(is.na(res_vec))<nrow(res_vec)])){

      bitlocus_mut = carry_out_compressed_sv(bitlocus, strsplit(res_vec$mut1, '_')[[1]])
      colnumbers_mut2 = as.numeric(colnames(bitlocus_mut))
      
      # mut 2
      m2_places = as.numeric(strsplit(res_vec$mut2, '_')[[1]][1:2]) 
      start_mean = round(sum(colnumbers_mut2[1:(m2_places[1] -1)]) + (colnumbers_mut2[m2_places[1]] / 2))
      end_mean = round(sum(colnumbers_mut2[1:(m2_places[2] -1)]) + (colnumbers_mut2[m2_places[2]] / 2))
      len_range = c(
        (end_mean-end_pm)-(start_mean+start_pm), # Shortest possible
        (end_mean+end_pm)-(start_mean-start_pm) # Largest possible
      )
      
      res_out[nrow, 'mut2_len'] = mean(len_range)
      res_out[nrow, 'mut2_len_pm'] = len_range[2] - mean(len_range)
      
      if ('mut3' %in% colnames(res_vec[,colSums(is.na(res_vec))<nrow(res_vec)])){
        bitlocus_mut_3 = carry_out_compressed_sv(bitlocus_mut, strsplit(res_vec$mut2, '_')[[1]])
        colnumbers_mut3 = as.numeric(colnames(bitlocus_mut_3))
        
        # mut 2
        m3_places = as.numeric(strsplit(res_vec$mut3, '_')[[1]][1:2]) 
        start_mean = round(sum(colnumbers_mut3[1:(m3_places[1] -1)]) + (colnumbers_mut3[m3_places[1]] / 2))
        end_mean = round(sum(colnumbers_mut3[1:(m3_places[2] -1)]) + (colnumbers_mut3[m3_places[2]] / 2))
        len_range = c(
          (end_mean-end_pm)-(start_mean+start_pm), # Shortest possible
          (end_mean+end_pm)-(start_mean-start_pm) # Largest possible
        )
        
        res_out[nrow, 'mut3_len'] = mean(len_range)
        res_out[nrow, 'mut3_len_pm'] = len_range[2] - mean(len_range)
      }
    }
  }
  
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
