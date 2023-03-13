#' decide_loop_continue_symmetry
#' 
#' In the tree search, we sometimes want to stop exploring an outcome if the dotplot is very unsymmetrical.
#' Symmetricity (is that a word?) is one of the simplest proxy for deciding if a dotplot is high or low score. 
#' @param bitl_f a condensed bitlocus as a matrix
#' @param symm_cutoff a value between 0 and 1 for lowest acceptable symmetry (default: 0.8)
#' @param orig_symm remnant from debugging, don't touch please :) 
#' @export
decide_loop_continue_symmetry <- function(bitl_f, symm_cutoff = 0.80, orig_symm = 1){
  
  climb_up_cost = as.numeric(row.names(bitl_f))
  walk_right_cost = as.numeric(colnames(bitl_f))
  symmetry = min(sum(climb_up_cost), sum(walk_right_cost)) / max(sum(climb_up_cost), sum(walk_right_cost))
  
  # Run away if there are at least 5 columns, and we have less than 75% symmetry
  if ((symmetry < (orig_symm * symm_cutoff)) & (dim(bitl_f)[1] > 5)){
    return(F)
  } 
  
  return(T)
}



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
#' 
#' To speed up grid search,, we keep an inventory of bitloci which we have seen. 
#' For this, we simplify the matrix ('return diag values') and then turn this in to a hash, 
#' which we save and can use for comparisons later. 
#' 
#' @param bitlocus input bitlocus; a matrix
#' @pairs a data frame (?) containing known hashes. 
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
    pairs[npair, 'hash'] = rlang::hash(return_diag_values_new(carry_out_compressed_sv(bitlocus, pairs[npair,1:3]), fraction = 0.05) )
    #pairs[npair, 'hash'] = rlang::hash(carry_out_compressed_sv(bitlocus, pairs[npair,1:3]))
    #pairs[npair, 'hash'] = rlang::hash(as.numeric(carry_out_compressed_sv(bitlocus, pairs[npair,1:3])))
    
  }
  
  return(pairs)
}

#' Diagonalize
#' 
#' Takes a matrix, turns the corners into zero. 
#' @param matrix a matrix typically representing a bitlocus
#' @param factor value between 0 and 1 determining how much to erase from the matrix.
#' @export
diagonalize <- function(matrix, factor = 0.25){
  return(matrix *  
            (  
              (abs(row(matrix) - col(matrix))) < (ncol(matrix) * factor) 
              )    
         )

}

#' Diagonalize.
#' 
#' Return values near the diagonal of a matrix.
#' @param matrix a matrix typically representing a bitlocus
#' @param fraction, how 'thick' is the piece that we are reporting (0 = nothing, 1 = whole matrix)
#' @export
return_diag_values_new <- function(matrix, fraction=0.05){
  
  max_dist = round(max(5, fraction * nrow(matrix)))
  
  # Bottom left to top right
  #diag_idx = seq(1, min(dim(matrix))**2, min(dim(matrix)) + 1)
  diag_idx = seq(1, ((dim(matrix)[1]+1) * (min(dim(matrix))-1)) + 1, dim(matrix)[1]+1)
  
  idxs = outer(diag_idx, seq(-max_dist, max_dist, 1), FUN = "+")
  idxs_clean = idxs[(idxs > 0) & (idxs <= ((dim(matrix)[1]+1) * (min(dim(matrix))-1)) + 1)]
  idxs_clean_2 = (dim(matrix)[1] * dim(matrix)[2]) - idxs_clean
  idxs_clean_all = c(idxs_clean, idxs_clean_2)
  
  # Top right to bottom left
  return(matrix[idxs_clean_all])
}


#' transform_res_new_to_old
#' 
#' This is a helperfunction transforming between different ways of representing the results of a tree search.
#' Further documentation required - for now, just don't touch it :). 
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



#' sort_new_by_penalised
#' 
#' After computing results, we want to penalize results that needed more steps; e.g, if a depth 1 and a depth 3 result
#' gave the same result, we want the depth 1 result to appear first. This also helps avoiding over-correcting of the model, 
#' if it e.g. uses 3 steps to correct a small alignment artifact. Each step therefore gets a penalty that is 1.5 times the size
#' of the smallest dot of a bitlocus. 
#' 
#' @param bitlocus_f matrix, representing a bitlocus
#' @param res_df_f  a data frame containing results from a tree search run. 
#' 
#' @export
sort_new_by_penalised <- function(bitlocus_f, res_df_f){
  
  # Just to be on the safe side here
  res_df_f$eval = as.numeric(res_df_f$eval)
  
  min_dotsize = min(as.numeric(row.names(bitlocus_f)), as.numeric(colnames(bitlocus_f)))
  step_penalty = (min_dotsize / sum(as.numeric(colnames(bitlocus_f))) * 100) * 1.5
  res_df_f$eval_penalised = as.numeric(res_df_f$eval) - (as.numeric(res_df_f$depth) * step_penalty)
  
  res_df_f = res_df_f[order(as.numeric(res_df_f$eval_penalised), decreasing=T ),]
  
  return(res_df_f)
}



#' Which percentage do we call 'solved'? 
#' To account for alignment incongruencies, we say that a puzzle is solved
#' If it deviates at most by 1 'unit'. 
#' 
#' I want to have a second look at this funciton... is it clean / reasonable?
#' 
#' @export
find_threshold <- function(bitlocus_f, nmut){
  
  min_dotsize = min(as.numeric(row.names(bitlocus_f)), as.numeric(colnames(bitlocus_f)))
  step_penalty = (min_dotsize / sum(as.numeric(colnames(bitlocus_f))) * 100) 
  
  threshold = 100 - (step_penalty * nmut)
  
  return(threshold)
}


#' Determine if a bitlocus is good for use in the first place. 
#' If it has many dots around its edges, it is not good and we
#' need a larger window.  
#' 
#' I want to have a second look at this function... is it clean / reasonable?
#' 
#' @export
is_cluttered <- function(bitlocus_f, clutter_limit_per_border = 5){
  is_large = (dim(bitlocus_f)[1] > 5 & dim(bitlocus_f)[2] > 5)
  
  if (!is_large){
    return(F)
  }
  
  
  clutter_firstrow = sum(bitlocus_f[1                 ,2:(dim(bitlocus_f)[2]-1)] != 0)
  clutter_lastrow = sum(bitlocus_f[dim(bitlocus_f)[1],2:(dim(bitlocus_f)[2]-1)] != 0)
  clutter_firstcol = sum(bitlocus_f[2:(dim(bitlocus_f)[1]-1),1] != 0)
  clutter_lastcol = sum(bitlocus_f[2:(dim(bitlocus_f)[1]-1),dim(bitlocus_f)[2]] != 0)
  clutter_all = c(clutter_firstrow, clutter_lastrow, clutter_firstcol, clutter_lastcol)
  
  if (any(clutter_all >= clutter_limit_per_border )){
    print('Cluttered alignment. Adding this info output.')
    if (exists('log_collection')){
      log_collection$cluttered_boundaries <<- T
    }
    return(T)
  }
  
  return(F)
}


#' Determine if a PAF DOTPLOT is good for use in the first place. 
#' If it has many dots around its edges, it is not good and we
#' need a larger window.  
#' 
#' I want to have a second look at this function... is it clean / reasonable?
#' 
#' @export
is_cluttered_paf <- function(paf, clutter_limit_per_border = 5){
  
  clutter1 = sum(paf$qstart == 0)
  clutter2 = sum(paf$qend == max(paf$qend))
  clutter3 = sum(paf$tstart == 0)
  clutter4 = sum(paf$tend == max(paf$tend))
  clutter_all = c(clutter1, clutter2, clutter3, clutter4)
  
  # Second line of evidence for a cluttered alignment.
  seems_simplistic = 
    (length(unique(paf$qstart)) < (nrow(paf) / 2)) | 
    (length(unique(paf$tstart)) < (nrow(paf) / 2)) | 
    (length(unique(paf$qend))   < (nrow(paf) / 2)) | 
    (length(unique(paf$qstart)) < (nrow(paf) / 2)) 
    

  
  if (any(clutter_all >= clutter_limit_per_border ) & seems_simplistic){
    print('Cluttered alignment. Adding this info output.')
    if (exists('log_collection')){
      log_collection$cluttered_boundaries <<- T
    }
    return(T)
  }
  
  return(F)
  
}


# What is the initial symmetry of the bitlocus? 
#' 
#' @params bitl matrix representing a bitlocus
#' @export
calc_symm <- function(bitl){
  climb_up_cost = as.numeric(row.names(bitl))
  walk_right_cost = as.numeric(colnames(bitl))
  orig_symm = min(sum(climb_up_cost), sum(walk_right_cost)) / max(sum(climb_up_cost), sum(walk_right_cost))
  
  return(orig_symm)
}


#' reduce_depth_if_needed
#' 
#' Called at the beginning of a tree search. If the dotplot is unreasonably large, 
#' we pragmatically reduce the search depth here to prevent overly long runtimes. 
#' @export
reduce_depth_if_needed <- function(bitlocus, increase_only, maxdepth){
  # Prepare bitlocus that we will be working on. 
  bitl = flip_bitl_y_if_needed(bitlocus)
  
  # Decide if we should go forward. 
  n_pairs = dim(find_sv_opportunities(bitl))[1]
  if ((n_pairs > 200) & (maxdepth == 3) & (increase_only==F)){
    print('Uh oh that is a bit large. Reducing depth to 2. Try to avoid producing such large alignments.')
    maxdepth = 2
  }
  
  # [W, 30th Nov 2022]
  # DELME WHEN DONE. 
  # should be 800! Moved to 200 for testing purposes only!!
  if ((n_pairs > 600) & (increase_only==F)){
    print('Huge Alignment! Going for depth 1. ')
    maxdepth = 1
  }
  
  return(maxdepth)
}



#' run_steepest_descent_one_layer_down
#' 
#' Helperfunction to run a tree search in steepest_descent mode. 
#' @export
run_steepest_descent_one_layer_down <- function(bitlocus, res_df, starting_depth = 1){
  
  # Get the top 10 results with depth 2 and individual hashes
  res_df2 = res_df[(res_df$depth == starting_depth), ]
  res_df2_sort = res_df2[order(res_df2$eval, decreasing = T), ]
  res_df2_sort = res_df2_sort[!duplicated(res_df2_sort$hash),]
  res_d2_top10 = na.omit(res_df2_sort[1:10,])
  
  # Keep track of maxres, which decides if a new result is worth saving.
  maxres = max(res_d2_top10$eval)

  # Iterate over all of those 'best' depth=2 results
  for (n_res in row.names(res_d2_top10)){
    
    thisres = res_d2_top10[n_res,]
    mut1 = strsplit(strsplit(thisres$mut_path, '\\+')[[1]][[2]], '_')[[1]]
    mut1_df = data.frame(p1=mut1[1],p2=mut1[2],sv=mut1[3])
    
    if (starting_depth == 2){
      mut2 = strsplit(strsplit(thisres$mut_path, '\\+')[[1]][[3]], '_')[[1]]
      mut2_df = data.frame(p1=mut2[1],p2=mut2[2],sv=mut2[3])
      
      # Execute the two mutations
      bitlocus_new = carry_out_compressed_sv(carry_out_compressed_sv(bitlocus, mut1), mut2)
    } else if (starting_depth == 1){
      # Execute one mutation
      bitlocus_new = carry_out_compressed_sv(bitlocus, mut1)
    }
    
    
    # Run depth = 1 run for this one. ('history' tells it that its from depth 2.)
    vis_list = (dfs(bitlocus_new, maxdepth = starting_depth + 1, increase_only = F, earlystop = 1e10, history=thisres))

    
    res_df_single = vis_list[[2]]
    
    res_df = rbind(res_df, res_df_single[res_df_single$eval >= maxres,])
    
    # Are we happy with the best result? 
    solve_th = find_threshold(bitlocus, starting_depth + 1)
    conclusion_found = (max(as.numeric(res_df$eval)) >= solve_th)
    
    if (conclusion_found){
      print("Conclusion found!!")
    }
    
  }
  
  return(res_df)
  
}


#' annotate_res_out_with_positions_lens
#' 
#' A helper function that takes a results dataframe and an original bitlocus, 
#' and translates the results (e.g. '4-10-inv') into coordinates (e.g. 450000 - 600000 - inv).
#' This is quite limited and rudimentary at this point. 
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
