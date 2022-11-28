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
    pairs[npair, 'hash'] = rlang::hash(as.numeric(keep_only_diagonal_of_bitlocus(carry_out_compressed_sv(bitlocus, pairs[npair,1:3]))))
    #pairs[npair, 'hash'] = rlang::hash(as.numeric(carry_out_compressed_sv(bitlocus, pairs[npair,1:3])))
    
  }
  
  return(pairs)
}


diagnonalize <- function(matrix, factor = 0.25){
  return(matrix *  
            (  
              (abs(row(matrix) - col(matrix))) < (ncol(matrix) * factor)  
              )    
         )
}


keep_only_diagonal_of_bitlocus <- function(bitlocus_f){
  
  max_dist = max(5, 0.25 * nrow(bitlocus_f))
  inds = which(bitlocus_f != 0, arr.ind=T)
  inds_die = inds[abs(inds[,1] - inds[,2]) > max_dist,]
  bitlocus_f[inds_die] = 0
  
  # for (i in 2:(nrow(bitlocus_f)-1)){
  #   for (j in 2:(ncol(bitlocus_f)-1)){
  #     if (sum(bitlocus_f[(i-1):(i+1),(j-1):(j+1)] != 0) < 2){
  #       mask[i,j] = 0
  #     }
  #   }
  # }
  
  return(bitlocus_f)
}


#' Transformation block. Need to document / figure out what exactly it does...
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
#' I want to have a second look at this funciton... is it clean / reasonable?
#' @export
sort_new_by_penalised <- function(bitlocus_f, res_df_f){
  
  #
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


#' Determine if a dotplot is good for use in the first place. 
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
  
  clutter1 = sum(bitlocus_f[1                 ,2:(dim(bitlocus_f)[2]-1)] != 0)
  clutter2 = sum(bitlocus_f[dim(bitlocus_f)[1],2:(dim(bitlocus_f)[2]-1)] != 0)
  clutter3 = sum(bitlocus_f[2:(dim(bitlocus_f)[1]-1),1] != 0)
  clutter4 = sum(bitlocus_f[2:(dim(bitlocus_f)[1]-1),dim(bitlocus_f)[2]] != 0)
  clutter_all = c(clutter1, clutter2, clutter3, clutter4)
  
  if (any(clutter_all >= clutter_limit_per_border )){
    return(T)
  }
  
  return(F)
}

# What is the initial symmetry of the bitlocus? 
#' 
#' @export
calc_symm <- function(bitl){
  climb_up_cost = as.numeric(row.names(bitl))
  walk_right_cost = as.numeric(colnames(bitl))
  orig_symm = min(sum(climb_up_cost), sum(walk_right_cost)) / max(sum(climb_up_cost), sum(walk_right_cost))
  
  return(orig_symm)
}

#' 
#' @export
reduce_depth_if_needed <- function(bitlocus, increase_only, maxdepth){
  # Prepare bitlocus that we will be working on. 
  bitl = flip_bitl_y_if_needed(bitlocus)
  
  # Decide if we should go forward. But it's not a good place here. 
  n_pairs = dim(find_sv_opportunities(bitl))[1]
  if ((n_pairs > 200) & (maxdepth == 3) & (increase_only==F)){
    print('Uh oh that is a bit large. Reducing depth to 2. Try to avoid producing such large alignments.')
    maxdepth = 2
  }
  
  if ((n_pairs > 500) & (increase_only==F)){
    print('Huge Alignment! Going for depth 1. ')
    maxdepth = 1
  }
  
  return(maxdepth)
}
