





library(hash)
library(microbenchmark)
library(devtools)
devtools::load_all()

# bitlocusfull = readRDS('~/Desktop/gridmatrix_sample')
#bitlocusfull = readRDS('~/Desktop/blub.tsv')
bitlocusfull = readRDS('~/Desktop/n15')

bitlocus = bitlocusfull#[15:35, 15:35]

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


annotate_pairs_with_hash <- function(bitlocus, pairs){
  
  # # border cases
  # if (dim(pairs)[1] == 0){
  #   return(pairs)
  # }
  
  pairs$hash = 'NA'
  for (npair in seq_along(pairs$hash)){
    pairs[npair, 'hash'] = rlang::hash(as.numeric(carry_out_compressed_sv(bitlocus, pairs[npair,1:3])))
  }
  
  return(pairs)
}


decide_loop_continue <- function(bitl_f, symm_cutoff = 0.80){
  
  climb_up_cost = as.numeric(row.names(bitl_f))
  walk_right_cost = as.numeric(colnames(bitl_f))
  symmetry = min(sum(climb_up_cost), sum(walk_right_cost)) / max(sum(climb_up_cost), sum(walk_right_cost))
  
  # Run away if there are at least 5 columns, and we have less than 75% symmetry
  if ((symmetry < symm_cutoff) & (dim(bitl_f)[1] > 5)){
    return(F)
  } 
  
  return(T)
}

dfsutil <- function(visited, pair, mutator, depth, maxdepth = 3, pairhistory=NULL){
  
  # Calc score of a node
  aln_score = calc_coarse_grained_aln_score(mutator)

  # Make an entry to the output df. 
  df_output <<- na.omit(rbind(df_output, c(pair[1,'hash'], depth, pairhistory, aln_score)))

  # Update the visited hash to reflect that we have seen and processed that node. 
  visited[[pair$hash]] = c(depth, pairhistory, aln_score)
  
  # Criteria to continue. 
  #   1) this is not the reference node
  #   2) If our node is perfect, no need to continue on that branch.
  #   3) If our node returns NA, we are somewhere weird
  #   4) on increasing depth, we have a cutoff crit for aln score. 
  
  # Make criterion 4:
  score_ref = as.numeric(df_output[df_output$mut_path == '1_1_ref',]$eval)
  min_score = c(score_ref * 0.6, score_ref * 0.8, score_ref)
  continue = (((pair[1,3] != 'ref') & 
               (aln_score < 100) & 
               !is.na(aln_score) & 
               (aln_score > min_score[depth]) &
               (depth < maxdepth)))

  # We caclualte the next batch of children
  if (continue){
    pairs = find_sv_opportunities(mutator)
    pairs = annotate_pairs_with_hash(mutator, pairs)
    
    
    for (npair in seq_along(pairs$hash)){
      # if(paste(pairhistory, paste0(pairs[npair,1:3], collapse = '_'), sep='+') == '16_20_dup+16_20_del'){
      #   browser()
      # }
      node_is_novel_bool = is.null(visited[[pairs[npair, 'hash']]])
      if (node_is_novel_bool){
        
        bitl_mut = carry_out_compressed_sv(mutator, pairs[npair,1:3])
        
        node_passes_symmetry_crit = decide_loop_continue(bitl_mut)
        if (node_passes_symmetry_crit){
          dfsutil(visited=visited, 
                  pair=pairs[npair,], 
                  mutator=bitl_mut, 
                  depth=depth+1, 
                  maxdepth=maxdepth,
                  pairhistory = paste(pairhistory, paste0(pairs[npair,1:3], collapse = '_'), sep='+'))
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
        df_output[df_output$hash == pairs[npair, 'hash'],c(2,3,4)] <<-  c(depth+1, 
                                         paste(pairhistory, paste0(pairs[npair,1:3], collapse = '_'), sep='+'), 
                                         visited[[pairs[npair, 'hash']]][3])
        
      }
    
    }
  }
  return(visited)
}


dfs <- function(bitlocus, maxdepth = 3, depth = 1){

  # Prepare bitlocus that we will be working on. 
  bitl = flip_bitl_y_if_needed(bitlocus)
  # Initialize a hash. 
  visited = hash()
  
  # Find children of the 'ref' node
  pairs = find_sv_opportunities(bitl)
  # Also make a 'ref' node to be calculated as well. 
  pairs = rbind(c(1,1,'ref'), pairs)
  
  # For each pair, calculate what the hash function of the child will be. 
  pairs = annotate_pairs_with_hash(bitl, pairs)

  # For every child...
  for (npair in seq_along(pairs$hash)){
    #browser()
    # Inform user where we are
    # if (npair == 54){
    #   browser()
    # }
    print(paste0('Calculating branch ', npair, ' out of ', dim(pairs)[1]))
    # if (paste(pairs[npair,1:3], collapse='_') == '8_15_inv'){
    #   browser()
    # } 
    node_is_novel_bool = is.null(visited[[pairs[npair, 'hash']]])
    

    # browser()
    if (node_is_novel_bool){
      bitl_mut = carry_out_compressed_sv(bitl, pairs[npair,1:3])
      
      # Check if symmetry of that node is acceptable. 
      node_passes_symmetry_crit = decide_loop_continue(bitl_mut)

      # If yes (or if this is the reference node, aka npair =1), continue.
      if (node_passes_symmetry_crit | (npair == 1)){
        visited = dfsutil(visited=visited, 
                          pair=pairs[npair,], 
                          mutator=bitl_mut, 
                          depth=depth, 
                          maxdepth=maxdepth,
                          pairhistory=paste(pairs[npair,1:3], collapse='_'))
      } 
      
    # If node is NOT novel, BUT has only been seen at higher depth,
    # update the path and depth. 
    } else if ((!node_is_novel_bool) &
              (visited[[pairs[npair, 'hash']]][1] > depth)){
        visited[[pairs[npair, 'hash']]][c(1,2)] = c(depth, 
                                                    paste(pairs[npair,1:3], collapse='_'))
        df_output[df_output$hash == pairs[npair, 'hash'],c(2,3,4)] <<- 
                                        c(depth, 
                                         paste(pairs[npair,1:3], collapse='_'), 
                                         visited[[pairs[npair, 'hash']]][3])
        
      }
  }
  print('sup')
  print('dfd')
  
  return(visited)
  
  
}

counter <<- 0
vis = NULL
df_output = data.frame(matrix(nrow = 1, ncol=4)) 

colnames(df_output) = c('hash','depth','mut_path','eval')
#vis = as.data.frame.hash(
vis = (dfs(bitlocus, maxdepth = 2))
#as.data.frame.hash(aa)
df_output = na.omit(as.data.frame(df_output))
df_output$eval = as.numeric(df_output$eval)

