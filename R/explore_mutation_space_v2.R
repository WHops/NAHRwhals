#' Alternative wrapper for calling dfs several times.
#'
#' This function runs a tree search with varying depths to find the best mutation for the given bit locus.
#' If the dotplot is unreasonably large, it reduces the search depth to prevent long runtimes.
#' If a conclusion is not found, it runs a steepest descent bonus round.
#'
#' @param bitlocus A bit locus to be worked on
#' @param maxdepth_input An integer indicating the maximum search depth
#' @param earlystop A number indicating the maximum evaluations
#' @param solve_th Threshold for considering a mutation 'solved', in percent (def. 98)
#' @return A data frame that stores the search results
#'
#' @export
solve_mutation <- function(bitlocus, maxdepth_input, earlystop = Inf, solve_th = 98, compression=1000, is_cluttered_already_paf=F) {
  #
  # browser()
  # bitlocus = readRDS('bitme_large')
  # maxdepth_input = 3
  # earlystop = Inf
  # ref_mut_pair = data.frame(p1='12', p2='25', sv='inv')
  # bitlocus = carry_out_compressed_sv(bitlocus, ref_mut_pair)
  est_highest <<- 0
  n_eval_total <<- 0
  n_eval_calc <<- 0
  n_hash_excluded <<- 0
  n_discontinued <<- 0

  # reject bitloci with cut-off alignments at the borders ('clutter')
  if (is_cluttered(bitlocus)){#} | is_cluttered_already_paf) {
    print("Alignments overlap with borders. Make frame larger!")
    res_out <- data.frame(eval = 0, mut1 = "ref")
    res_out <- annotate_res_out_with_positions_lens(res_out, bitlocus)

    return(res_out)
  }

  # Optional: adjust depth
  maxdepth <- reduce_depth_if_needed(bitlocus, increase_only = F, maxdepth_input)

  # Run iterations of increasing depth.
  # Always run the whole thing, but also don't go
  # deeper once a good solution is found.
  conclusion_found <- F
  current_depth <- 1
  
  # If the thingy is SO big that even the first stage will take very long,
  # Then we want to run only one slim run in which the 'ref' gets its eval. 
  if (maxdepth == 0){
    current_depth = 0
  }

  # The tree is explored in an iterative manner, in which all mutations from the one depth are
  # explored. If one depth does not yield a solution, the next depth is explored, until
  # the largest depth (3) is reached.
  # In cases where >1M mutations are expected in depth 3, (corresponding to 100 mutations on depth 1),
  # the tree is explored until depth 2, and only the 100 highest-scoring depth-2 mutations are followed up
  # in layer three. This approach cuts computation time drastically, but may miss mutation chains behind
  # local minima in rare instances. Additionally, it is not guaranteed that all alternative paths to an
  # optimal result are identified with this.
  while ((!conclusion_found) & (current_depth <= maxdepth)) {
    #print(paste0("Running depth layer: ", current_depth))

    # Run the dfs machinery (main work horse)
    vis_list <- (dfs(bitlocus, maxdepth = current_depth, increase_only = F, earlystop = earlystop))
    # Extract res_df
    res_df <- vis_list[[2]]

    # Are we happy with the best result?
    # solve_th = params$eval_th# find_threshold(bitlocus, current_depth)
    conclusion_found <- (max(as.numeric(res_df$eval)) >= solve_th)

    #if (conclusion_found) {
    #  print("SV found")
    #}

    # Higher depth for next iteration!
    current_depth <- current_depth + 1
  }

  # OPTIONAL: THIS IS THE STEEPEST DESCENT OF LAYER 3
  # IF max depth was reduced due to size, and if no conclusion was found,
  # then we go into a steepest-descent bonus round.
  if ((current_depth == 2) & (maxdepth < maxdepth_input) & (conclusion_found == F)) {
    res_df <- run_steepest_descent_one_layer_down(bitlocus, res_df, starting_depth = current_depth - 1)

    if (conclusion_found == F) {
      res_df <- run_steepest_descent_one_layer_down(bitlocus, res_df, starting_depth = current_depth)
    }
  }

  if ((current_depth == 3) & (maxdepth < maxdepth_input) & (conclusion_found == F)) {
    res_df <- run_steepest_descent_one_layer_down(bitlocus, res_df, starting_depth = current_depth - 1)
  }

  # Prepare output for returns
  res_df_sort <- sort_new_by_penalised(bitlocus, res_df, compression)
  res_out <- transform_res_new_to_old(res_df_sort)
  res_out$eval <- as.numeric(res_out$eval)

  res_out <- annotate_res_out_with_positions_lens(res_out, bitlocus)


  # Make an entry to the output logfile #
  if (exists("log_collection")) {
    log_collection$depth <<- maxdepth
    log_collection$mut_simulated <<- dim(res_df)[1]
    log_collection$mut_tested <<- dim(res_df)[1]
  }


  return(res_out)
}
