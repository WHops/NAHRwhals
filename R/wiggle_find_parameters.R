#' Find parameters by wiggling.
#' This produces unstable and unpredictable results; use is discouraged but
#' we are leaving the function alive for now. Is called when compression_params$auto_find is set to True (default: False)
#'
#' @export
find_minlen_compression_params_wiggle <-
  function(inpaf,
           compression_params) {
    
    # (unfortunately) There is a random factor in this at the moment.
    # Would be better to go without randomness in the future.
    # set.seed(1)
    
    # Which mode are we operating in?
    if (compression_params$mode == 'precise') {
      quantile_preferred = 0.0
    } else if (compression_params$mode == 'compressed') {
      quantile_preferred = 0.1
    }
    
    # Load paf to get a size approximation of the alignment.
    paf_loaded = read.table(inpaf)
    approx_size = (max(paf_loaded$V8) - min(paf_loaded$V8))
    
    # Initialize dataframe to store minlen-compression-success info.
    res_interesting = data.frame()
    res = data.frame()
    
    # Keep count how many blocks of tries we need. Every round gets
    # relaxed criteria. 
    nround = 0
    
    # Do this until one round has yielded  at least 3 viable grids.
    while (dim(res_interesting)[1] < 2){

      # If we haven't found anything after e.g. 20 x 5 runs, return something suboptimal.
      if  (dim(res)[1] >= (compression_params$n_tests * compression_params$n_max_testchunks)){
        resp = res[res$n_rows_cols < Inf,]
        return(list(minlen=resp[resp$n_rows_cols == max(resp$n_rows_cols),]$minlen,
                    compression=resp[resp$n_rows_cols == min(resp$n_rows_cols),]$compression))
      } 
      
      # Increase lower cutoff with every round. We steer towards more coarse grids.
      min_ = compression_params$baseline_log_minsize_min + (2 * nround)
      max_ = max(compression_params$baseline_log_minsize_max)
      
      # Run the block of n tests (typically around 20, 50 or 100)
      res = rbind(
        res,
        run_n_tests(
          inpaf,
          min_,
          max_,
          res,
          compression_params$n_tests,
          compression_params$max_n_alns,
          compression_params$max_size_col_plus_rows
        )
      )
      # The interesting stuff is what was 'successful' (no incongruencies)
      # and is not too big. 
      res_interesting = res[(res$success == T) &
                              (res$n_rows_cols <= compression_params$max_size_col_plus_rows), ]
      
      # If we don't have enough good results, raise criteria now and see if 
      # it works then. This prevents us from rerunning everything. 
      if (dim(res_interesting)[1] < 2) {
        compression_params$max_n_alns = compression_params$max_n_alns * 1.25
        compression_params$max_size_col_plus_rows = compression_params$max_size_col_plus_rows * 1.25
        nround = nround + 1
        res_interesting = res[(res$success == T) &
                                (res$n_rows_cols <= compression_params$max_size_col_plus_rows), ]
        
      }
    }
      
    # If we are here, it means we have at least 3 viable grids. Take the one 
    # with the lowest collapse penalty.
    res_to_take_idx = (which.min(
      abs(res_interesting$collapspenalty - 
            quantile(res_interesting$collapspenalty, quantile_preferred)
            )))[[1]]
    take = res_interesting[res_to_take_idx, ]
    
    print(paste0(
      'Choosing as ideal wiggle values: ',
      take$minlen,
      ' and ',
      take$compression
    ))
    
    return(list(
      minlen = take$minlen,
      compression = take$compression
    ))
}



#' TODO: describe
#' @export
run_n_tests <-
  function(inpaf,
           min_,
           max_,
           res,
           n_tests,
           n_max_alns,
           max_size_col_plus_rows) {

    upper_bound_minlen = 2 ** max_
    smallest_fail_minlen = 1e10
    lower_bound_minlen = 2 ** min_
    for (i in 1:n_tests) {
      
      print(paste0('Compression wiggle: ', i, ' of ', n_tests))
      
      if (upper_bound_minlen < lower_bound_minlen) {
        lower_bound_minlen = 2 ** min_
      }
      minlen = as.integer(runif(1, min = lower_bound_minlen, max = upper_bound_minlen))
      compression = as.integer(runif(1, min = lower_bound_minlen, max = minlen))
      print(paste0("Trying parameters: ", minlen, ' and ', compression))
      paf = load_and_prep_paf_for_gridplot_transform(
        inpaf,
        minlen = minlen,
        compression = compression,
        quadrantsize = 1e5
      )
      
      print('hi')
      # If not too many alignments, construct a grid.
      # If the grid also not too many dimensions, we are happy.
      n_alns_in_paf_ok = (dim(paf)[1] > 0) &
        (dim(paf)[1] <= n_max_alns)
      print(dim(paf))

      if (n_alns_in_paf_ok) {
        gxy = make_xy_grid(paf, n_additional_bounces = 2)
        if ((length(gxy[[1]]) + length(gxy[[2]])) >= max_size_col_plus_rows) {
          print('Failing: Exceeding dotplot dimensions.')
          if (minlen > lower_bound_minlen) {
            lower_bound_minlen = minlen
            print(paste0('Updating lower_bound minlen to ', minlen))
          }
        } else if (gxy[[3]] == F) {
          print('Failing: Condensed dotplot with incongruencies.')
        }
        res = rbind(
          res,
          data.frame(
            minlen = minlen,
            compression = compression,
            success = gxy[[3]],
            n_alns = dim(paf)[1],
            n_rows_cols = length(gxy[[1]]) + length(gxy[[2]]),
            collapspenalty = Inf
          )
        )
      } else if (dim(paf)[1] > n_max_alns){
        print('Too many alignments!')
        if (minlen > lower_bound_minlen) {
          lower_bound_minlen = minlen
          print(paste0('Updating lower_bound minlen to ', minlen))
          res = rbind(res,
                      data.frame(
                        minlen = minlen,
                        compression = compression,
                        success = F,
                        n_alns = dim(paf)[1],
                        n_rows_cols = Inf,
                        collapspenalty = Inf
                      ))
          
        }
      } else if (dim(paf)[1] == 0){
        print('No alignments! Continue searching in larger areas')
        upper_bound_minlen = minlen
        res = rbind(res,
                    data.frame(
                      minlen = minlen,
                      compression = compression,
                      success = F,
                      n_alns = dim(paf)[1],
                      n_rows_cols = Inf,
                      collapspenalty = Inf
                      
                    ))
      } else {
      }
      
      
      res_success = res[(res$success == T) &
                          (res$n_rows_cols <= max_size_col_plus_rows), ]
      if (dim(res_success)[1] > 0) {
        upper_bound_minlen = min(res_success[res_success$success == T, 'minlen']) * 1.2
        print(paste0(
          'Updating upper_bound minlen to ',
          upper_bound_minlen
        ))
      }
      
    }
    
    res$collapspenalty = (res$minlen) + (res$compression)
    
    return(res)
  }
