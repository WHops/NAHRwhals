#' Find parameters by wiggling.
#' 
#' @export
find_minlen_compression_params_wiggle <- function(inpaf, n_tests = 50, n_max_alns = 100, mode='precise', max_size_col_plus_rows=100){
  
  if (mode == 'precise'){
    quantile_preferred = 0.0
  } else if (mode == 'compressed'){
    quantile_preferred = 0.1
  }
  
  paf_loaded = read.table(inpaf)
  approx_size = (max(paf_loaded$V8) - min(paf_loaded$V8))
  
  baseline_minsize_min = max(log2(256), log2(approx_size/1000))
  baseline_minsize_max = max(log2(16384), log2(approx_size/10))

  res_interesting = data.frame()
  res = data.frame()
  
  nround = 0
  
  while (dim(res_interesting)[1] < 2){
    print('running...')
    # mode can be 'precise' or 'compressed'

    min_ = baseline_minsize_min + nround
    max_ = baseline_minsize_max 
    
    print(2**min_)
    print(2**max_)
    res = rbind(res, run_n_tests(inpaf, min_, max_, res, n_tests, n_max_alns, max_size_col_plus_rows))
    
    res_interesting = res[(res$success == T) & (res$n_rows_cols <= max_size_col_plus_rows),]
    
    if (dim(res_interesting)[1] < 2){
      
      n_max_alns = n_max_alns * 1.25
      max_size_col_plus_rows = max_size_col_plus_rows * 1.25
      nround = nround + 1
      res_interesting = res[(res$success == T) & (res$n_rows_cols <= max_size_col_plus_rows),]
      
  }
    res_to_take_idx = (which.min(abs(res_interesting$collapspenalty-quantile(res_interesting$collapspenalty, quantile_preferred))))[[1]]
    take = res_interesting[res_to_take_idx,]
    
    print(paste0('Choosing as ideal wiggle values: ', take$minlen, ' and ', take$compression ))
    return(list(minlen=take$minlen, compression=take$compression, res=res))
    
  }
}
  #else {
  #   print('Unable to find non-failing parameter. Try again with more tests or higher intervals. Using default values 1000, 1000.')
  #   return(list(minlen=1000, compression=1000, res=res))
  #   
  #}
#}


run_n_tests <- function(inpaf, min_, max_, res, n_tests, n_max_alns, max_size_col_plus_rows){
  
  counter = 0
  smallest_success_minlen = 2**max_
  for (i in 1:n_tests){
    minlen = as.integer(runif(1, min = 2**min_, max = smallest_success_minlen))
    compression = as.integer(runif(1, min = 2**min_, max = minlen))
    
    print(minlen)
    print(compression)
    # 
    # compression = min(x1,x2)
    # minlen = max(x1,x2)
    
    counter = counter + 1
    paf = load_and_prep_paf_for_gridplot_transform(inpaf, minlen= minlen, compression = compression, quadrantsize = 1e5)
    
    # If not too many alignments, construct a grid.
    # If the grid also not too many dimensions, we are happy. 
    n_alns_in_paf_ok = (dim(paf)[1] > 0) & (dim(paf)[1] <= n_max_alns)
    if (n_alns_in_paf_ok){
      gxy = make_xy_grid(paf, n_additional_bounces = 2)
      if ((length(gxy[[1]]) + length(gxy[[2]])) >= max_size_col_plus_rows){
        print('Failing: Exceeding dotplot dimensions.')
      } else if (gxy[[3]] == F){
        print('Failing: Condensed dotplot with incongruencies.')
      }
      res = rbind(res, data.frame(minlen=minlen, compression=compression, success=gxy[[3]], n_alns=dim(paf)[1], n_rows_cols=length(gxy[[1]]) + length(gxy[[2]])))
    } else {
      print('Too many alignments!')
    }
    
    
    res_success = res[(res$success == T) & (res$n_rows_cols <= max_size_col_plus_rows),]
    if (dim(res_success)[1] > 0){
      smallest_success_minlen = min(res_success[res_success$success==T,'minlen']) * 1.5
      print('New smallest minlen:')
      print(smallest_success_minlen)
    }
    
    
    print(paste0('Compression wiggle: ', counter, ' of ', n_tests ))
  }
  
  res$collapspenalty = ((res$minlen) + (res$compression)) # * (max(res$minlen, res$compression) / min(res$minlen, res$compression))

  return(res)
}
