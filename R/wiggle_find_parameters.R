#' Find parameters by wiggling.
#' 
#' @export
find_minlen_compression_params_wiggle <- function(inpaf, n_tests = 50, n_max_alns = 100, mode='precise'){
  
  # mode can be 'precise' or 'compressed'
  if (mode == 'precise'){
    quantile_preferred = 0.0
  } else if (mode == 'compressed'){
    quantile_preferred = 0.5
  }
  
  counter = 0
  res = data.frame()
  
  for (i in 1:n_tests){
    x1 = as.integer(2**(runif(1, min = 8, max = 16)))
    x2 = as.integer(2**(runif(1, min = 8, max = 16)))
    
    compression = min(x1,x2)
    minlen = max(x1,x2)
    
    counter = counter + 1
    paf = load_and_prep_paf_for_gridplot_transform(inpaf, minlen= minlen, compression = compression, quadrantsize = 1e5)
    
    if ((dim(paf)[1] > 0) &  (dim(paf)[1] < n_max_alns)){
      gxy = make_xy_grid(paf, n_additional_bounces = 2)
      res = rbind(res, data.frame(minlen=minlen, compression=compression, success=gxy[[3]], n_alns=dim(paf)[1]))
    }
    print(paste0('Compression wiggle: ', counter, ' of ', n_tests ))
  }
  
  res$collapspenalty = (res$minlen + res$compression)# * (max(res$minlen, res$compression) / min(res$minlen, res$compression))
  res_interesting = res[res$success == T,]
  
  if (dim(res_interesting)[1] > 0){
    
    
    res_to_take_idx = (which.min(abs(res_interesting$collapspenalty-quantile(res_interesting$collapspenalty, quantile_preferred))))[[1]]
    take = res_interesting[res_to_take_idx,]
    
    print(paste0('Choosing as ideal wiggle values: ', take$minlen, ' and ', take$compression ))
    return(list(minlen=take$minlen, compression=take$compression, res=res))
    
  } else {
    print('Unable to find non-failing parameter. Try again with more tests or higher intervals. Using default values 1000, 1000.')
    return(list(minlen=1000, compression=1000, res=res))
    
  }
}