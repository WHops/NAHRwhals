#' filter_res
#' TODO: description
#' @export
filter_res <- function(res, threshold){
  
  # First, some border cases
  if (dim(res)[1] <= 1){
    return(res)
  }
  
  res$eval = as.numeric(res$eval)
  
  if (max(res$eval) < threshold){
    return(res)
  }
  
  # Needs any exceptions? 

  # The things R makes you do sometimes ...
  if ((dim(res)[2]-1) == 1){
    res$n_mut = 1 - (is.na(res[,2]) + (res[,2]=='ref'))
  } else {
    res$n_mut = (dim(res)[2]-1) - (rowSums(is.na(res[,2:dim(res)[2]]) + (res[,2:dim(res)[2]]  == 'ref') %>% replace(is.na(.), F)))
  }
  res$eval_accept = res$eval >= threshold
  res_sort = res[order(-res$eval_accept, res$n_mut, -res$eval ),]
  
  res_sort[,c('n_mut','eval_accept')] = NULL
  return(res_sort)
}

res = data.frame(eval=c(99.688, 99.688), mut1=c('ref','ref'))
filter_res(res, 98)

res = data.frame(eval=c(99.688, 99.688, 98, 98.1), mut1=c('ref','ref', 'inv', 'del'))
filter_res(res, 98)


res = data.frame(eval=c(91, 95.5, 94.4, 97, 97.0), mut1=c('ref','ref', 'inv', 'del', 'dup'), mut2=c(NA,NA,'inv',NA, NA))
filter_res(res, 98)

res = data.frame(eval=c(100, 82), mut1=c('2_inv','ref'), mut2 = c('dup',NA))
filter_res(res, 98)
