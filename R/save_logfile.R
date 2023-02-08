# Save logfile
#' I apologize for the existence of this terrible function. 
#' @export
save_to_logfile <- function(log, res, logfile, params, alt_x = F){
  res_ref = res[res$mut1=='ref','eval']
  res_max = res[1,'eval']
  n_res_max = dim(res[res$eval == max(res$eval),])[1]
  
  
  maxres =  res[
    which.max(
      rowSums(
        is.na(
          res[
            res$eval==max(res$eval),names(res) %in% c('mut1','mut2', 'mut3'), drop=F
            ]
          )
        )
      )
    ,]
  
  if (is.null(log$compression)){
    log$compression = 0
  }
  stopifnot(ncol(maxres) >= 2)

  
  # Find the result with the fewest steps.
  # If we have more than one mutation modeled
  if (ncol(maxres)>2){
    
    nmut = sum(c('mut1','mut2','mut3') %in% colnames(maxres[,colSums(is.na(maxres))<nrow(maxres)]))
    mut_cols = paste0('mut', 1:nmut)
    
    mut_max = paste0(maxres[,mut_cols], collapse='+')
    
  } else if (ncol(maxres) == 2){
    mut_max = paste0(maxres[1,2])
  }
  
  ###### TBD ###########
  to_append = data.frame(
    log$chr,
    log$start,
    log$end,
    log$samplename,
    as.numeric(log$end) - as.numeric(log$start),
    log$xpad,
    res_ref,
    res_max,
    n_res_max,
    mut_max,
    log$mut_simulated,
    log$mut_tested,
    log$depth,
    log$compression,
    log$exceeds_x,
    log$exceeds_y,
    log$grid_inconsistency,
    log$flip_unsure,
    log$cluttered_boundaries,
    maxres$mut1_start,
    maxres$mut1_end,
    maxres$mut1_pos_pm,
    maxres$mut1_len,
    maxres$mut1_len_pm,
    maxres$mut2_len,
    maxres$mut2_len_pm,
    maxres$mut3_len,
    maxres$mut3_len_pm
  )
  
  colnames_ = data.frame(
    'seqname',
    'start',
    'end',
    'sample',
    'width_orig',
    'xpad',
    'res_ref',
    'res_max',
    'n_res_max',
    'mut_max',
    'mut_simulated',
    'mut_tested',
    'search_depth',
    'grid_compression',
    'exceeds_x',
    'exceeds_y',
    'grid_inconsistency',
    'flip_unsure',
    'cluttered_boundaries',
    'mut1_start',
    'mut1_end',
    'mut1_pos_pm',
    'mut1_len',
    'mut1_len_pm',
    'mut2_len',
    'mut2_len_pm',
    'mut3_len',
    'mut3_len_pm'
  )

  if (alt_x){
    to_append = cbind(params$seqname_x_alt, params$start_x_alt, params$end_x_alt, to_append)
    colnames_ = cbind('seqname_alt', 'start_alt', 'end_alt', colnames_)
  }
  # If file doesn't exist, write colnames.
  if (is.na(file.size(logfile)) | (file.size(logfile) == 0)) {
    write.table(colnames_,file=logfile,append=TRUE, sep='\t', col.names = F, row.names = F, quote=F)
  }
  
  write.table(to_append,file=logfile,append=TRUE, sep='\t', col.names = F, row.names = F, quote=F)
  
  print('Logfile written.')
}


#' init_log_with_def_values
#' @export
init_log_with_def_values <- function(){
  
  log_collection = list()
  log_collection$exceeds_x = F
  log_collection$exceeds_y = F
  log_collection$grid_inconsistency = F
  log_collection$flip_unsure = F
  log_collection$cluttered_boundaries = F
  log_collection$mut_simulated = 0
  log_collection$mut_tested = 0
  
  return(log_collection)
}
