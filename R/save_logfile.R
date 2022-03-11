# Save logfile
#' I apologize for the existence of this terrible function. 
#' Everything about it is very bad. I hope I will make it better
#' one day. 
#' 
#' @export
save_to_logfile <- function(log, res, logfile){
  
  res_ref = res[res$mut1=='ref','eval']
  res_max = res[1,'eval']
  n_res_max = dim(res[res$eval == max(res$eval),])[1]

  # Terrible line here. Terrible...
  mut_max = paste0(res[1,2:ncol(res)][,((res[1,2:ncol(res)] != "NA") & (!is.na((res[1,2:ncol(res)]))))], collapse='+')
  
  print('updated')
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
    log$depth,
    log$compression,
    log$exceeds_x,
    log$exceeds_y,
    log$grid_inconsistency,
    log$flip_unsure
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
    'search_depth',
    'grid_compression',
    'exceeds_x',
    'exceeds_y',
    'grid_inconsistency',
    'flip_unsure'
  )

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
  log_collection$mut_simulated = 0
  log_collection$mut_tested = 0
  
  return(log_collection)
}
