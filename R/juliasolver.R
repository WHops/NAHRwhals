#' @export
solve_mutation_julia_wrapper <- function(params, mat, gridlines_x, inmat_path, inlens_path, juliares_path){
  
  mat <- flip_bitl_y_if_needed(mat)
  write.table(mat, inmat_path, col.names = F, row.names=F, quote=F, sep='\t')
  write.table(t(as.numeric(row.names(mat))), inlens_path, col.names = F, row.names=F, quote=F, sep='\t')
  write.table(t(as.numeric(colnames(mat))), inlens_path, append=T, col.names = F, row.names=F, quote=F, sep='\t')
  
  run_solver_julia(params, inmat_path, inlens_path, juliares_path)
  return(format_julia_output(juliares_path, gridlines_x, params$depth))
}

#' @export
run_solver_julia <- function(params, inmat_path, inlens_path, outpath){
  # Run julia
  julia_path = params$julia_bin  #'/Users/hoeps/.juliaup/bin/julia'
  solver_path = params$solverscript #'~/Programs/mutation-matrices/src/solver.jl'
  maxdepth = params$depth
  maxdup = params$maxdup
  minreport = params$minreport
  init_width = params$init_width

  julia_cmd = paste0(julia_path, " ", solver_path, ' -d ', maxdepth, ' -u ', maxdup, ' -w ', init_width, ' -r ', minreport, ' ', inmat_path, ' ', inlens_path, ' ', outpath)
  system(julia_cmd)
}

#' @export
format_julia_output <- function(juliares_path, gridlines_x, depth){
  # Load jout
  ncols = max(count.fields(juliares_path, sep='\t'))
  jout = read.table(juliares_path, sep='\t', header=F, fill = NA, col.names = paste0('V', 1:ncols))

  n_mut <- (ncol(jout) - 2) / 3

  # Block 1: determine start, end, len coordinates
  jout_meta <- jout
  for(i in 1:n_mut) {
    
    # Start and End indices columns for each mutation triplet
    start_col <- 1 + (3 * i)
    end_col <- 2 + (3 * i) 
    
    # Extract start and end positions from gridlines_x
    jout_meta[paste0('mut', i, '_start')] <- gridlines_x[jout[, start_col]]
    jout_meta[paste0('mut', i, '_end')] <- gridlines_x[jout[, end_col]]
    jout_meta[paste0('mut', i, '_len')] <- jout_meta[paste0('mut', i, '_end')] - jout_meta[paste0('mut', i, '_start')]
  }
  
  # Block 2: concatenate mutation names
  jout_new <- jout[, 1:2, drop = F]
  colnames(jout_new) <- c('eval', 'mut_max')
  jout_new$eval = as.numeric(jout_new$eval) * 100 

  # Iterate over the number of 'mut' columns and generate the concatenated columns
  for(i in 1:n_mut) {
    # Get the next three columns, concatenate them using the custom function, and add to jout_new
    jout_new[paste('mut', i, sep='')] <- apply(jout[ , (3*(i-1) + 3):(3*i + 2)], 1, concat_triplet)
  }
  
  # Stitch together Blocks 1 and 2
  jout_combine = cbind(jout_new, jout_meta[, (2 + (3*n_mut) + 1):ncol(jout_meta)])
  
  # Sort ascending by first column and ascending by the number of NA per line
  jout_combine <- jout_combine[order(-jout_combine$eval, -rowSums(is.na(jout_combine))), ]
  
  if (exists("log_collection")) {
    log_collection$depth <<- depth
    log_collection$mut_simulated <<- dim(jout_combine)[1]
    log_collection$mut_tested <<- dim(jout_combine)[1]
  }
  
  jout_combine[jout_combine$mut_max == 0, 'mut1'] = 'ref'
  if (!'ref' %in% jout_combine$mut1){
    jout_combine = rbind(jout_combine, c(1, 0, 'ref', rep(NA,ncol(jout_combine)-3)))
  }
  jout_combine$mut_max = NULL
  
  return(jout_combine)
}



#' @export
concat_triplet <- function(row) {
  if(any(is.na(row))) {
    return(NA)
  }
  concatenated <- paste0(row[c(2,3,1)], collapse='_')
  # Remove all spaces
  concatenated <- gsub(" ", "", concatenated)
  concatenated <- gsub("move_", "", concatenated)
  return(concatenated)
}
