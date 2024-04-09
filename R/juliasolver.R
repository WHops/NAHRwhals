#' @export
solve_mutation_julia_wrapper <- function(params, mat, gridlines_x, inmat_path, inlens_path, juliares_path, julia_solver_path = NULL){

  mat <- flip_bitl_y_if_needed(mat)
  write.table(mat, inmat_path, col.names = F, row.names=F, quote=F, sep='\t')
  write.table(t(as.numeric(row.names(mat))), inlens_path, col.names = F, row.names=F, quote=F, sep='\t')
  write.table(t(as.numeric(colnames(mat))), inlens_path, append=T, col.names = F, row.names=F, quote=F, sep='\t')

  possible_moves = find_sv_opportunities(mat)
  possible_dups = possible_moves[possible_moves$sv == 'dup',]
  longest_dup = max(possible_dups$p2 - possible_dups$p1)

  # Berzerk matrices are those that can explode due to the number of possible duplications.
  # These conditions will trigger very rarely. 
  breakcondition1 = dim(find_sv_opportunities(mat))[1] > 400
  breakcondition2 = (dim(find_sv_opportunities(mat))[1] > 200) &
                    (longest_dup > 40) &
                    (nrow(possible_dups) > 40)

  if (breakcondition1 | breakcondition2){
    # Emergency break here for absoluvely berzerk large matrices. 
    
    # Inform the user what just happened and how we react
    message(paste0('Warning. Julia detected a matrix prone to duplication-explosion. 
    Setting the defalt parameter for the BFS init_width to 1/10 of default. 
    This is unlikely to have an effect on the result.'))
    params$init_width = params$init_width/10
  }
  run_solver_julia(params, inmat_path, inlens_path, juliares_path, julia_solver_path)
  return(format_julia_output(juliares_path, gridlines_x, params$depth))
}

#' @export
run_solver_julia <- function(params, inmat_path, inlens_path, outpath, julia_solver_path){
  # Run julia
  julia_path = 'julia'

  if (is.null(params$julia_solver_path)){
    julia_solver_path = system.file('extdata', 'scripts', 'scripts_nw_main', 'solver.jl', package='nahrwhals')
  }

  maxdepth = params$depth
  maxdup = params$maxdup
  minreport = params$minreport
  init_width = params$init_width
# 
  julia_cmd = paste0(julia_path, " ", julia_solver_path, ' -d ', maxdepth, ' -u ', maxdup, ' -w ', init_width, ' -r ', minreport, ' -R 1000 ', inmat_path, ' ', inlens_path, ' ', outpath)
  if (params$silent == F){
    system(julia_cmd)
  } else {
    run_silent(julia_cmd)
  }
}

#' @export
format_julia_output <- function(juliares_path, gridlines_x, depth){
  # Load jout
  ncols = max(count.fields(juliares_path, sep='\t'))
  jout = read.table(juliares_path, sep='\t', header=F, fill = NA, col.names = paste0('V', 1:ncols))
  n_mut <- floor((ncol(jout) - 2) / 3)
  if (n_mut == 0){
    if (exists("log_collection")) {
      log_collection$depth <<- depth
      log_collection$mut_simulated <<- 0
      log_collection$mut_tested <<- 0
      
      return(data.frame(eval=jout[1], nmut = jout[2], mut1 = 'ref'))
    }
  }
  # Block 1: determine start, end, len coordinates
  jout_meta <- jout
  jout_meta[jout_meta == ''] = NA
  
  for (j in 1:nrow(jout_meta)){
    
    if (jout_meta[j,'V2'] == 0){
      next()
    }
    procrow = process_row(jout_meta[j,])
    df = turn_mut_max_into_svdf(procrow)
    for (i in 1:nrow(df)) {
      
      # Extract start and end positions from gridlines_x
      jout_meta[j,paste0('mut', i, '_start')] <- gridlines_x[df[i,'start']]
      jout_meta[j,paste0('mut', i, '_end')] <- gridlines_x[df[i,'end']]
      jout_meta[j,paste0('mut', i, '_len')] <- gridlines_x[df[i,'end']] - gridlines_x[df[i,'start']] 
    }
  }
  # Block 2: concatenate mutation names
  jout_new <- jout[, 1:2, drop = F]
  colnames(jout_new) <- c('eval', 'mut_max')
  jout_new$eval = as.numeric(jout_new$eval) 

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

  jout_combine$mut_max = NULL
  
  # Clean up
  #system(paste0('rm ', juliares_path))

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

#' turn_mut_max_into_svdf
#' A function to correct the coordinates of the mutations in the mut_max column
#' @export
turn_mut_max_into_svdf <- function(t, correction=T, return_string=F){
  if (t == 'ref'){
    return(F)
  }
  
  
  strsplit(t, split='n+')
  sv_list = strsplit(t, split="\\+")[[1]]
  svdf_clump = data.frame(sv_list)
  svdf = data.frame(do.call('rbind', strsplit(as.character(svdf_clump$sv_list),'_',fixed=TRUE)))
  colnames(svdf) = c('start', 'end', 'sv')
  
  svdf$start = as.numeric(svdf$start)
  svdf$end = as.numeric(svdf$end)
  svdf$depth = 1: dim(svdf)[1]
  
  
  if (correction == F){
    return(svdf)
  }
  
  if (nrow(svdf) == 1){
    if (return_string){
      out_string = paste0(svdf[1,1:3], sep='', collapse='_')
      return(out_string)
    }
  return (svdf)
  }

  # Iterate over each sv
  for (nsv in 1:(nrow(svdf)-1)){
    sv = svdf[nsv,'sv']
    sv_start = svdf[nsv,'start']
    sv_end = svdf[nsv,'end']
    depth = svdf[nsv,'depth']

    svdf_allstarts = svdf$start
    svdf_allends = svdf$end
    svdf_alldepths = svdf$depth

    # if it's an inv, complicated changes have to happen.
    if (sv == 'inv'){
           
      #start breakpoint gets mirrored
      svdf[(svdf_allstarts > sv_start) & 
           (svdf_allstarts < sv_end) & 
           (svdf_alldepths > depth),'newstart'] =
        sv_end - 
          (svdf[(svdf_allstarts > sv_start) & 
                (svdf_allstarts < sv_end) & 
                (svdf_alldepths > depth),'start']
                - sv_start)

      # end breakpoint
      svdf[(svdf_allends > sv_start) & 
           (svdf_allends < sv_end) & 
           (svdf_alldepths > depth),'newend'] =
        sv_end - 
          (svdf[(svdf_allends > sv_start) & 
                (svdf_allends < sv_end) & 
                (svdf_alldepths > depth),'end'] 
                - sv_start)
      
      newend_no_nas = !is.na(svdf$newend)
      newstart_no_nas = !is.na(svdf$newstart)
      svdf[newend_no_nas,'end'] = svdf[newend_no_nas,'newend']
      svdf[newstart_no_nas,'start'] = svdf[newstart_no_nas,'newstart']
      svdf[,c('newend', 'newstart')] = NULL

      # If in any row the start is larger than the end, we have to flip them.
      flip = svdf$start > svdf$end
      svdf[flip, c('start', 'end')] = svdf[flip, c('end', 'start')]

    }
    
    # If it's del/dup, change the del coordinates.
    if (sv == 'del'){
      
      #start breakpoints after the beginning of del get a plus
      svdf[(svdf_allstarts > sv_start) & (svdf_alldepths > depth),'start'] =
        svdf[(svdf_allstarts > sv_start) & (svdf_alldepths > depth),'start'] + (sv_end - sv_start)
      
      #end breakpoints after the beginning of del get a plus
      svdf[(svdf_allends > sv_start) & (svdf_alldepths > depth),'end'] =
        svdf[(svdf_allends > sv_start) & (svdf_alldepths > depth),'end'] + (sv_end - sv_start)
      
    } else if (sv == 'dup'){
      
      # start breakpoints after the dup get a minus
      svdf[(svdf_allstarts > sv_end) & (svdf_alldepths > depth),'start'] =
        svdf[(svdf_allstarts > sv_end) & (svdf_alldepths > depth),'start'] - (sv_end - sv_start)
      
      # end breakpoints after thhe dup get a minus
      svdf[(svdf_allends > sv_end) & (svdf_alldepths > depth),'end'] =
        svdf[(svdf_allends > sv_end) & (svdf_alldepths > depth),'end'] - (sv_end - sv_start)
    }
  }

  if (return_string){
    out_string = paste0(svdf[1,1:3], sep='', collapse='_')
    if (nrow(svdf) == 1){
      return(out_string)
    }
    for (i in 2:nrow(svdf)){
      out_string = paste(out_string, paste0(svdf[i,1:3], sep='', collapse='_'), sep='+')
    }
    return(out_string)
  }
  return (svdf)
}

# Function to process each row
process_row <- function(row) {
  # Skip first two columns and iterate over every three columns
  triplets <- sapply(seq(3, ncol(row), by = 3), function(i) {
    # Check if the values are NA and concatenate the triplet
    if (!is.na(row[i])) {
      move_part <- gsub("move_", "", row[i]) # Remove the "move_" prefix
      paste0(row[i+1], "_", row[i+2], "_", move_part) # Reorder the elements
    }
  })
  
  # Remove NA values and concatenate the triplets
  return(gsub("move_", "", paste(na.omit(unlist(triplets)), collapse = "+")))
}
