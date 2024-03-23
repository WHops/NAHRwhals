split_genotype_windows_above_max_size <- function(df, window_size = 5e6, overlap = 2.5e6) {
  # Initialize an empty dataframe for results
  results <- data.frame(chr = character(), start = numeric(), end = numeric(), stringsAsFactors = FALSE)
  
  for (i in 1:nrow(df)) {
    chr <- df$chr[i]
    start <- df$start[i]
    end <- df$end[i]
    
    if ((end - start) > window_size) {
      while (start < end) {
        new_end <- min(start + window_size, end)
        results <- rbind(results, data.frame(chr = chr, start = start, end = new_end))
        start <- start + window_size - overlap
      }
    } else {
      results <- rbind(results, df[i,])
    }
  }
  
  return(results)
}

#' @export
scan_for_windows <- function(ref_fa, asm_fa, threads = 1, outfile = NULL, sizes = c(200000, 1000000, 5000000 ), ref_masked_regions = 'none') {

    sizes <- sort(sizes, decreasing = TRUE)

    # Initialize a variable to keep track of the cumulative bed file name
    cumulative_bed <- NULL

    # outdirectory is the base of outfile
    outdirectory <- dirname(outfile)
    if (!dir.exists(outdirectory)) {
        dir.create(outdirectory, recursive = TRUE)
    }
    for (i in 1:length(sizes)) {
  
        size <- sizes[i]
        test_list <- wga_write_interval_list(ref_fa, asm_fa, outdirectory, outfile, size, 10000, ref_masked_regions, threads)
        system(paste0('cp ', outdirectory, '/list_cut_final.bed ', paste0(outdirectory, '/list_cut_final_', i, '.bed')))
        Sys.sleep(0.5)
        # If it's not the first size, perform intersection with the cumulative result
        if (!is.null(cumulative_bed)) {
            sys_command <- paste0('bedtools intersect -b ', cumulative_bed, ' -a ', paste0(outdirectory, '/list_cut_final_', i, '.bed'), ' -F 0.5 -v > ', paste0(outdirectory, '/list_cut_final_', i, '_noredund.bed'))
            system(paste0(sys_command))
            Sys.sleep(0.5)
            # Concatenate the non-redundant bed file of the current size with the cumulative bed file and sort
            sys_command_2 <- paste0("cat ", paste0(outdirectory, '/list_cut_final_', i, '_noredund.bed'), ' ', cumulative_bed, " | bedtools sort -i - > ", paste0(outdirectory, '/temp_list_cut_final.bed'))
            system(paste0(sys_command_2))
            Sys.sleep(0.5)
            # Now move the sorted results to the final file
            system(paste0('mv ', paste0(outdirectory, '/temp_list_cut_final.bed '), paste0(outdirectory, '/res_list_cut_final.bed')))
            Sys.sleep(0.5)
            # Update the cumulative bed file
            cumulative_bed <- paste0(outdirectory, '/res_list_cut_final.bed')
        } else {
        # For the first size, the current bed file is the cumulative file
          cumulative_bed <- paste0(outdirectory, '/list_cut_final_', i, '.bed')
        }
    }   
    system(paste0('mv ', cumulative_bed, ' ', outfile))

    windows_df = read.table(outfile, header = FALSE, sep = '\t', col.names = c('chr', 'start', 'end'))
    windows_df = split_genotype_windows_above_max_size(windows_df)
    write.table(windows_df, outfile, sep='\t', row.names = F, col.names = F, quote = F)
    return(windows_df)
    
}

