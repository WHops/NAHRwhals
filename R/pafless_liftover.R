# Pafless liftover tools

#' Doc pending
#' @export
create_mmi_if_doesnt_exists <- function(params){
  
  if (file.exists(params$genome_y_fa_mmi)){
    print('Found existing minimap2 index ".mmi" file. Skipping re-calculation.')
    return()
  }
  
  minimap2_bin = '/Users/hoeps/Programs/minimap2-v2/minimap2-2.24/minimap2'  

  print('No minimap2 index ".mmi" file of the assembly fasta found. Creating one now (takes around 1 minute.)')
  
  minimap2_indexing_command = paste0(minimap2_bin, ' -k 28 -w 255 --idx-no-seq -H -d ', params$genome_y_fa_mmi, ' ', params$genome_y_fa) 
  system(minimap2_indexing_command)
  
}


#' Doc pending
#' @export
make_params_conversionpaf <- function(params, outlinks){
  
  minimap2_bin = '/Users/hoeps/Programs/minimap2-v2/minimap2-2.24/minimap2'  
  
  # Write hg38 to file
  extract_subseq_bedtools(params$genome_x_fa,
                          params$seqname_x,
                          params$start_x,
                          params$end_x,
                          outlinks$genome_x_fa_subseq)
  
  random_tag = as.character(runif(1, 1e10, 1e11))
  tmp_conversionpaf  = paste0('tmp_conversionpaf', random_tag, '.paf')

  minimap2_mapping_command = paste0(minimap2_bin, ' ', params$genome_y_fa_mmi, ' ', outlinks$genome_x_fa_subseq, ' > ', tmp_conversionpaf)
  system(minimap2_mapping_command)
  
  paf = read.table(tmp_conversionpaf)

  paf$start_add <- NA
  paf$end_add <- NA
  
  # loop through each row and extract the two numbers
  for (i in 1:nrow(paf)) {
    string <- paf$V1[i]
    # use regular expressions to extract the two numbers
    numbers <- as.numeric(unlist(stringr::str_extract_all(string, "\\d+")))
    paf$start_add[i] <- numbers[length(numbers)-1]
    paf$end_add[i] <- numbers[length(numbers)]
  }
  
  paf$V3 = paf$V3 + paf$start_add
  paf$V4 = paf$V4 + paf$end_add
  paf[,c('start_add', 'end_add' )] = NULL
  
  paf$V1 = params$seqname_x
  
  write.table(paf, file=tmp_conversionpaf, sep='\t', quote = F, row.names = F, col.names = F)

  params$conversionpaf_link = tmp_conversionpaf
  
  return(params)
}
