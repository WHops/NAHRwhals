#' Creates a minimap2 index (.mmi) file for a genome assembly if it doesn't exist.
#' @param params A list containing parameters for the function.
#' @export
create_mmi_if_doesnt_exists <- function(params) {
  if (file.exists(params$genome_y_fa_mmi)) {
    print('Found existing minimap2 index ".mmi" file. Skipping re-calculation.')
    return()
  }


  print('No minimap2 index ".mmi" file of the assembly fasta found. Creating one now (takes around 1 minute for a whole genome assembly.)')

  minimap2_indexing_command <- paste0(params$minimap2_bin, " -k 28 -w 255 --idx-no-seq -H -d ", params$genome_y_fa_mmi, " ", params$genome_y_fa)
  system(minimap2_indexing_command)
}


#' DOC PENDING.
#' @param params A list of parameters specifying the input and output file paths and other options
#' @param outlinks A list of output file paths
#'
#' @return params, a list.
#'
#' @author Wolfram Hoeps
#' @export
make_params_conversionpaf <- function(params, outlinks) {
  # Write hg38 to file
  extract_subseq_bedtools(
    params$genome_x_fa,
    params$seqname_x,
    params$start_x,
    params$end_x,
    outlinks$genome_x_fa_subseq,
    params
  )
  random_tag <- as.character(runif(1, 1e10, 1e11))
  tmp_conversionpaf <- paste0("tmp_conversionpaf", random_tag, ".paf")
  minimap2_mapping_command <- paste0(params$minimap2_bin, " ", params$genome_y_fa_mmi, " ", outlinks$genome_x_fa_subseq, " > ", tmp_conversionpaf)
  print("Attempting to locate input sequence homolog in y assembly... ")
  system(minimap2_mapping_command)
  
  if (file.info(tmp_conversionpaf)$size == 0){
    print('Sequence not found!')
    params$conversionpaf_link <- NA
    return(params)
  }
  
  paf <- read.table(tmp_conversionpaf, fill = T)

  # The following will help eliminate the 'second' translation step and thus cut runtime in half. 
  # shred_seq_bedtools(outlinks$genome_x_fa_subseq, paste0(outlinks$genome_x_fa_subseq,'_chunked.fa'), params$chunklen, params)
  # 
  # minimap2_mapping_command <- paste0(params$minimap2_bin, " ", params$genome_y_fa_mmi, " ", paste0(outlinks$genome_x_fa_subseq,'_chunked.fa'), " > ", tmp_conversionpaf)
  # print("Attempting to locate input sequence homolog in y assembly... ")
  # system(minimap2_mapping_command)
  # 
  # correct_paf(tmp_conversionpaf, paste0(tmp_conversionpaf, '_filter.paf'))
  # compress_paf_fnct(inpaf_link = paste0(tmp_conversionpaf, '_filter.paf'), outpaf_link =paste0(tmp_conversionpaf, '_filter_stitch.paf') , inparam_chunklen = params$chunklen)
  # paf2 = read_and_prep_paf(paste0(tmp_conversionpaf, '_filter.paf'))
  #                   
  # write.table(paf2, file='~/Desktop/paffast.paf', col.names=T, row.names=F, quote=F, sep='\t')                  
  #                   
  
  
  
  paf$start_add <- NA

  # loop through each row and extract the two numbers
  for (i in 1:nrow(paf)) {
    string <- paf$V1[i]
    # use regular expressions to extract the two numbers
    numbers <- as.numeric(unlist(stringr::str_extract_all(string, "\\d+")))
    paf$start_add[i] <- numbers[length(numbers) - 1]
  }

  paf$V3 <- paf$V3 + paf$start_add
  paf$V4 <- paf$V4 + paf$start_add
  paf[, c("start_add")] <- NULL

  paf$V1 <- params$seqname_x

  write.table(paf, file = tmp_conversionpaf, sep = "\t", quote = F, row.names = F, col.names = F)

  params$conversionpaf_link <- tmp_conversionpaf

  return(params)
}
