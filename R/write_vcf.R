# @title Preprocess input data
# @description Remove alt samples, replace mut_max by 'UNK' where res_max < threshold, extract sample_id and phase from sample
# @param input_file A character string specifying the path to the input file
# @param res_max_threshold A numeric value specifying the res_max threshold for assigning 'UNK' to mut_max
# @author Wolfram Hoeps
# @export
# @return A preprocessed data frame
preprocess_data <- function(input_file, res_max_threshold = 0.98){
  df <- read.delim(input_file)
  df <- df[!grepl("\\.alt$", df$sample),]
  df[df$res_max < res_max_threshold,]$mut_max <- 'UNK'
  df <- df %>%
    dplyr::mutate(sample_id = sub("(\\.|_)([^0-9]).*$", "", df$sample), 
                  phase = ifelse(grepl("h[ap]{0,2}2", df$sample), 2, 1)) 
  return(df)
}

# @title Generate genotype strings
# @description Generate the genotype strings for each unique start and mutation
# @param df A data frame from the preprocess_data function
# @param sample_list A character vector of unique sample ids
# @param haplotype_list A numeric vector specifying the haplotype options
# @author Wolfram Hoeps
# @export
# @return A data frame with all the lines for the VCF file
generate_gt_strings <- function(df, sample_list, haplotype_list){
  unique_starts <- unique(df$start)
  all_lines <- data.frame()
  
  for (start in unique_starts){
    df_subset <- df[df$start == start,]
    df_subset$mut_maxsimple <- gsub("[0-9_]+", "", df_subset$mut_max)
    
    for (mutation in unique(df_subset$mut_maxsimple)){
      
      if (mutation %in% c('ref', 'UNK', '.')){
        next()
      }
      
      df_temp <- df_subset
      df_temp[df_temp$mut_maxsimple == mutation, 'GT'] <- 1
      df_temp[df_temp$mut_maxsimple != mutation, 'GT'] <- 0
      df_temp[df_temp$mut_maxsimple == 'UNK', 'GT'] <- '.'
      
      all_combinations <- expand.grid(sample_id = sample_list, phase = haplotype_list)
      df_complete <- merge(all_combinations, df_temp, by = c("sample_id", "phase"), all.x = TRUE)
      df_complete$GT[is.na(df_complete$GT)] <- '.'
      
      GT_combined <- dplyr::summarise(
        dplyr::group_by(df_complete, sample_id),
        GT_combined = paste0(GT, collapse = "|")
      )
      
      GTS = data.frame(t(GT_combined$GT_combined), stringsAsFactors = FALSE)
      colnames(GTS)=GT_combined$sample_id
      
      finished_line <- create_vcf_line(df_temp, GTS)
      all_lines <- rbind(all_lines, finished_line)
    }
  }
  
  return(all_lines)
}

# @title Create each line of the VCF
# @description Create each line of the VCF with the genotype string
# @param df A data frame from the generate_gt_strings function
# @param GTS A data frame with the genotype strings
# @author Wolfram Hoeps
# @export
# @return A data frame with a single VCF line
create_vcf_line <- function(df, GTS){
  df_line_preGT <- data.frame(
    CHROM = df[1,'seqname'],
    POS = df[1,'start'],
    ID = paste0('<', 'NWhal.', df[1,'mut_maxsimple'], '_', df[1,'seqname'], 
                '-', df[1,'start'], '-', df[1,'end'], '>'),
    REF = 'N',
    ALT = paste0('<', df[1,'mut_maxsimple'], '>'),
    QUAL = '.',
    FILTER = 'PASS',
    INFO = paste0('SVDESC=', df[1,'mut_maxsimple'], ';', 
                  'SVDEPTH=', str_count(df[1,'mut_maxsimple'], '\\+')+1 , ';',
                  'SVLEN=', df[1,'width_orig'] , ';', 
                  'END=', df[1,'start']),
    FMT = 'GT'
  )
  
  finished_line <- cbind(df_line_preGT, GTS)
  return(finished_line)
}

# @title Generate the header
# @description Generate the header for the VCF file
# @param all_lines A data frame from the generate_gt_strings function
# @param GTS A data frame with the genotype strings
# @author Wolfram Hoeps
# @export
# @return A character vector with the VCF header
generate_header <- function(all_lines, GTS){
  contigs_header <- paste0("##contig=<ID=", unique(all_lines$CHROM), ">")
  header <- c("##fileformat=VCFv4.2",
              contigs_header,
              "##INFO=<ID=SVDESC,Number=1,Type=String,Description=\"SV description\">",
              "##INFO=<ID=SVDEPTH,Number=1,Type=Integer,Description=\"SV depth, i.e number of consecutive SVs the variant represents\">",
              "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"SV length\">",
              "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">",
              "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
              paste0("#", paste(c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", 
                                  colnames(GTS)), collapse="\t"))
  )
  return(header)
}

# @title Write VCF file
# @description Create VCF file from input data and save it
# @param input_file A character string specifying the path to the input file
# @param output_file A character string specifying the path to the output VCF file
# @param res_max_threshold A numeric value specifying the res_max threshold for assigning 'UNK' to mut_max
# @author Wolfram Hoeps
# @export
write_vcf <- function(input_file, output_file, res_max_threshold = 0.98, sort = F){
  
  # Preprocess the data
  df <- preprocess_data(input_file, res_max_threshold)
  
  # Generate unique lists for samples and haplotypes
  sample_list <- unique(df$sample_id)
  haplotype_list <- c(1,2)
  
  # Generate genotype strings
  all_lines <- generate_gt_strings(df, sample_list, haplotype_list)
  
  if (sort == T){
    all_lines = all_lines[order(all_lines$CHROM, all_lines$POS),]
  }
  
  # Generate the VCF header
  header <- generate_header(all_lines, all_lines)
  
  # Convert data to character and generate VCF body
  body <- apply(all_lines, 1, function(x) paste(x, collapse="\t"))
  
  # Combine header and body
  vcf_data <- c(header, body)
  
  # Write to VCF file
  writeLines(vcf_data, output_file)
}

