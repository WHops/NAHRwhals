#!/usr/local/bin/Rscript


#' A helperfunction (2/2) for the process of turning a bed to tsv. 
#'   
#' @param sds [data frame] an intermediate sds file between bed and tsv. Produced by 
#' make_sd_partners, the other helperfunction.
#' @return a tsv dataframe
#' 
#' @author Wolfram Höps
#' @rdname format_bed_to_tsv
#' @export
add_and_fix_columns <- function(sds){
  colnames_final = c(
    "chrom", "chromStart", "chromEnd", "name",
    "score","strand", "otherChrom","otherStart", "otherEnd",
    "otherSize", "uid", "posBasesHit", "testResult",
    "verdict", "chits","ccov","alignfile","alignL","indelN",
    "indelS","alignB","matchB","mismatchB", "transitionsB",
    "transversionsB","fracMatch","fracMatchIndel","jcK", "k2K"
  )
  
  # To fill with NA
  cols_fill_na = c('testResult','verdict','chits','ccov','alignfile')
  sds[,cols_fill_na] = 'N/A'
  
  # To fill with 0
  cols_fill_0 = c('score','otherSize','posBasesHit','indelN','indelS','transitionsB',
                  'transversionsB','jcK','k2K')
  sds[,cols_fill_0] = 0
  
  # To calculate
  sds$alignL = sds$chromEnd - sds$chromStart
  sds$alignB = sds$chromEnd - sds$chromStart
  sds$matchB = (sds$chromEnd - sds$chromStart) * sds$fracMatch
  sds$mismatchB = (sds$chromEnd - sds$chromStart) * (1-sds$fracMatch)
  sds$fracMatchIndel = sds$fracMatch
  
  # To fill custom
  sds$name = paste0(sds$otherChrom, ':', sds$chromStart)
  
  sds = sds[,colnames_final]
  return(sds)
}

#' A helperfunction (1/2) for the process of turning a bed to tsv. 
#'   
#' @param sds_raw [data frame] a loaded bedfile with column names annotated.
#' @return a dataframe that is in some intermediate state between bed and tsv and
#' which needs to be plugged into add_and_fix_columns (another helperfunction) next.
#' 
#' 
#' @author Wolfram Höps
#' @rdname format_bed_to_tsv
#' @export
make_sd_partners <- function(sds_raw){
  
  # Create the partners
  sds_partners <- transform(
    sds_raw, 
    
    chrom = otherChrom,
    otherChrom = chrom,
    
    chromStart = otherStart,
    otherStart = chromStart,
    
    chromEnd = otherEnd,
    otherEnd = chromEnd
  )
  
  sds = zipFastener(sds_raw, sds_partners, 1)
  
  return(sds)
}


#' A helperfunction for the helperfunction. Does some technical dataframe zipping stuff. 
#'   
#' @param df1 [data frame] a dataframe
#' @param df2 [data frame] a dataframe
#' @param along [1,2] axis along which to zip fasten merge whatever
#' @return merged dataframe
#' 
#' 
#' @author Unk - someone from the internet.
#' @rdname format_bed_to_tsv
#' @export
zipFastener <- function(df1, df2, along=2)
{
  # parameter checking
  if(!is.element(along, c(1,2))){
    stop("along must be 1 or 2 for rows and columns
                                              respectively")
  }
  # if merged by using zip feeding along the columns, the
  # same no. of rows is required and vice versa
  if(along==1 & (ncol(df1)!= ncol(df2))) {
    stop ("the no. of columns has to be equal to merge
               them by zip feeding")
  }
  if(along==2 & (nrow(df1)!= nrow(df2))) {
    stop ("the no. of rows has to be equal to merge them by
               zip feeding")
  }
  # zip fastener preperations
  d1 <- dim(df1)[along]
  d2 <- dim(df2)[along]
  i1 <- 1:d1           # index vector 1
  i2 <- 1:d2 + d1      # index vector 2
  # set biggest dimension dMax
  if(d1==d2) {
    dMax <- d1
  } else if (d1 > d2) {
    length(i2) <- length(i1)    # make vectors same length, 
    dMax <- d1                  # fill blanks with NAs   
  } else  if(d1 < d2){
    length(i1) <- length(i2)    # make vectors same length,
    dMax <- d2                  # fill blanks with NAs   
  }
  
  # zip fastener operations
  index <- as.vector(matrix(c(i1, i2), ncol=dMax, byrow=T))
  index <- index[!is.na(index)]         # remove NAs
  
  if(along==1){
    colnames(df2) <- colnames(df1)   # keep 1st colnames                  
    res <- rbind(df1,df2)[ index, ]  # reorder data frame
  }
  if(along==2) res <- cbind(df1,df2)[ , index]           
  return(res)
}



#' QC of the input bedfile. 
#'   
#' @param sds_raw [data frame] a loaded bedfile with column names annotated.
#' @return an error if there is something to complain about. (Päckla Schelln
#' is schnell aufgmacht :P)
#' 
#' @author Wolfram Höps
#' @rdname format_bed_to_tsv
#' @export
run_qc_bed_df <- function(sds_raw){
  stopifnot("Number of fields !=  9. Check fields (tab separation)" = 
              dim(sds_raw)[2] == 9)
  stopifnot("unique SD identifiers are not unique" = 
              dim(sds_raw)[1] == length(unique(sds_raw$uid)))
  stopifnot("Not all pairs are equally sized" = 
              unique((sds_raw$chromEnd - sds_raw$chromStart) == (sds_raw$otherEnd - sds_raw$otherStart)) == T)
  stopifnot("Partner #1: all starts should be smaller than ends" = 
              unique((sds_raw$chromEnd - sds_raw$chromStart) > 0) == T)
  stopifnot("Partner #2: all starts should be smaller than ends" = 
              unique((sds_raw$otherEnd - sds_raw$otherStart) > 0) == T)
}

#' Convert.
#' @rdname convert_bed_to_tsv
#' @export
convert_bed_to_tsv <- function(inbed, outtsv){

  #inbed = "/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/data/sds10y.bed"
  
  sds_raw = read.table(inbed, sep='\t')
  colnames(sds_raw) = c('chrom','chromStart',
                        'chromEnd', 'uid', 'otherChrom',
                        'otherStart','otherEnd', 'strand', 'fracMatch')
  
  # Run QC
  run_qc_bed_df(sds_raw)
  
  # Transform
  sds_partnered = make_sd_partners(sds_raw)
  sds_save = add_and_fix_columns(sds_partnered)
  
  # Save
  write.table(sds_save, file=outtsv, sep='\t', col.names = F, row.names = F, quote = F)
}

# runs only when script is run by itself
if (sys.nframe() == 0){
  
  # Define input
  inbed = commandArgs(trailingOnly=TRUE)[1]
  outtsv = commandArgs(trailingOnly=TRUE)[2]
  convert_bed_to_tsv(inbed, outtsv)
 
}




