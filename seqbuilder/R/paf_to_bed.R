# Turn paf to SD bed



# Go from a paf file (fresh out of minimap2) to a 
# filtered paf file (only nondiagonal) and then
# to a bedfile. The Bedfile will have the wrong colnames
# but it doesn't matter because it's saved without them. 
paf_write_bed <- function(inpaf_link, outsdbed_link){
  
  paf_sd = paf_to_sd_paf(inpaf_link)
  bed_sd = paf_to_bed(paf_sd)
  
  if (is.null(outsdbed_link)){
    return(bed_sd)
  } else if (!is.null(outsdbed_link)){
    write.table(bed_sd, outsdbed_link, sep='\t', col.names=F, row.names=F, quote=F)
  }
}


# Filter paf to sd relevant entries
paf_to_sd_paf <- function(inpaf_link){
  
  inpaf = read.table(inpaf_link, sep='\t')
  colnames_paf = c('qname','qlen','qstart','qend',
                   'strand','tname','tlen','tstart',
                   'tend','nmatch','alen','mapq')
  
  colnames(inpaf)[1:length(colnames_paf)] = colnames_paf
  
  # For safety: sort entries by qstart
  inpaf = inpaf[order(inpaf$qstart),]
  #paf_sd = inpaf
  paf_sd = inpaf[(
    (inpaf$qstart != inpaf$tstart) &
    (inpaf$qend != inpaf$tend) &
    (inpaf$qstart > inpaf$tstart)
      ),]
  
  return(paf_sd)

}


# Convert to bed
paf_to_bed <- function(inpaf){
  
  # In case colnames are not yet there
  colnames_paf = c('qname','qlen','qstart','qend',
                   'strand','tname','tlen','tstart',
                   'tend','nmatch','alen','mapq')
  
  colnames(inpaf)[1:length(colnames_paf)] = colnames_paf
  
  
  inpaf$id = inpaf$nmatch / inpaf$alen
  bed_sd = inpaf[,c('qname','qstart','qend','qname','qname', 'tstart','tend','strand','id')]
  bed_sd$qname = sub("_.*","", inpaf$qname)
  bed_sd$qname.2 = sub("_.*","", inpaf$qname)
  
  return(bed_sd)
}
