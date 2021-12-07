# Turn paf to SD bed




paf_write_bed <- function(inpaf_link, outsdbed_link){
  
paf_sd = paf_to_sd_paf(inpaf_link)
bed_sd = paf_to_bed(paf_sd)

write.table(bed_sd, outsdbed_link, sep='\t', col.names=F, row.names=F, quote=F)
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
  
  paf_sd = inpaf[(
    (inpaf$qstart != inpaf$tstart) & 
    (inpaf$qend != inpaf$tend) &
    (inpaf$qstart < inpaf$tstart)
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
