# Fasta extraction and asssembly tools

#' extract_subseq
#' helperfunction to extract a subfasta by coordinates
#' 
#' Is awfully slow unfortunately. We should find a way to only read
#' the chromosome of interest. 
#' @author Wolfram Hoeps
#' @rdname alignment
#' @export
extract_subseq <- function(infasta, seqname, start, end, outfasta){

  # Read the whole infasta
  seq = Biostrings::readDNAStringSet(infasta)
  
  # Subset it to the sequence name of interest
  if (!is.null(seqname)){
    subseq = Biostrings::subseq(seq[[seqname]], start=start, end=end)
  } else {
    subseq = Biostrings::subseq(seq, start=start, end=end)
  }
  # Write fasta. Imported function from seqbilder_functions.R
  writeFasta(data.frame(name = 'seq', seq = as.character(subseq)), outfasta)

}


#' @description A core function. Give me a fasta with a sequence of interest (e.g. from hg38),
#' and give me an assembly fasta, and I will return the fasta of a corresponding
#' sequence in the assembly. sas
#'
#' @param targetfasta [character/link] link to the 'target' single-sequence fasta (sometimes reference, e.g. chm13.)
#' @param queryfasta [character/link] link to the 'query' single-sequence fasta.
#' @return nothing. But output files written.
#'
#' @author Wolfram Höps
#' @rdname alignment
#' @export
liftover_segment_direct <- function(targetfasta, queryfasta){
  
  # Queryfasta: e.g. a piece of hg38 sequence.
  # Targetfasta: typically a de-novo alignment. 
  
  queryfasta = outfasta
  targetfasta = "/Users/hoeps/PhD/projects/huminvs/genomes/hifi-asm/HG00512_hgsvc_pbsq2-ccs_1000-hifiasm.h1-un.fasta"
  outpaf_tmp = "./tpm_out.paf"
  # Make an alignment with mmap2.
  # We align the sequence we want to find to the assembly (targetfasta)
  run_minimap2(queryfasta, targetfasta, outpaf_tmp)
  
  
}

# outpaf_tmp = '/Users/hoeps/phd/projects/nahrcall/nahrchainer/seqbuilder/res/oa1.chunk'
# outpaf_tmp2 = '/Users/hoeps/phd/projects/nahrcall/nahrchainer/seqbuilder/res/oa1.chunk.filter'
# 
# filter_paf_to_main_region(outpaf_tmp, outpaf_tmp2, inlen = 200000)

# To document
filter_paf_to_main_region <- function(paflink, outpaflink){
  
  paf = read.table(paflink, sep='\t', fill=T, row.names=NULL)
    
  colnames_paf = c('qname','qlen','qstart','qend',
                   'strand','tname','tlen','tstart',
                   'tend','nmatch','alen','mapq')
  colnames(paf)[1:length(colnames_paf)] = colnames_paf
  
  insequence_len = sum(na.omit(as.numeric(dplyr::slice(dplyr::group_by(paf,qname),1)$qlen)))
  
  paf_winners = dplyr::filter(dplyr::group_by(paf, qname), nmatch == max(nmatch))
  
  
  chr_winners = tail(names(sort(table(paf_winners$tname))), 1)
  middle_median = median(((paf_winners[paf_winners$tname == chr_winners,]$tend + 
                        paf_winners[paf_winners$tname == chr_winners,]$tstart)/2))
  
  start_winners = middle_median - (1.5*(insequence_len/2))
  end_winners = middle_median + (1.5*(insequence_len/2))
  
  paf_cut = paf[((paf$tname == chr_winners) &
                 (paf$tstart >= start_winners) &
                 (paf$tend <= end_winners)),]

  # Transform for complete madness.
  paf_cut = transform(paf_cut, 
                      qname = tname, 
                      tname = qname, 
                      qlen = tlen, 
                      tlen = qlen, 
                      qstart = tstart,
                      tstart = qstart,
                      qend = tend,
                      tend = qend)
  
  
  write.table(paf_cut, file=outpaflink, sep='\t', row.names=F, col.names=F, quote=F)

  
}



# hg38fa = "/Users/hoeps/PhD/projects/huminvs/genomes/hg38/hg38.fa"
# outfasta = "~/Desktop/trash.fa"
# seqname = 'chr11'
# start = 13104252
# end = 13122521
# 
# extract_subseq(hg38fa, seqname, start, end, outfasta)
# 


# TRASH #

# extract_subseq_bedtools <- function(infasta, seqname, start, end, outfasta){
#   
#   cmd = (paste0("bedtools getfasta -fi ", 
#                 infasta, 
#                 " -bed ", 
#                 "<(echo \"", seqname, "\t", start, "\t", end, "\")", 
#                 " > ", outfasta))
#   
#   system(paste0("echo ", cmd, ' >', ' ./tmp.sh'))
#   
#}

# infasta = "/Users/hoeps/PhD/projects/huminvs/genomes/hg38/hg38.fa"
# outfasta = "~/Desktop/trash.fa"
# seqname = 'chr11'
# start = 13104252
# end = 13122521
# 
# extract_subseq(infasta, seqname, start, end, outfasta)

