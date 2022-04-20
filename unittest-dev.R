# Whoeps, 11th Apr 2022
# Sort-of, custom unittest


library(devtools)

hg38_fa = "/Users/hoeps/PhD/projects/huminvs/genomes/hg38/hg38.fa"
aln_dir = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/alns/'
conersionpaf_dir = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/liftover_custom/'

test_list = 'unittest-dev.txt'

assembly_fastas = c(
  'HG00512_h2' = paste0(aln_dir, 'HG00512_hgsvc_pbsq2-ccs_1000-hifiasm.h2-un.fasta'),
  'HG00731_h1' = paste0(aln_dir, 'HG00731_hgsvc_pbsq2-ccs_1000-hifiasm.h1-un.fasta'),
  'HG00731_h2' = paste0(aln_dir, 'HG00731_hgsvc_pbsq2-ccs_1000-hifiasm.h2-un.fasta'),
  'HG00733_h1' = paste0(aln_dir, 'HG00733_hgsvc_pbsq2-ccs_1000-hifiasm.h1-un.fasta'),
  'NA12878_h1' = paste0(aln_dir, 'NA12878_giab_pbsq2-ccs_1000-hifiasm.h1-un.fasta'),
  'NA12878_h2' = paste0(aln_dir, 'NA12878_giab_pbsq2-ccs_1000-hifiasm.h2-un.fasta'),
  'NA19240_h2' = paste0(aln_dir, 'NA19240_hgsvc_pbsq2-ccs_1000-hifiasm.h2-un.fasta')
)

conversion_pafs = c(
  'HG00512_h2' = paste0(conersionpaf_dir, 'HG00512_hgsvc_pbsq2-ccs_1000-hifiasm.h2-un_hg38.paf'),
  'HG00731_h1' = paste0(conersionpaf_dir, 'HG00731_hgsvc_pbsq2-ccs_1000-hifiasm.h1-un_hg38.paf'),
  'HG00731_h2' = paste0(conersionpaf_dir, 'HG00731_hgsvc_pbsq2-ccs_1000-hifiasm.h2-un_hg38.paf'),
  'HG00733_h1' = paste0(conersionpaf_dir, 'HG00733_hgsvc_pbsq2-ccs_1000-hifiasm.h1-un_hg38.paf'),
  'NA12878_h1' = paste0(conersionpaf_dir, 'NA12878_giab_pbsq2-ccs_1000-hifiasm.h1-un_hg38.paf'),
  'NA12878_h2' = paste0(conersionpaf_dir, 'NA12878_giab_pbsq2-ccs_1000-hifiasm.h2-un_hg38.paf'),
  'NA19240_h2' = paste0(conersionpaf_dir, 'NA19240_hgsvc_pbsq2-ccs_1000-hifiasm.h2-un_hg38.paf')
)


determine_xpad <- function(start, end){
  # Play with padding values
  if ((end - start) < 1000){
    xpad = 3
  } else if ((end - start) < 100000){
    xpad = 3
  } else {
    xpad = 2
  }
  return(xpad)
}

determine_chunklen_compression <- function(start, end){
  if ((end - start) > 500 * 1000){
    chunklen = 10000
    compression = 1000
  } else if (((end - start)) < 50 * 1000){
    chunklen = 1000
    compression = 100 # min(100, ((end-start)/10))
  } else { 
    chunklen = 1000
    compression = 100
  }
  
  return(c(chunklen, compression))
}

devtools::load_all()
tests = read.table(test_list, header=T, sep='\t')

for (nrow in 1:10){#dim(tests)[1]){
  nrun = 1
  interval = tests[nrow,]
  chunk_comp = determine_chunklen_compression(interval$start, interval$end)
  chunklen = chunk_comp[1]
  compression = chunk_comp[2]
  #compression = 10§§
    wrapper_aln_and_analyse(interval$chr,
                            interval$start,
                            interval$end,
                            hg38_fa,
                            assembly_fastas[[interval$sample]],
                            conversion_pafs[[interval$sample]],
                            samplename = paste0(interval$sample, '-', interval$category, '-', nrun),
                            chunklen = chunklen,
                            sd_minlen = compression,
                            compression = compression,
                            depth = 2, 
                            xpad = determine_xpad(interval$start, interval$end),
                            logfile = 'res/unittest_v3.tsv',
                            debug=F)
}
  

