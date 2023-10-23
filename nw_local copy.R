


#!/bin/Rscript
library(argparse)
library(devtools)
library(parallel)

devtools::load_all()
aln_dir = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/alns/'

assembly_fastas = c(
  #'HG00512_h2' = paste0(aln_dir, 'gzip/HG00512_hgsvc_pbsq2-ccs_1000-hifiasm.h2-un.fasta.gz'),
  'HG00512_h2' = paste0(aln_dir, 'HG00512_hgsvc_pbsq2-ccs_1000-hifiasm.h2-un.fasta'),
  
  'HG00513_h1' = paste0(aln_dir, 'HG00513_hgsvc_pbsq2-ccs_1000-hifiasm.h1-un.fasta'),
  'HG02818_h2' = paste0(aln_dir, 'HG02818_hgsvc_pbsq2-ccs_1000-hifiasm.h2-un.fasta'),
  'HG00731_h1' = paste0(aln_dir, 'HG00731_hgsvc_pbsq2-ccs_1000-hifiasm.h1-un.fasta'),
  'HG00731_h2' = paste0(aln_dir, 'HG00731_hgsvc_pbsq2-ccs_1000-hifiasm.h2-un.fasta'),
  'HG00732_h1' = paste0(aln_dir, 'HG00732_hgsvc_pbsq2-ccs_1000-hifiasm.h1-un.fasta'),
  'HG00733_h1' = paste0(aln_dir, 'HG00733_hgsvc_pbsq2-ccs_1000-hifiasm.h1-un.fasta'),
  'NA12878_h1' = paste0(aln_dir, 'NA12878_giab_pbsq2-ccs_1000-hifiasm.h1-un.fasta'),
  'NA12878_h2' = paste0(aln_dir, 'NA12878_giab_pbsq2-ccs_1000-hifiasm.h2-un.fasta'),
  'NA19240_h2' = paste0(aln_dir, 'NA19240_hgsvc_pbsq2-ccs_1000-hifiasm.h2-un.fasta'),
  'HC02666' = '/Users/hoeps/PhD/projects/chrY/asms/HC02666.HIFIRW.ONTUL.na.chrY.fasta',
  'HG03248' = '/Users/hoeps/PhD/projects/chrY/asms/HG03248.HIFIRW.ONTUL.na.chrY.fasta',
  'HG02018_h1' = paste0(aln_dir, 'HG02018.asm.dip.hap1.p_ctg.fasta'),
  'HG02018_h2' = paste0(aln_dir, 'HG02018.asm.dip.hap2.p_ctg.fasta'),
  'GM19129_h1' = paste0(aln_dir, 'GM19129.asm.dip.hap1.p_ctg.fasta'),
  'HG02282_h1' = paste0(aln_dir, 'HG02282.asm.hic.hap1.p_ctg.fasta'),
  'HG02953_h1' = paste0(aln_dir, 'HG02953.asm.hic.hap1.p_ctg.fasta'),
  'HG03452_h1' = paste0(aln_dir, 'HG03452.asm.hic.hap1.p_ctg.fasta'),
  
  'panTro' = paste0(aln_dir, 'panTro6.fasta'),
  'ponAbe' = "/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/alns/ape/ponAbe3.fa",#paste0(aln_dir, 'ponAbe3.fasta'),
  
  
  'T2T' = paste0(aln_dir, 'T2T-CHM13v2.0.fasta'),
  'hg38' = NA
)


run_nw_once_specified <- function(row, test_list, ref_fa, asm_fa, anntrack, minimap2_bin, samplename_y, samplename_x){
  #tests = read.table(test_list, sep='\t')
  colnames(tests) = c('chr', 'start', 'end')
  seqname_x = 'chr1'#tests[row, 'chr']
  start_x =   862092#as.numeric(tests[row, 'start'])
  end_x =     1761081#as.numeric(tests[row, 'end'])
  print(asm_fa)
  nahrwhals(genome_x_fa = ref_fa,
                genome_y_fa = asm_fa,
                seqname_x = seqname_x,
                start_x = start_x,
                end_x =   end_x,
                samplename_y = samplename_y,
                samplename_x = samplename_x,
                self_plots = F,#AN
                anntrack = anntrack,
                plot_only = F,
                minimap2_bin = minimap2_bin,
                noclutterplots = T,
                logfile = paste0('res/res_', strsplit(asm_fa, '/')[[1]][length((strsplit(asm_fa, '/')[[1]]))], '.tsv')
  
  )
  
}

run_nw_once <- function(row, test_list, ref_fa, asm_fa, anntrack, minimap2_bin, samplename_y, samplename_x){
  tests = read.table(test_list, sep='\t')
  colnames(tests) = c('chr', 'start', 'end')
  seqname_x = tests[row, 'chr']
  start_x =   as.numeric(tests[row, 'start'])
  end_x =     as.numeric(tests[row, 'end'])
  print(asm_fa)
  try(nahrwhals(genome_x_fa = ref_fa,
                genome_y_fa = asm_fa,
                seqname_x = seqname_x,
                start_x = start_x,
                end_x =   end_x,
                samplename_y = samplename_y,
                samplename_x = samplename_x,
                self_plots = F,#AN
                anntrack = anntrack,
                plot_only = F,
                minimap2_bin = minimap2_bin,
                noclutterplots = F,
                logfile = paste0('res/res_', strsplit(asm_fa, '/')[[1]][length((strsplit(asm_fa, '/')[[1]]))], '.tsv')
  )
  )
  
}
#run_nw_once_specified(1, test_list, ref_fa, asm_fa, anntrack, minimap2_bin, 'NA12878_h2', 'T2T')

# Hardcoded values moved out of the function
hg38_fa = "/Users/hoeps/PhD/projects/huminvs/genomes/hg38/centro_lab/hg38_masked.fa"
t2t_fa = "/Users/hoeps/PhD/projects/huminvs/genomes/T2T-CHM13v2.0/T2T-CHM13v2.0.fa"
ref_fa = t2t_fa
#ref_fa = "/hg38_fa/hoeps/PhD/projects/nahrcall/nahrchainer/data/tomato/download?path=S.lycopersicum.Heinz1706.genomic.fa"
#asm_fa = "/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/tomato/download?path=S.lycopersicum.M82.genomic.fa"
#ref_fa = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/athaliana/Col-0.fasta'
#asm_fa = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/athaliana/C24.chr.all.v2.0.fasta'
asm_fa = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/alns/NA12878_giab_pbsq2-ccs_1000-hifiasm.h1-un.fasta'
anntrack = F#/Users/hoeps/PhD/projects/huminvs/analyses_paper/data/genes/hg38/for_ntk/genes_hg38.bed"
minimap2_bin = '/Users/hoeps/opt/anaconda3/envs/snakemake/bin/minimap2'
samplename_y = 'NA12878_h1' 
samplename_x = 'T2T'

threads = 5

# Inform the user what is going on, including some parameters etc. 
message("Running nahrwhals with the following parameters:")
message("ref_fa: ", ref_fa)
message("asm_fa: ", asm_fa)
message("anntrack: ", anntrack)
message("minimap2_bin: ", minimap2_bin)
message("samplename_y: ", samplename_y)
message("samplename_x: ", samplename_x)
message("threads: ", threads)
test_list = ''

asm_fas = c(
  '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/athaliana/An-1.chr.all.v2.0.fasta',
  '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/athaliana/C24.chr.all.v2.0.fasta',
  '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/athaliana/Cvi.chr.all.v2.0.fasta',
  '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/athaliana/Eri.chr.all.v2.0.fasta',
  '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/athaliana/Kyo.chr.all.v2.0.fasta',
  '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/athaliana/Ler.chr.all.v2.0.fasta',
  '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/athaliana/Sha.chr.all.v2.0.fasta'
)

asm_names = c('An-1', 'C24','Cvi','Eri','Kyo','Ler','Sha')

# for (i in 1:length(asm_names)){
#   browser()
#   asm_fa = asm_fas[i]
#   asm_name = asm_names[i]
#   test_list = wga_write_interval_list(ref_fa, asm_fa, paste0('wga_test_ara_',asm_name), 1000000, 100000, threads)
#   run_nw_once(i, test_list, ref_fa, asm_fa, anntrack, minimap2_bin, asm_names[i], samplename_x)
# }
# 
# ref_fa = asm_fas[1]
# asm_fa = asm_fas[2]
# test_list = wga_write_interval_list(ref_fa, asm_fa, paste0('wga_test_ara'), 1000000, 100000, threads)
# tests = read.table(test_list, sep='\t')
# inform_about_tests(test_list) 
# 
# # Inform the user how many tests we are running
# message("Running ", nrow(tests), " tests")

# for (i in 1:length(asm_names)){
#   asm_fa = asm_fas[i]
#   asm_name = asm_names[i]
t2t_mask = '~/PhD/projects/huminvs/genomes/T2T-CHM13v2.0/censat/censat_noct_100kb_merge_2.bed'


for (sample in c('NA12878_h2', 'HG00733_h1', 'HG02018_h1', 'HG02018_h2', 'hg38', 'NA19240_h2')){
  
  asm_fa = assembly_fastas[sample]
  test_list = wga_write_interval_list(ref_fa, asm_fa, paste0('wga_t2t_',sample), 1000000, 10000, t2t_mask, threads)
  tests = read.table(test_list, sep='\t')
  genome_file = paste0('wga_t2t_',sample, '/ref.genome')
  
  #test_list = 'wga_hg38_NA12878_h1/list_cut_final.bed'
  #make_karyogram(test_list, genome_file, specified_text = 'T2T_NA12878_h1')
  
  #for (i in 22:nrow(tests)){
  #   samplename_y = 'NA12878_h1'
  #   print(i)
  #   run_nw_once(i, test_list, ref_fa, asm_fa, anntrack, minimap2_bin, samplename_y, samplename_x)
  #}
  # 
  # for (i in 10:20){
  #   run_nw_once(i, test_list, ref_fa, asm_fa, anntrack, minimap2_bin, 'NA12878_h2', 'T2T')
  # }
  
  sample = samplename_y
  results <- mclapply(20:nrow(tests), function(idx) run_nw_once(idx, test_list, ref_fa, asm_fa, anntrack, minimap2_bin, samplename_y, samplename_x), mc.cores = threads)
}
# }
# 
# 
# 
# for (i in 1:length(asm_fas)){
#   print(asm_fas[1])
#   asm_fa = asm_fas[i]
#   run_nw_once(1, test_list, ref_fa, asm_fa, anntrack, minimap2_bin, asm_names[i], samplename_x)
#   
# }


