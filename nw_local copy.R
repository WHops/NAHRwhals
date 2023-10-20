#!/bin/Rscript
library(argparse)
library(devtools)
library(parallel)

devtools::load_all()

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
                logfile = paste0('res/res_', strsplit(asm_fa, '/')[[1]][length((strsplit(asm_fa, '/')[[1]]))], '.tsv')
  )
  )
  
}

i = 2
as = asm_fas[i]
#run_nw_once(1, test_list, ref_fa, as, anntrack, minimap2_bin, asm_names[i], 'An-1')


# Hardcoded values moved out of the function
hg38_fa = "/Users/hoeps/PhD/projects/huminvs/genomes/hg38/centro_lab/hg38_masked.fa"
t2t_fa = "/Users/hoeps/PhD/projects/huminvs/genomes/T2T_v1.1/chm13.draft_v1.1.fasta"
ref_fa = hg38_fa
#ref_fa = "/hg38_fa/hoeps/PhD/projects/nahrcall/nahrchainer/data/tomato/download?path=S.lycopersicum.Heinz1706.genomic.fa"
#asm_fa = "/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/tomato/download?path=S.lycopersicum.M82.genomic.fa"
#ref_fa = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/athaliana/Col-0.fasta'
#asm_fa = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/athaliana/C24.chr.all.v2.0.fasta'
asm_fa = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/alns/NA12878_giab_pbsq2-ccs_1000-hifiasm.h1-un.fasta'
anntrack = F#/Users/hoeps/PhD/projects/huminvs/analyses_paper/data/genes/hg38/for_ntk/genes_hg38.bed"
minimap2_bin = '/Users/hoeps/opt/anaconda3/envs/snakemake/bin/minimap2'
samplename_y = 'NA12878_h1' 
samplename_x = 'hg38'

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
#test_list = wga_write_interval_list(ref_fa, asm_fa, paste0('wga_hg38_NA12878_h1'), 1000000, 10000, threads)
tests = read.table(test_list, sep='\t')
#genome_file = 'wga_test_ara_An-1//ref.genome'

test_list = 'wga_hg38_NA12878_h1/list_cut_final.bed'
#make_karyogram(test_list, genome_file, specified_text = 'CP0')

for (i in 22:nrow(tests)){
  samplename_y = 'NA12878_h1'
  print(i)
  run_nw_once(i, test_list, ref_fa, asm_fa, anntrack, minimap2_bin, samplename_y, samplename_x)
}
#results <- mclapply(1:nrow(tests), function(idx) run_nw_once(idx, test_list, ref_fa, asm_fa, anntrack, minimap2_bin, samplename_y, samplename_x), mc.cores = threads)
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


