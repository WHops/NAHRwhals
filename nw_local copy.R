#!/bin/Rscript
library(argparse)
library(devtools)
library(parallel)

devtools::load_all()

run_nw_once <- function(row, test_list, ref_fa, asm_fa, anntrack, minimap2_bin, samplename_y, samplename_x){
  tests = read.table(test_list, sep='\t')
  colnames(tests) = c('chr', 'start', 'end')
  
  seqname_x = tests[row, 'chr']
  start_x =   tests[row, 'start']
  end_x =     tests[row, 'end']
  
  nahrwhals(genome_x_fa = ref_fa,
            genome_y_fa = asm_fa,
            seqname_x = seqname_x,
            start_x = start_x,
            end_x =   end_x,
            samplename_y = samplename_y,
            samplename_x = samplename_x,
            self_plots = T,
            anntrack = anntrack,
            plot_only = F,
            minimap2_bin = minimap2_bin
  )
}

# Hardcoded values moved out of the function
#hg38_fa = "/Users/hoeps/PhD/projects/huminvs/genomes/hg38/centro_lab/hg38_masked.fa"
##t2t_fa = "/Users/hoeps/PhD/projects/huminvs/genomes/T2T_v1.1/chm13.draft_v1.1.fasta"
#ref_fa = t2t_fa
ref_fa = "/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/tomato/download?path=S.lycopersicum.Heinz1706.genomic.fa"
asm_fa = "/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/tomato/download?path=S.lycopersicum.M82.genomic.fa"
#asm_fa = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/alns/NA12878_giab_pbsq2-ccs_1000-hifiasm.h1-un.fasta'
anntrack = "None"#/Users/hoeps/PhD/projects/huminvs/analyses_paper/data/genes/hg38/for_ntk/genes_hg38.bed"
minimap2_bin = '/Users/hoeps/opt/anaconda3/envs/snakemake/bin/minimap2'
samplename_y = 'M82'
samplename_x = 'Heinz1706'

threads = 8

# Inform the user what is going on, including some parameters etc. 
message("Running nahrwhals with the following parameters:")
message("ref_fa: ", ref_fa)
message("asm_fa: ", asm_fa)
message("anntrack: ", anntrack)
message("minimap2_bin: ", minimap2_bin)
message("samplename_y: ", samplename_y)
message("samplename_x: ", samplename_x)
message("threads: ", threads)
test_list = wga_write_interval_list(ref_fa, asm_fa, 'wga_tom_v1', threads)
tests = read.table(test_list, sep='\t')

# Inform the user how many tests we are running
message("Running ", nrow(tests), " tests")
results <- mclapply(1:nrow(tests), function(idx) run_nw_once(idx, test_list, ref_fa, asm_fa, anntrack, minimap2_bin, samplename_y, samplename_x), mc.cores = threads)

