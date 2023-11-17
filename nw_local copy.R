


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
  seqname_x = tests[row, 'chr']
  start_x =   as.numeric(tests[row, 'start'])
  end_x =     as.numeric(tests[row, 'end'])
  #print(asm_fa)
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
            debug=T,
            logfile = paste0('res/res_', strsplit(asm_fa, '/')[[1]][length((strsplit(asm_fa, '/')[[1]]))], '.tsv')
            
  )
  
}




#run_nw_once_specified(1, test_list, ref_fa, asm_fa, anntrack, minimap2_bin, 'Y', 'X')

# Hardcoded values moved out of the function
hg38_fa = "/Users/hoeps/PhD/projects/huminvs/genomes/hg38/centro_lab/hg38_chr1_22_XY.fa"
t2t_fa = "/Users/hoeps/PhD/projects/huminvs/genomes/T2T-CHM13v2.0/T2T-CHM13v2.0.fa"
ref_fa = t2t_fa
#ref_fa = "/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/tomato/download?path=S.lycopersicum.Heinz1706.genomic.fa"
#asm_fa = "/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/tomato/download?path=S.lycopersicum.M82.genomic.fa"
#ref_fa = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/athaliana/Col-0.fasta'
#asm_fa = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/athaliana/C24.chr.all.v2.0.fasta'
#asm_fa = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/alns/NA12878_giab_pbsq2-ccs_1000-hifiasm.h1-un.fasta'
anntrack = F#/Users/hoeps/PhD/projects/huminvs/analyses_paper/data/genes/hg38/for_ntk/genes_hg38.bed"
minimap2_bin = '/Users/hoeps/opt/anaconda3/envs/snakemake/bin/minimap2'
samplename_y = 'HG00733_h1' 
samplename_x = 'T2T'
asm_fa = assembly_fastas[samplename_y]

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

ara_fas = c(
  'An-1' = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/athaliana/An-1.chr.all.v2.0.fasta',
  'C24' = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/athaliana/C24.chr.all.v2.0.fasta',
  'Cvi' = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/athaliana/Cvi.chr.all.v2.0.fasta',
  'Eri' = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/athaliana/Eri.chr.all.v2.0.fasta',
  'Kyo' = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/athaliana/Kyo.chr.all.v2.0.fasta',
  'Ler' = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/athaliana/Ler.chr.all.v2.0.fasta',
  'Sha' = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/athaliana/Sha.chr.all.v2.0.fasta'
)

ara_names = c('An-1', 'C24','Cvi','Eri','Kyo','Ler','Sha')

tomato_fas = c(
  
)

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


#for (sample in ara_names[5:length(ara_names)]){
  ref_fa = t2t_fa
  asm_fa = hg38_fa
  sample = 'hg38'
  test_list = wga_write_interval_list(ref_fa, asm_fa, paste0('wga_t2t_',sample), 1000000, 10000, t2t_mask, threads)
  tests = read.table(test_list, sep='\t')
  genome_file = paste0('wga_t2t_',sample, '/ref.genome')
  
  #test_list = 'wga_hg38_NA12878_h1/list_cut_final.bed'
  #extra_list = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper_2/hg38_t2t/yang_et_al/data/validated_chm13_sort.bed'
  extra_list = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper_2/hg38_t2t/yang_et_al/data/task_d/SDRs.bed'
  make_karyogram(test_list, genome_file, extra_list, specified_text = 'Stuff')
  
  #for (i in 22:nrow(tests)){
  #   samplename_y = 'NA12878_h1'
  #   print(i)
  #   run_nw_once(i, test_list, ref_fa, asm_fa, anntrack, minimap2_bin, samplename_y, samplename_x)
  #}
  # 
  # for (i in 10:20){
  #   run_nw_once(i, test_list, ref_fa, asm_fa, anntrack, minimap2_bin, 'NA12878_h2', 'T2T')
  # }
  samplename_x = 't2t'
  samplename_y = sample#
  devtools::load_all()
  
#   
  run_nw_once <- function(row, test_list, ref_fa, asm_fa, anntrack, minimap2_bin, samplename_y, samplename_x, conversionpaf_link){

    tests = read.table(test_list, sep='\t')
    colnames(tests) = c('chr', 'start', 'end')
    seqname_x = tests[row, 'chr']
    start_x =   as.numeric(tests[row, 'start'])
    end_x =     as.numeric(tests[row, 'end'])
    print('#######################')
    print(row)
    print('#######################')

    #print(asm_fa)
    try(nahrwhals(genome_x_fa = ref_fa,
               genome_y_fa = asm_fa,
               seqname_x = seqname_x,
               start_x = start_x,
               end_x =   end_x,
               samplename_y = samplename_y,
               samplename_x = samplename_x,
               compression = 1000,
               minlen = 1000,
               max_size_col_plus_rows = 500,
               max_n_alns = 10000,
               self_plots = F,#AN
               anntrack = anntrack,
               plot_only = F,
               minimap2_bin = minimap2_bin,
               noclutterplots = F,
               use_paf_library = T,
               debug=F,
               init_width = 10,
               conversionpaf_link = conversionpaf_link,
               logfile = paste0('res/res_', strsplit(asm_fa, '/')[[1]][length((strsplit(asm_fa, '/')[[1]]))], '.tsv')
    )
    )

  }



  #ntest = 180
  date()
  results <- mclapply(150:200, function(idx) run_nw_once(idx, test_list, ref_fa, asm_fa, anntrack, minimap2_bin, samplename_y, samplename_x, 'wga_t2t_hg38/fullaln.paf_10kbp_chunked_corrected.paf'), mc.cores = 10)
  date()
  #
#   startsample = 0
#   step = 10
#   n_samples = 20#nrow(tests)
#   n_batches = ceiling(n_samples / step)-1
#   for (n_batch in 0:n_batches){
#     start = startsample + (n_batch * step) + 1
#     end = startsample + (n_batch * step) + step
#     print('next')
#     print(start)
#     print(end)
#     results <- mclapply(start:end, function(idx) run_nw_once(idx, test_list, ref_fa, asm_fa, anntrack, minimap2_bin, samplename_y, samplename_x, 'wga_t2t_hg38/fullaln.paf_10kbp_chunked_corrected.paf'), mc.cores = threads)
#     print('hi back')
#   }
# 
#   
# 
# for (idx in 1:10){
#   print('#######################')
#   print(idx)
#   run_nw_once(idx, test_list, ref_fa, asm_fa, anntrack, minimap2_bin, samplename_y, samplename_x, 'wga_t2t_hg38/fullaln.paf_10kbp_chunked_corrected.paf')
# }
#   
  
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


