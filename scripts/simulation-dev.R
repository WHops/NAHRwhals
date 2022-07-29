# Whoeps, 11th Apr 2022
# Sort-of, custom unittest

make_random_bed_2sd_pairs <- function(seqsize, sdsize, border, similarity){
  
  
  # make seed bed
  system(paste0('echo "simcontig\t',0, '\t',format(sdsize, scientific=F),'\n',
                'simcontig\t',0, '\t',format(sdsize, scientific=F),'\n',
                'simcontig\t',0, '\t',format(sdsize, scientific=F),'\n',
                'simcontig\t',0, '\t',format(sdsize, scientific=F),
                '" > seed.bed'))

  # genome thingy
  system(paste0('echo "simcontig\t',format(seqsize, scientific=F),'" > genome.genome'))
  
  # shuffle
  print('shuffling')
  system(paste0('bedtools shuffle -i seed.bed -g genome.genome -noOverlapping -maxTries 1000000000000 > shuffled.bed'))
  print('shuffling ended')
  
  # determine orientations
  orientations = sample(c("+", "-"), 2, replace=T)
  # Who belongs together?
  pairs = sample(c(1,2,3,4), 4, replace=F)

  # Load here
  bedpos = read.table('shuffled.bed')
  bedpos$V2 = bedpos$V2 + border
  sd1 = min(bedpos[c(pairs[1],pairs[2]),'V2'])
  sd2 = max(bedpos[c(pairs[1],pairs[2]),'V2'])
  sd3 = min(bedpos[c(pairs[3],pairs[4]),'V2'])
  sd4 = max(bedpos[c(pairs[3],pairs[4]),'V2']) 
  
  outbed = data.frame(v1='simcontig', v2=c(sd1, sd3), v3=c(sd1+sdsize, sd3+sdsize), v4=c('SDA','SDB'), 
                      V5='simcontig', V6=c(sd2, sd4), v7=c(sd2+sdsize, sd4+sdsize), v8=orientations, v9=similarity/100)
  
  write.table(outbed, file='outbed.bed', quote=F, row.names=F, col.names=F, sep='\t')
  
  print(outbed)
  print('outbed saved.')
}



library(devtools)
library(dplyr)
devtools::load_all()
#set.seed(1234)

# 
for (sd_sim in c(90, 95, 99)){
for (seqlen in c(1000, 5000, 10000)){#}, 20000, 100000, 1000000)){
#for (sd_sim in c(99)){
#  for (seqlen in c(100000)){
    for (i in 1:100){
      
#sdfile_bed_1 = system.file('extdata', 'sds_twonest.bed', package='nahrtoolkit')
#sdfile_bed_1000 = system.file('extdata', 'sds_twonest.bed', package='nahrtoolkit')

      # Prepare default run parameters
      params_nahrtoolkit_default = list(
        debug = F,
        plot_only = F,
        xpad = 1,
        chunklen = determine_chunklen_compression(0, seqlen),
        plot_minlen = determine_plot_minlen(0, seqlen),
        auto_find_compression = T,
        mode = 'precise',
        minlen = 0, 
        compression = 0,
        max_n_alns = 80,
        n_tests = 10,
        n_max_testchunks = 5,
        max_size_col_plus_rows = 100,
        baseline_log_minsize_min = log2(1),
        baseline_log_minsize_max = max(log2(20000), log2((seqlen) / 10)),
        depth = 3,
        clean_after_yourself = T
      )
outfasta = './simulated_seq_twonest.fa'

make_random_bed_2sd_pairs(seqlen, seqlen/10, border = seqlen/10, sd_sim)
sdfile_bed = 'outbed.bed'# paste0('/Users/hoeps/PhD/projects/nahrcall/nahrchainer/trash/12jul/§sd_beds/sds_twonest_eq/sds_twonest_eq_', format(seqlen / 1000, scientific=F), 'k_', sd_sim,'.bed')


sdfile_tsv = './sds_twonest.tsv'
# Turn the easy-to-write .bed file into a .tsv file
nahrtoolkit::convert_bed_to_tsv(sdfile_bed, sdfile_tsv)

## MAKE THE SEQUENCE
nahrtoolkit::simulate_seq(seqlen * 1.2, sdfile_tsv, outfasta)





depth = 2

for (sdpair in c('SDa', 'SDb')){
  for (svtype in c('inv','del','dup')){

    if (depth == 1){
      svfile = paste0('/Users/hoeps/PhD/projects/nahrcall/nahrchainer/trash/12jul/sv_actions/', sdpair, '_', svtype, '.txt')

      outmutfasta = './mutation.fa'
      outmutsd = './mutation.tsv'

      x = tryCatch(nahrtoolkit::mutate_seq(outfasta, sdfile_tsv, svfile, outmutfasta, outmutsd), error = function(e) e)

      if (!(any(class(x) == 'error'))){
        print('I like it')
          wrapper_aln_and_analyse('chr0',
                                  0,
                                  seqlen*1.2,
                                  outfasta,
                                  outmutfasta,
                                  NULL,
                                  samplename = paste0(sdpair, '-', svtype, '-', sd_sim),
                                  params = params_nahrtoolkit_default,
                                  logfile = paste0('res/chr0-0-',format((seqlen*1.2), scientific=F),'/calls.tsv'))
      } else {
        print('I dont like it')
      }

    } else if (depth == 2){

      svfile = paste0('/Users/hoeps/PhD/projects/nahrcall/nahrchainer/trash/12jul/sv_actions/', sdpair, '_', svtype, '.txt')

      outmutfasta = './mutation.fa'
      outmutsd = './mutation.tsv'

      x = tryCatch(mutate_seq(outfasta, sdfile_tsv, svfile, outmutfasta, outmutsd), error = function(e) e)

      
      if (!(any(class(x) == 'error'))){

        wrapper_aln_and_analyse('chr0',
                                0,
                                seqlen*1.2,
                                outfasta,
                                outmutfasta,
                                NULL,
                                samplename = paste0(sdpair, '-', svtype, '-depth1-', sd_sim),
                                params = params_nahrtoolkit_default,
                                logfile = paste0('res/chr0-0-',format((seqlen*1.2), scientific=F),'/calls.tsv'))
        
        if (depth == 2){
          for (sdpair2 in c('SDa', 'SDb')){
            for (svtype2 in c('inv','del','dup')){
  
              svfile2 = paste0('/Users/hoeps/PhD/projects/nahrcall/nahrchainer/trash/12jul/sv_actions/', sdpair2, '_', svtype2, '.txt')
  
              outmutfasta2 = './mutation2.fa'
              outmutsd2 = './mutation2.tsv'
  
              y = tryCatch(nahrtoolkit::mutate_seq(outmutfasta, outmutsd, svfile2, outmutfasta2, outmutsd2), error = function(e) e)
  
              if (!(any(class(y) == 'error'))){
                wrapper_aln_and_analyse('chr0',
                                        0,
                                        seqlen*1.2,
                                        outfasta,
                                        outmutfasta2,
                                        NULL,
                                        samplename = paste0(sdpair, '-', svtype, '+', sdpair2, '-', svtype2, '-', sd_sim, '-', i),
                                        params = params_nahrtoolkit_default,
                                        logfile = paste0('res/chr0-0-',format(seqlen*1.2, scientific=F),'/calls.tsv'))
          }
        }
          }
        }
  }
    }
  }
}
}
}
}
# Make an exact dotplot, in case this is feasible.
# if (seqlen <= 25000){
#   exact_plot = nahrtoolkit::make_dotplot(outfasta, outmutfasta2, 15, save=F)
#   print(exact_plot  + ggplot2::labs(title='Simulated Sequence #3. Two nested SD pairs in indirect orientation.'))
# }








             
