# Whoeps, 11th Apr 2022
# Sort-of, custom unittest


library(devtools)
library(dplyr)
devtools::load_all()
seqlen = 25000
sdfile_bed = system.file('extdata', 'sds_twonest.bed', package='nahrtoolkit')

outfasta = './simulated_seq_twonest.fa'
sdfile_bed = 'trash/12jul/sds_twonest.bed'

sdfile_tsv = './sds_twonest.tsv'
# Turn the easy-to-write .bed file into a .tsv file
nahrtoolkit::convert_bed_to_tsv(sdfile_bed, sdfile_tsv)

## MAKE THE SEQUENCE
nahrtoolkit::simulate_seq(seqlen, sdfile_tsv, outfasta)
#> [1] "Done. Simulated sequence written to: ./simulated_seq_twonest.fa"

# Make an exact dotplot, in case this is feasible.
# if (seqlen <= 25000){
#   exact_plot = nahrtoolkit::make_dotplot(outfasta, outfasta, 15, save=F)
#   print(exact_plot  + ggplot2::labs(title='Simulated Sequence #3. Two nested SD pairs in indirect orientation.'))
# }

svfile_invA = system.file('extdata', 'SDa_inv.txt', package='nahrtoolkit')
svfile_delB = "/Users/hoeps/PhD/projects/nahrcall/nahrchainer/trash/12jul/SDb_del.txt"


outmutfasta = './simulated_seq_twonest_inv.fa'
outmutsd = './simulated_seq_twonest_inv.tsv'
outmutfasta2 = './simulated_seq_twonest_inv_del.fa'
outmutsd2 = './simulated_seq_twonest_inv_del.tsv'


nahrtoolkit::mutate_seq(outfasta, sdfile_tsv, svfile_invA, outmutfasta, outmutsd)
nahrtoolkit::mutate_seq(outmutfasta, outmutsd, svfile_delB, outmutfasta2, outmutsd2)

# Make an exact dotplot, in case this is feasible.
# if (seqlen <= 25000){
#   exact_plot = nahrtoolkit::make_dotplot(outfasta, outmutfasta2, 15, save=F)
#   print(exact_plot  + ggplot2::labs(title='Simulated Sequence #3. Two nested SD pairs in indirect orientation.'))
# }

params = list(
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
  n_tests = 20,
  n_max_testchunks = 5,
  max_size_col_plus_rows = 100,
  baseline_log_minsize_min = log2(1),
  baseline_log_minsize_max = max(log2(20000), log2((seqlen) / 10)),
  depth = 3,
  clean_after_yourself = T
)


wrapper_aln_and_analyse('chr0',
                        0,
                        seqlen,
                        outfasta,
                        outmutfasta2,
                        NULL,
                        samplename = 'test',#paste0(interval$sample, '-', interval$category, '-', nrun),
                        params = params,
                        logfile = 'res2/unittest_v4.tsv')



