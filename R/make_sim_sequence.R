# Make sample sequence with 10 complex SVs. 

# Simulate a sequence with SDs

# Non-export!
make_sim_seq_wrapper <- function(){
seqlen = 140000
sdfile_bed = "/Users/hoeps/PhD/projects/nahrcall/nahrchainer/test/test_100k/orig_fa/sds_100k.bed"
simfasta = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/test/test_100k/orig_fa/sim100k.fa'

sdfile_tsv = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/test/test_100k/orig_fa/sds_100k.tsv'
# Turn the easy-to-write .bed file into a .tsv file 
convert_bed_to_tsv(sdfile_bed, sdfile_tsv)

## MAKE THE SEQUENCE
simulate_seq(seqlen, sdfile_tsv, simfasta)

# Introduce mutations

for (n in as.character(1:11)){
  print(n)
  svfile = paste0( '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/test/test_100k/svs/sv',n, '.tsv')
  outmutfasta = paste0( '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/test/test_100k/mut_fa/mut',n, '.fa')
  outmutsd = paste0( '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/test/test_100k/mut_fa/mut',n, '.tsv')
  mutate_seq(simfasta, sdfile_tsv, svfile, outmutfasta, outmutsd)
}
}
