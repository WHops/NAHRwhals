

ref_fa = '/Users/hoeps/PhD/projects/huminvs/genomes/hg38/hg38.fa'
asm_fa = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/alns/NA12878_giab_pbsq2-ccs_1000-hifiasm.h1-un.fasta'

allpaf = '/Users/hoeps/PhD/projects/nahrcall/revisions/WGS_mode/breakpoint_lab/new/bps/data/wga_primary.paf'
out_dir = 'wga'
bedtools_bin = 'bedtools'
minimap2_bin = '/Users/hoeps/opt/anaconda3/envs/snakemake/envs/nahrwhalsAPR/bin/minimap2'

align_all_vs_all_using_minimap2(minimap2_bin, ref_fa, asm_fa, allpaf)

make_genome_file(in_fa_fai, genome)
extract_test_list_from_paf(allpaf, out_dir, genome, bedtools_bin)



# Parallel tests

f <- function(i) {
  lmer(Petal.Width ~ . - Species + (1 | Species), data = iris)
}

system.time(save1 <- lapply(1:100, f))
system.time(save2 <- mclapply(1:100, f))
