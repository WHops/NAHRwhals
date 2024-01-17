# run NW

devtools::load_all()

hg38_fa = "/mnt/c/Users/Z364220/projects/data/genomes/hg38/hg38.fa"
t2t_fa = "/mnt/c/Users/Z364220/projects/data/genomes/t2t/chm13v2.0.fa"
HG002_mat_fa = "/mnt/c/Users/Z364220/projects/data/genomes/other_asms/HG002/maternal.fasta"
HG002_pat_fa = "/mnt/c/Users/Z364220/projects/data/genomes/other_asms/HG002/paternal.fasta"

sim_fa = "/mnt/c/Users/Z364220/projects/SD_mapping/simulation/res/seq_sim/single_pair_100kbp.fa"

hg19_chr15_fa = "/mnt/c/Users/Z364220/Downloads/chr15.fa"

DNA11_hap2 = "/mnt/c/Users/Z364220/projects/SD_mapping/data/DNA11-26362.bp.hap2.p_ctg.fa"


# Cycle genome_y_fa through the different assemblies of HG002


genome_x_fa = sim_fa
genome_y_fa = sim_fa
seqname_x = 'sim-sequence'
start_x = 100000
end_x = 900000
#start_x = 30590000
#end_x = 31500000
samplename_x = 'sim_seq'
samplename_y = 'sim_seq'

nahrwhals(
    genome_x_fa = genome_x_fa,
    genome_y_fa = genome_y_fa,
    seqname_x = seqname_x,
    start_x = start_x,
    end_x = end_x,
    samplename_x = samplename_x,
    samplename_y = samplename_y,
    plot_only=T,
    compare_full_fastas <- T
)
