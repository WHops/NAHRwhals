#### RUN THE FOLLOWING CODE AS A MINIMAL WORKING EXAMPLE:

# Load nahrtoolkit
library(devtools)
devtools::load_all()


# Take an example input fasta shipped with the package
samplefasta_link = system.file('extdata', '10ktest.fa', package='nahrtoolkit')

# Define an output alignment file
samplepaf_link = paste0(samplefasta_link, '.sample.paf')

# Make alignment & dotplot
make_chunked_minimap_alnment(samplefasta_link, samplefasta_link, samplepaf_link, 
                             chunklen = 10000, minsdlen = 100, saveplot=F, 
                             hllink = samplepaf_link, hltype = 'paf', quadrantsize = 10000)

# Make the condensed dotplot
grid = wrapper_paf_to_bitlocus(samplepaf_link, minlen=0, compression = 1)

# Run a mutation search
gridmatrix = gridlist_to_gridmatrix(grid[[3]])

res = explore_mutation_space(gridmatrix, depth = 3)


### [OVER] ###


# Take an example input fasta shipped with the package
samplefasta_link = system.file('extdata', 'simulated_seq_twonest.fa', package='nahrtoolkit')
samplefasta_inv_link = system.file('extdata', 'simulated_seq_twonest_inv.fa', package='nahrtoolkit')

# Define an output alignment file
samplepaf_link = paste0(samplefasta_link, '.sample2.paf')

# Make alignment & dotplot
make_chunked_minimap_alnment(samplefasta_link, samplefasta_inv_link, samplepaf_link, 
                             chunklen = 1000, minsdlen = 100, saveplot=F, 
                             hllink = samplepaf_link, hltype = 'paf', quadrantsize = 10000)

# Make the condensed dotplot
grid_inv = wrapper_paf_to_bitlocus(samplepaf_link, minlen=0, compression = 10)

# Run a mutation search
gridmatrix = gridlist_to_gridmatrix(grid_inv[[3]])
unmatching_bases(gridmatrix, verbose=T)

res = explore_mutation_space(gridmatrix, depth = 3)

# Define an output alignment file
samplepaf_link = paste0(samplefasta_link, '.sample_ref.paf')

# Make alignment & dotplot
make_chunked_minimap_alnment(samplefasta_link, samplefasta_link, samplepaf_link, 
                             chunklen = 1000, minsdlen = 100, saveplot=F, 
                             hllink = samplepaf_link, hltype = 'paf', quadrantsize = 10000)


# Make the condensed dotplot
grid = wrapper_paf_to_bitlocus(samplepaf_link, minlen=0, compression = 100)
# Run a mutation search
unmatching_bases(gridmatrix, verbose=T)

gridmatrix = gridlist_to_gridmatrix(grid[[3]])
res = explore_mutation_space(gridmatrix, depth = 3)





# Run a mutation search
res = explore_mutation_space(gridmatrix, depth = 3)






# ref = '/Users/hoeps/PhD/projects/huminvs/revision/alu_analysis/data/seqs/chr12-12391906-INV-1713/chr12-12391906-INV-1713_hg38_pm200kb.fa'
# alt = '/Users/hoeps/PhD/projects/huminvs/revision/alu_analysis/data/seqs/chr12-12391906-INV-1713/chr12-12391906-INV-1713_HG00512_h1_rev_rectified.fa'
# 
# # targetfa, queryfa
# exact_plot = make_dotplot(ref, alt, 7, save=F, 
#                           targetrange=c(199845-1000, 201913+1000), 
#                           queryrange=c(214844-1000 , 216914+1000))
# exact_plot + ggplot2::labs(y='HG00512_h1, corresponding region',
#                   title='Dotplot of chr12-12391906-INV-1713')
# 
# exact_plot = make_dotplot(ref, alt, 4, save=F, 
#                           targetrange=c(199845, 199845+400), 
#                           queryrange=c(214844 , 214844+400))
# exact_plot
# 
# targetfa = ref
# 
# make_dotplot(ref, alt, 8, save=F, 
#              targetrange=c(199845, 199845+500), 
#              queryrange =c(183086, 183086+500))
# # BP1
# alt_junk = get_subseq(alt, c(214844, 214844+400))
# ref_pre = get_subseq(ref, c(199845, 199845+400))
# ref_post = get_subseq(ref, c(201913-400,201913))
#   
# ref_post_rev = as.character(Biostrings::reverseComplement(Biostrings::DNAString(ref_post)))
# 
# alt_junk
# ref_pre
# ref_post_rev
# 
# # BP2
# alt_junk2 = get_subseq(alt, c(216914-400, 216914))
# ref_pre = get_subseq(ref, c(199845, 199845+400))
# 
# ref_post = get_subseq(ref, c(201913-400, 201913))
# ref_pre_rev = as.character(Biostrings::reverseComplement(Biostrings::DNAString(ref_pre)))
# 
# alt_junk2
# ref_pre_rev
# ref_post
# 

fa_link = './simulated_seq_10kb_del.fa'
chm13_link = '/Users/hoeps/PhD/projects/huminvs/genomes/CHM13_T2T/fasta/chm13.draft_v1.1.fasta'
outpaf_link = as.character(runif(1,1e10,1e11))

make_chunked_minimap_alnment(targetfasta=chm13_link, queryfasta=fa_link, outpaf_link, 
                             outplot=NULL, chunklen = 10000, minsdlen = 2000, saveplot=F, 
                             hllink = samplepaf_link, hltype = 'paf', quadrantsize = 100000)

grid = wrapper_paf_to_bitlocus(outpaf_link, gp = 1)


samplefasta_link = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/vignettes/simulated_seq_10kb_4SDs.fa'
samplemutfasta_link = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/vignettes/simulated_seq_10kb_del.fa'
samplepaf_link = paste0('blub2.paf')
make_chunked_minimap_alnment(samplefasta_link, samplemutfasta_link, samplepaf_link, 
                             chunklen = 2000, minsdlen = 0, saveplot=F, 
                             hllink = samplepaf_link, hltype = 'paf', quadrantsize = 1000)
grid = wrapper_paf_to_bitlocus(samplepaf_link, minlen = 10, compression = 1)
# Run a mutation search
gridmatrix = gridlist_to_gridmatrix(grid[[3]])
unmatching_bases(gridmatrix, verbose=T)

res = explore_mutation_space(gridmatrix, depth = 3)


library(nahrtoolkit)
samplefasta_link = system.file('extdata', '10ktest.fa', package='nahrtoolkit')



# blub 
samplefasta_link = system.file('extdata', '10ktest.fa', package='nahrtoolkit')
samplepaf_link = paste0(samplefasta_link, 'sample.paf')
make_chunked_minimap_alnment(samplefasta_link, samplefasta_link, samplepaf_link, 
                             chunklen = 10000, minsdlen = 100, saveplot=F, 
                             hllink = samplepaf_link, hltype = 'paf', quadrantsize = 10000)


grid = wrapper_paf_to_bitlocus(samplepaf_link, minlen=0, compression = 1)[[3]]
gridmatrix = gridlist_to_gridmatrix(grid)
unmatching_bases(gridmatrix, verbose=T)

# Run a mutation search
res = explore_mutation_space(gridmatrix, depth = 3)
res[res$eval == max(res$eval),]


### [OVER] ###

# samplefasta_link = '/Users/hoeps/phd/projects/nahrcall/nahrchainer/seqbuilder/data/realexample.fa'
# #samplefasta_link = '/Users/hoeps/PhD/projects/huminvs/genomes/hg38/wbs.fa'
# outpaf_link = '/Users/hoeps/phd/projects/nahrcall/nahrchainer/seqbuilder/res/x2'
# make_chunked_minimap_alnment('hg38', samplefasta_link, outpaf_link, 
#                                outplot=NULL, chunklen = 10000, minsdlen = 2000, saveplot=F, 
#                                hllink = outpaf_link, hltype = 'paf')
# vgrid = wrapper_paf_to_bitlocus(outpaf_link, minlen = 1000, gp = 1000)






samplefasta_link = '/Users/hoeps/phd/projects/nahrcall/nahrchainer/data/realexample.fa'
#samplefasta_link = '/Users/hoeps/PhD/projects/huminvs/genomes/hg38/wbs.fa'
outpaf_link = '/Users/hoeps/phd/projects/nahrcall/nahrchainer/res/x6223'
# make_chunked_minimap_alnment('/Users/hoeps/PhD/projects/huminvs/genomes/hifi-asm/HG00512_hgsvc_pbsq2-ccs_1000-hifiasm.h1-un.fasta', samplefasta_link, outpaf_link, 
#                              chunklen = 1000, minsdlen = 2000, saveplot=F, 
#                              hllink = F, hltype = F)
make_chunked_minimap_alnment(samplefasta_link, samplefasta_link, outpaf_link, 
                             chunklen = 1000, minsdlen = 1000, saveplot=F, 
                             hllink = F, hltype = F)


grid = wrapper_paf_to_bitlocus(outpaf_link, minlen=1000, compression = 1000)[[3]]
gridmatrix = gridlist_to_gridmatrix(grid)

# Run a mutation search
res = explore_mutation_space(gridmatrix, depth = 3)
res[res$eval == max(res$eval),]
  #'/Users/hoeps/PhD/projects/huminvs/genomes/hifi-asm/HG00512_hgsvc_pbsq2-ccs_1000-hifiasm.h1-un.fasta'
# Problems/Errors spotted: 

# 1) If outpaf exists, problems appear. SOMETIMES. 

# 2) Merging algorithm dies if there are too many alignments. 
samplefasta_link = '/Users/hoeps/PhD/projects/huminvs/genomes/hg38/wbs.fa'
outpaf_link = '/Users/hoeps/phd/projects/nahrcall/nahrchainer/res/outpaf9'
make_chunked_minimap_alnment(samplefasta_link, samplefasta_link, outpaf_link, 
                             chunklen = 10000, minsdlen = 500, saveplot=F, 
                             hllink = outpaf_link, hltype = 'paf', quadrantsize = 200000)
grid = wrapper_paf_to_bitlocus(outpaf_link, minlen = 10000, compression = 10000)
gridmatrix = gridlist_to_gridmatrix(grid[[3]])

res = explore_mutation_space(gridmatrix, depth = 2)
res[res$eval == max(res$eval),]

# 3) Alignments CAN break
samplefasta_link = '/Users/hoeps/PhD/projects/huminvs/genomes/hg38/s5.fa'
outpaf_link = '/Users/hoeps/phd/projects/nahrcall/nahrchainer/res/outpafsds1323354'
make_chunked_minimap_alnment(samplefasta_link, samplefasta_link, outpaf_link, 
                             chunklen = 1000, minsdlen = 10000, saveplot=F, 
                             hllink = outpaf_link, hltype = 'paf', quadrantsize = 200000)
grid = wrapper_paf_to_bitlocus(outpaf_link, minlen = 5000, compression = 5000)
gridmatrix = gridlist_to_gridmatrix(grid[[3]])

res = explore_mutation_space(gridmatrix, depth = 2)
res[res$eval == max(res$eval),]

### Stuff from Introduction
outfasta = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/vignettes/simulated_seq_twonest.fa'
outmutfasta = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/vignettes/simulated_seq_twonest_inv.fa'
outmutfasta2 = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/vignettes/simulated_seq_twonest_inv_dup.fa'
samplepaf_link = 'blub'
samplepaf_link2 = 'blub2'
samplepaf_link3 = 'blub3'

make_chunked_minimap_alnment(outfasta, outfasta, samplepaf_link, 
                             chunklen = 1000, minsdlen = 0, saveplot=F, 
                             hllink = samplepaf_link, hltype = 'paf', quadrantsize = 1000)
make_chunked_minimap_alnment(outfasta, outmutfasta, samplepaf_link2, 
                             chunklen = 1000, minsdlen = 0, saveplot=F, 
                             hllink = samplepaf_link, hltype = 'paf', quadrantsize = 1000)
make_chunked_minimap_alnment(outfasta, outmutfasta2, samplepaf_link3, 
                             chunklen = 1000, minsdlen = 0, saveplot=F, 
                             hllink = samplepaf_link, hltype = 'paf', quadrantsize = 1000)
grid = wrapper_paf_to_bitlocus(samplepaf_link, minlen=0, compression = 10)
grid = wrapper_paf_to_bitlocus(samplepaf_link2, minlen=0, compression = 10)
grid = wrapper_paf_to_bitlocus(samplepaf_link3, minlen=0, compression = 10)


# 4) Plot groundtruth of SDs (should work, but we need to asjust coordinate systems so that the SD
# Annotations from groundtruth start at zero too...) 


# Future: 

# a) implement a tree search
# b) implement a quality metric simple alignment test
# c) implement a 'region search': region in, and need to find correct region on chr. 
library(devtools)
devtools::load_all()
ref_fa = "/Users/hoeps/PhD/projects/huminvs/genomes/hg38/hg38.fa"
hg38_fa = ref_fa
hg38fa = hg38_fa
chm13fa = '/Users/hoeps/PhD/projects/huminvs/genomes/CHM13_T2T/fasta/chm13.draft_v1.1.fasta'
aln_fa = '/Users/hoeps/PhD/projects/huminvs/genomes/hifi-asm/HG00512_hgsvc_pbsq2-ccs_1000-hifiasm.h1-un.fasta'
conversionpaf_link = "/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/liftover_custom/HG00512_hgsvc_pbsq2-ccs_1000-hifiasm.h1-un_hg38.paf"

aln_fa = '~/Desktop/alns/NA12878_giab_pbsq2-ccs_1000-hifiasm.h2-un.fasta'
conversionpaf_link = "/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/liftover_custom/NA12878_giab_pbsq2-ccs_1000-hifiasm.h2-un_hg38.paf"

aln_fa = '~/Desktop/alns/HG00731_hgsvc_pbsq2-ccs_1000-hifiasm.h2-un.fasta'
conversionpaf_link = "/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/liftover_custom/HG00731_hgsvc_pbsq2-ccs_1000-hifiasm.h2-un_hg38.paf"


aln_fa = '~/Desktop/alns/NA19240_hgsvc_pbsq2-ccs_1000-hifiasm.h2-un.fasta'
conversionpaf_link = "/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/liftover_custom/NA19240_hgsvc_pbsq2-ccs_1000-hifiasm.h2-un_hg38.paf"


outfasta_hg38 = "/Users/hoeps/phd/projects/nahrcall/nahrchainer/res/hg38_sub.fa"
outfasta_aln = "/Users/hoeps/phd/projects/nahrcall/nahrchainer/res/aln_sub.fa"

outpaf_link = '/Users/hoeps/phd/projects/nahrcall/nahrchainer/res/122342526.paf'

  

rhemac10_fa = '/Users/hoeps/PhD/projects/huminvs/genomes/rheMac10/rheMac10.fa'
conversionpaf_link = '/Users/hoeps/PhD/projects/huminvs/genomes/rheMac10/rheMac10_hg38.paf'

wrapper_aln_and_analyse('chr1',
                        108348182 - 00000,
                        108546190 + 00000,
                        hg38_fa,
                        hg38_fa,
                        conversionpaf_link,
                        runname = 'hg38-rhesus-5',
                        sd_minlen = 100, 
                        compression = 100,
                        depth = 1)

seqname = 'chr7'
start = 142532924
end = 144507494








#hg38
seqname = 'chr1'
start = 5000000
end =   5500000


library(devtools)
devtools::load_all()

seqname = 'chr1'
start = 108348182
end = 108546191

# seqname = 'chrX'
# 
# 
# start = 136505547 
# end =   137131190 
outpaf_link = as.character(runif(1,1e10,1e11))

coords_liftover = liftover_coarse(seqname, start, end, conversionpaf_link, lenfactor = 1)
print(end - start)
print(coords_liftover$lift_end - coords_liftover$lift_start)
extract_subseq_bedtools(hg38fa, seqname, start, end, outfasta_hg38)
extract_subseq_bedtools(aln_fa, coords_liftover$lift_contig, coords_liftover$lift_start, coords_liftover$lift_end, outfasta_aln)

make_chunked_minimap_alnment(outfasta_hg38, outfasta_aln, outpaf_link,
                             chunklen = 1000, minsdlen = 2000, saveplot=F,
                             hllink = outpaf_link, hltype = 'paf')#, wholegenome = F)

make_chunked_minimap_alnment(outfasta_hg38, outfasta_hg38, outpaf_link,
                             chunklen = 1000, minsdlen = 2000, saveplot=F,
                             hllink = F, hltype = F, quadrantsize = 250000)#, wholegenome = F)
make_chunked_minimap_alnment(outfasta_aln, outfasta_aln, outpaf_link,
                             chunklen = 1000, minsdlen = 2000, saveplot=F,
                             hllink = F, hltype = F)#, wholegenome = F)


grid = wrapper_paf_to_bitlocus(outpaf_link, minlen = 1000, compression = 1000)#[[3]]
gridmatrix = gridlist_to_gridmatrix(grid)



res = explore_mutation_space(gridmatrix, depth = 1)

res[res$eval == max(res$eval),][1,]

gm2 = carry_out_compressed_sv(carry_out_compressed_sv(gridmatrix, c(3,33,'inv')), c(31,55,'del'))
plot_matrix(log(abs(gm2)+1)+1)
calc_coarse_grained_aln_score(gm2, verbose=T)

#gm1 = carry_out_compressed_sv(gridmatrix, c('3', '8', 'inv'))
#gm2 = carry_out_compressed_sv(gm1, c('9', '12', 'del'))
#plot_matrix(sign(gridmatrix))
#plot_matrix(sign(gm1))
#plot_matrix(sign(gm2))
# Run a mutation search
res = explore_mutation_space(gridmatrix, depth = 1)


profvis(explore_mutation_space(gridmatrix, depth = 2))
plot_matrix(gm2)
### [OVER] ###




outpaf_link = as.character(runif(1, 1e10,1e11))
make_chunked_minimap_alnment(outfasta_hg38, outfasta_hg38, outpaf_link,
                            chunklen = 1000, minsdlen = 2000, saveplot=F,
                            hllink = F, hltype = F)#, wholegenome = F)

make_dotplot(outfasta_hg38, outfasta_hg38, 15, save=F)

seqname = 'chr7'
start = 74869950
end = 75058098
coords_liftover = liftover_coarse(seqname, start, end, conversionpaf_link, factor = 5)
print(end - start)
print(coords_liftover$lift_end - coords_liftover$lift_start)
extract_subseq_bedtools(hg38fa, seqname, start, end, outfasta_hg38)
extract_subseq_bedtools(aln_fa, coords_liftover$lift_contig, coords_liftover$lift_start, coords_liftover$lift_end, outfasta_aln)
outpaf_link = 'deleteme.tmp'

make_chunked_minimap_alnment(outfasta_hg38, outfasta_aln, outpaf_link,
                             chunklen = 10000, minsdlen = 2000, saveplot=F,
                             hllink = F, hltype = F, wholegenome = F)

seqname = 'chr8'
start = 2023351
end = 2508385
coords_liftover = liftover_coarse(seqname, start, end, conversionpaf_link, lenfactor = 1)
print(end - start)
print(coords_liftover$lift_end - coords_liftover$lift_start)
extract_subseq_bedtools(hg38fa, seqname, start, end, outfasta_hg38)
extract_subseq_bedtools(aln_fa, coords_liftover$lift_contig, coords_liftover$lift_start, coords_liftover$lift_end, outfasta_aln)
outpaf_link = as.character(runif(1,1e10,1e11))

make_chunked_minimap_alnment(outfasta_hg38, outfasta_aln, outpaf_link,
                             chunklen = 10000, minsdlen = 2000, saveplot=F,
                             hllink = F, hltype = F, wholegenome = F)



seqname = 'chr19'
start = 54700000
end =   54900000
coords_liftover = liftover_coarse(seqname, start, end, conversionpaf_link, factor = 0)
print(end - start)
print(coords_liftover$lift_end - coords_liftover$lift_start)
extract_subseq_bedtools(hg38fa, seqname, start, end, outfasta_hg38)
extract_subseq_bedtools(aln_fa, coords_liftover$lift_contig, coords_liftover$lift_start, coords_liftover$lift_end, outfasta_aln)
outpaf_link = as.character(runif(1,1e10,1e11))

make_chunked_minimap_alnment(outfasta_hg38, outfasta_hg38, outpaf_link,
                             chunklen = 500, minsdlen = 2000, saveplot=F,
                             hllink = F, hltype = F, wholegenome = F)


seqname = 'chr19'
start = 53557693
end =   55538585
coords_liftover = liftover_coarse_attempt2(seqname, start, end, conversionpaf_link, lenfactor = 1)
print(end - start)
print(coords_liftover$lift_end - coords_liftover$lift_start)
extract_subseq_bedtools(hg38fa, seqname, start, end, outfasta_hg38)
extract_subseq_bedtools(aln_fa, coords_liftover$lift_contig, coords_liftover$lift_start, coords_liftover$lift_end, outfasta_aln)
outpaf_link = as.character(runif(1,1e10,1e11))

make_chunked_minimap_alnment(outfasta_hg38, outfasta_aln, outpaf_link,
                             chunklen = 1000, minsdlen = 2000, saveplot=F,
                             hllink = F, hltype = F)#, wholegenome = F)

outpaf_link = as.character(runif(1,1e10,1e11))
make_chunked_minimap_alnment(outfasta_hg38, outfasta_hg38, outpaf_link,
                             chunklen = 10000, minsdlen = 2000, saveplot=F,
                             hllink = F, hltype = F)



wrapper_dotplot_with_alignment(seqname, start, end, hg38fa, aln_fa, 
                               '~/Desktop/test2/hg38.fa', '~/Desktop/test2/aln.fa',
                               conversionpaf_link, '~/Desktop/test2/run.paf')

wrapper_dotplot_with_alignment <- function(seqname, start, end, genome_x_fa, genome_y_fa, subseqfasta_x, 
                                           subseqfasta_y, conversionpaf_link, outpaf_link, 
                                           chunklen = 1000, factor = 0.5){
  
  print('1')  
  # Get coords in assembly
  coords_liftover = liftover_coarse(seqname, start, end, conversionpaf_link, factor = 1)
  
  print('2')
  # Gimme fasta
  extract_subseq_bedtools(genome_x_fa, seqname, start, end, subseqfasta_x)
  extract_subseq_bedtools(genome_y_fa, coords_liftover$lift_contig, coords_liftover$lift_start, coords_liftover$lift_end, subseqfasta_y)
  
  print('3')
  #outpaf_link = as.character(runif(1,1e10, 1e11))
  plot = make_chunked_minimap_alnment(subseqfasta_y, subseqfasta_x, outpaf_link, 
                                      chunklen = 1000, minsdlen = 2000, saveplot=T, 
                                      hllink = F, hltype = F, wholegenome = T)
  return(plot)
  
}



  
HG00732sd1 = "/Users/hoeps/PhD/projects/huminvs/revision/recurrent_invs_bps/chr8_fine/fasta/SDs/HG00732_SD1.fa"
HG00732sd2 = "/Users/hoeps/PhD/projects/huminvs/revision/recurrent_invs_bps/chr8_fine/fasta/SDs/HG00732_SD2.fa"

HG00731sd1 = "/Users/hoeps/PhD/projects/huminvs/revision/recurrent_invs_bps/chr8_fine/fasta/SDs/HG00731_SD1.fa"
HG00731sd2 = "/Users/hoeps/PhD/projects/huminvs/revision/recurrent_invs_bps/chr8_fine/fasta/SDs/HG00731_SD2.fa"

afrh1sd1 = '/Users/hoeps/PhD/projects/huminvs/revision/recurrent_invs_bps/chr8_fine/fasta_afr/SDs/HG03125h1_SD1.fa'
afrh1sd2 = '/Users/hoeps/PhD/projects/huminvs/revision/recurrent_invs_bps/chr8_fine/fasta_afr/SDs/HG03125h1_SD2.fa'

afrh2sd1 = '/Users/hoeps/PhD/projects/huminvs/revision/recurrent_invs_bps/chr8_fine/fasta_afr/SDs/HG03125h2_SD1.fa'
afrh2sd2 = '/Users/hoeps/PhD/projects/huminvs/revision/recurrent_invs_bps/chr8_fine/fasta_afr/SDs/HG00732h2_SD2.fa'

afrh1inv = '/Users/hoeps/PhD/projects/huminvs/revision/recurrent_invs_bps/chr8_fine/fasta_afr/SDs/HG03125h1_INV.fa'
afrh2inv = '/Users/hoeps/PhD/projects/huminvs/revision/recurrent_invs_bps/chr8_fine/fasta_afr/SDs/HG03125h2_INV.fa'

  
outpaf_link = './blub5'


make_chunked_minimap_alnment(targetfasta = HG00732sd1,
                             queryfasta = afrh2sd2,
                             outpaf_link,
                             chunklen = 1000, minsdlen = 2000, saveplot=F,
                             hllink = F, hltype = F)#, wholegenome = F)

make_chunked_minimap_alnment(targetfasta = HG00731sd1,
                             queryfasta = HG00731sd2,
                             outpaf_link,
                             chunklen = 500, minsdlen = 2000, saveplot=F,
                             hllink = F, hltype = F)#, wholegenome = F)


make_dotplot(HG00732sd2, HG00731sd2, 100, save=F,
                          targetrange=c(0, 20000),
                          queryrange=c(0 , 20000))


  # inside = cpaf_i[(cpaf_i$qstart <= start) & (cpaf_i$qend >= end),]
  # if (dim(inside)[1] == 1){
  #   print('Found intact contig.')
  #   
  #   if (inside$strand == '+'){
  #     liftover_start = (start - inside$qstart) + inside$tstart
  #     liftover_end = inside$tend - (inside$qend - end)
  #   } else {
  #     liftover_start = (inside$qend - end) + inside$tstart
  #     liftover_end = inside$tend - (start - inside$qstart)
  #   }
  #   
  #   
  #   return(c(inside$tname, liftover_start, liftover_end))
  # }
#}

  
  
  
  
  
# # 
# # # Get hg38 fasta
# make_chunked_minimap_alnment(chm13fa, outfasta_aln, outpaf_link,
#                              chunklen = 1000, minsdlen = 2000, saveplot=F,
#                              hllink = F, hltype = F, wholegenome = T)
# 
# plot = wrapper_dotplot_with_alignment(seqname, start, end, hg38fa, chm13fa, 'hg38:chr1:12991825-13228346.fa', 'aln_seq.fa', 'deleteme.paf')
# 
# 
# 










# aaa = pafdotplot_make(
#   "/Users/hoeps/phd/projects/nahrcall/nahrchainer/seqbuilder/res/oa1.filter",
#   "/Users/hoeps/phd/projects/nahrcall/nahrchainer/seqbuilder/res/oa1.pdf",
#   keep_ref = 10000,
#   plot_size = 10,
#   hllink = F,
#   hltype = "paf",
#   minsdlen = 2000,
#   save = F
# )
