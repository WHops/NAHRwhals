

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

fa_link = '/Users/hoeps/chr10:45.fa'
chm13_link = '/Users/hoeps/PhD/projects/huminvs/genomes/CHM13_T2T/fasta/chm13.draft_v1.1.fasta'


make_chunked_minimap_alnment(targetfasta=chm13_link, queryfasta=fa_link, './chr10:45-2.fa.paf', 
                             outplot=NULL, chunklen = 10000, minsdlen = 2000, saveplot=F, 
                             hllink = samplepaf_link, hltype = 'paf', quadrantsize = 100000, wholegenome = T)

grid = wrapper_paf_to_bitlocus('./chr1ff6:28471289ddfrrd2-228d2632765221.fa.paf', gp = 1)


samplefasta_link = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/vignettes/simulated_seq_10kb_4SDs.fa'
samplemutfasta_link = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/vignettes/simulated_seq_10kb_del.fa'
samplepaf_link = paste0('blub2.paf')
make_chunked_minimap_alnment(samplefasta_link, samplemutfasta_link, samplepaf_link, 
                             outplot=NULL, chunklen = 2000, minsdlen = 0, saveplot=F, 
                             hllink = samplepaf_link, hltype = 'paf', quadrantsize = 1000)
grid = wrapper_paf_to_bitlocus(samplepaf_link, gp = 1)



library(nahrtoolkit)
samplefasta_link = system.file('extdata', '10ktest.fa', package='nahrtoolkit')



# blub 
samplefasta_link = system.file('extdata', '10ktest.fa', package='nahrtoolkit')
samplepaf_link = paste0(samplefasta_link, 'sample.paf')
make_chunked_minimap_alnment(samplefasta_link, samplefasta_link, samplepaf_link, 
                             outplot=NULL, chunklen = 10000, minsdlen = 100, saveplot=F, 
                             hllink = samplepaf_link, hltype = 'paf', quadrantsize = 10000)


grid = wrapper_paf_to_bitlocus(samplepaf_link, gp = 1)


# samplefasta_link = '/Users/hoeps/phd/projects/nahrcall/nahrchainer/seqbuilder/data/realexample.fa'
# #samplefasta_link = '/Users/hoeps/PhD/projects/huminvs/genomes/hg38/wbs.fa'
# outpaf_link = '/Users/hoeps/phd/projects/nahrcall/nahrchainer/seqbuilder/res/x2'
# make_chunked_minimap_alnment('hg38', samplefasta_link, outpaf_link, 
#                                outplot=NULL, chunklen = 10000, minsdlen = 2000, saveplot=F, 
#                                hllink = outpaf_link, hltype = 'paf')
# vgrid = wrapper_paf_to_bitlocus(outpaf_link, minlen = 1000, gp = 1000)






samplefasta_link = '/Users/hoeps/phd/projects/nahrcall/nahrchainer/seqbuilder/data/realexample.fa'
#samplefasta_link = '/Users/hoeps/PhD/projects/huminvs/genomes/hg38/wbs.fa'
outpaf_link = '/Users/hoeps/phd/projects/nahrcall/nahrchainer/seqbuilder/res/x6'
make_chunked_minimap_alnment('/Users/hoeps/PhD/projects/huminvs/genomes/hifi-asm/HG00512_hgsvc_pbsq2-ccs_1000-hifiasm.h1-un.fasta', samplefasta_link, outpaf_link, 
                             outplot=NULL, chunklen = 1000, minsdlen = 2000, saveplot=F, 
                             hllink = F, hltype = F)
vgrid = wrapper_paf_to_bitlocus(outpaf_link, minlen = 1000, gp = 1000)
#'/Users/hoeps/PhD/projects/huminvs/genomes/hifi-asm/HG00512_hgsvc_pbsq2-ccs_1000-hifiasm.h1-un.fasta'
vv# Problems/Errors spotted: 

# 1) If outpaf exists, problems appear. SOMETIMES. 

# 2) Merging algorithm dies if there are too many alignments. 
samplefasta_link = '/Users/hoeps/PhD/projects/huminvs/genomes/hg38/wbs.fa'
outpaf_link = '/Users/hoeps/phd/projects/nahrcall/nahrchainer/seqbuilder/res/outpaf9'
make_chunked_minimap_alnment(samplefasta_link, samplefasta_link, outpaf_link, 
                             outplot=NULL, chunklen = 10000, minsdlen = 500, saveplot=F, 
                             hllink = outpaf_link, hltype = 'paf', quadrantsize = 200000)
grid = wrapper_paf_to_bitlocus(outpaf_link, minlen = 10000, gp = 10000)


# 3) Alignments CAN break
samplefasta_link = '/Users/hoeps/PhD/projects/huminvs/genomes/hg38/s5.fa'
outpaf_link = '/Users/hoeps/phd/projects/nahrcall/nahrchainer/seqbuilder/res/outpafsds132335'
make_chunked_minimap_alnment(samplefasta_link, samplefasta_link, outpaf_link, 
                             outplot=NULL, chunklen = 10000, minsdlen = 10000, saveplot=F, 
                             hllink = outpaf_link, hltype = 'paf', quadrantsize = 200000)
grid = wrapper_paf_to_bitlocus(outpaf_link, minlen = 10000, gp = 10000)


### Stuff from Introduction
outfasta = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/vignettes/simulated_seq_twonest.fa'
outmutfasta = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/vignettes/simulated_seq_twonest_inv.fa'
outmutfasta2 = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/vignettes/simulated_seq_twonest_inv_dup.fa'
samplepaf_link = 'blub'
samplepaf_link2 = 'blub2'
samplepaf_link3 = 'blub3'

make_chunked_minimap_alnment(outfasta, outfasta, samplepaf_link, 
                             outplot=NULL, chunklen = 1000, minsdlen = 0, saveplot=F, 
                             hllink = samplepaf_link, hltype = 'paf', quadrantsize = 1000)
make_chunked_minimap_alnment(outfasta, outmutfasta, samplepaf_link2, 
                             outplot=NULL, chunklen = 1000, minsdlen = 0, saveplot=F, 
                             hllink = samplepaf_link, hltype = 'paf', quadrantsize = 1000)
make_chunked_minimap_alnment(outfasta, outmutfasta2, samplepaf_link3, 
                             outplot=NULL, chunklen = 1000, minsdlen = 0, saveplot=F, 
                             hllink = samplepaf_link, hltype = 'paf', quadrantsize = 1000)
grid = wrapper_paf_to_bitlocus(samplepaf_link, gp = 1)
grid = wrapper_paf_to_bitlocus(samplepaf_link2, gp = 1)
grid = wrapper_paf_to_bitlocus(samplepaf_link3, gp = 1)


# 4) Plot groundtruth of SDs (should work, but we need to asjust coordinate systems so that the SD
# Annotations from groundtruth start at zero too...) 


# Future: 

# a) implement a tree search
# b) implement a quality metric simple alignment test
# c) implement a 'region search': region in, and need to find correct region on chr. 

hg38fa = "/Users/hoeps/PhD/projects/huminvs/genomes/hg38/hg38.fa"
chm13fa = '/Users/hoeps/PhD/projects/huminvs/genomes/CHM13_T2T/fasta/chm13.draft_v1.1.fasta'
aln_fa = '/Users/hoeps/PhD/projects/huminvs/genomes/hifi-asm/HG00512_hgsvc_pbsq2-ccs_1000-hifiasm.h1-un.fasta'

outfasta_hg38 = "/Users/hoeps/phd/projects/nahrcall/nahrchainer/seqbuilder/res/hg38_sub.fa"
outfasta_aln = "/Users/hoeps/phd/projects/nahrcall/nahrchainer/seqbuilder/res/aln_sub.fa"

outpaf_link = '/Users/hoeps/phd/projects/nahrcall/nahrchainer/seqbuilder/res/1234.paf'

#hg38
seqname = 'chr1'
start = 12991825
end =   13528346
extract_subseq(hg38fa, seqname, start, end, outfasta_hg38)

# #HG00512
# seqname = 'cluster10_h1tg000002l'
# start = 3600000
# end = 4400000
# extract_subseq(aln_fa, seqname, start, end, outfasta_aln)
# 
# # Get hg38 fasta
make_chunked_minimap_alnment(outfasta_hg38, outfasta_aln, outpaf_link,
                             chunklen = 10000, minsdlen = 2000, saveplot=F,
                             hllink = F, hltype = F, wholegenome = F)

plot = wrapper_dotplot_with_alignment(seqname, start, end, hg38fa, chm13fa, 'hg38:chr1:12991825-13228346.fa', 'aln_seq.fa', 'deleteme.paf')













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
