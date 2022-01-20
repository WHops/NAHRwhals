library(nahrtoolkit)


seqname = 'chr7'
start = 74769950
end = 76058098
hg38fa = "/Users/hoeps/PhD/projects/huminvs/genomes/hg38/hg38.fa"
aln_fa = "/Users/hoeps/PhD/projects/huminvs/genomes/hifi-asm/HG00512_hgsvc_pbsq2-ccs_1000-hifiasm.h1-un.fasta"
subseqfasta_x =  'res/hg38.fa'
subseqfasta_y =  'res/aln.fa'
conversionpaf_link = "/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/liftover_custom/HG00512_hgsvc_pbsq2-ccs_1000-hifiasm.h1-un_hg38.paf"
outpaf_link = 'res/res.paf'

nahrtoolkit::wrapper_dotplot_with_alignment_fast(seqname, start, end, hg38fa, aln_fa, 
                               subseqfasta_y, subseqfasta_x, 
                               conversionpaf_link, outpaf_link)
