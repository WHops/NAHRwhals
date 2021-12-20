# blub 

samplefasta_link = '/Users/hoeps/phd/projects/nahrcall/nahrchainer/seqbuilder/data/realexample.fa'
#samplefasta_link = '/Users/hoeps/PhD/projects/huminvs/genomes/hg38/wbs.fa'
outpaf_link = '/Users/hoeps/phd/projects/nahrcall/nahrchainer/seqbuilder/res/outpaf6'
make_chunked_minimap_alnment(samplefasta_link, samplefasta_link, outpaf_link, 
                               outplot=NULL, chunklen = 100, minsdlen = 2000, saveplot=F, 
                               hllink = outpaf_link, hltype = 'paf')


# Problems/Errors spotted: 

# 1) If outpaf exists, problems appear. SOMETIMES. 

# 2) Merging algorithm dies if there are too many alignments. 
samplefasta_link = '/Users/hoeps/PhD/projects/huminvs/genomes/hg38/wbs.fa'
outpaf_link = '/Users/hoeps/phd/projects/nahrcall/nahrchainer/seqbuilder/res/outpaf9'
make_chunked_minimap_alnment(samplefasta_link, samplefasta_link, outpaf_link, 
                             outplot=NULL, chunklen = 50000, minsdlen = 500, saveplot=F, 
                             hllink = outpaf_link, hltype = 'paf')
                             

# 3) Alignments CAN break
samplefasta_link = '/Users/hoeps/PhD/projects/huminvs/genomes/hg38/s5.fa'
outpaf_link = '/Users/hoeps/phd/projects/nahrcall/nahrchainer/seqbuilder/res/outpaf9'
make_chunked_minimap_alnment(samplefasta_link, samplefasta_link, outpaf_link, 
                             outplot=NULL, chunklen = 5000, minsdlen = 5000, saveplot=F, 
                             hllink = outpaf_link, hltype = 'paf')

# 4) Plot groundtruth of SDs (should work, but we need to asjust coordinate systems so that the SD
# Annotations from groundtruth start at zero too...) 


# Future: 

# a) implement a tree search
# b) implement a quality metric simple alignment test
# c) implement a 'region search': region in, and need to find correct region on chr. 
