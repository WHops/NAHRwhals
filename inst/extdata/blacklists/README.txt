# Whoeps, 7th April 2024

Creation of those blacklists: 


### t2t ###

wget https://hgdownload.soe.ucsc.edu/gbdb/hs1/censat/censat.bb
bigBedToBed censat.bb censat.bed
grep -v ct_ censat.bed | bedtools merge -i - -d 100000 | awk '($3-$2 > 1000000)' - > t2t_blacklist.bed

### hg38 ###

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz

gunzip gap.txt.gz centromeres.txt.gz

cut -f 2-5 gap.txt | awk '($3-$2 > 50000)' - | bedtools sort -i - > gaps_50kbp.bed
cut -f 2-5 centromeres.txt > hg38_centromere_cleaned.bed
cat gaps_50kbp.bed hg38_centromere_cleaned.bed | bedtools sort -i - | bedtools merge -d 1000 -i - > hg38_blacklist.bed

