Rscript seqbuilder_wrapper.R -l 10000 -s ../data/sds.bed -o ../res/blub.fa -p none

Rscript seqbuilder_wrapper.R -l 10000 -s ../data/sds10y.tsv -o blub



../../../../../bbmap/shred.sh in=10k.fa out=10_chunk2k.fa length=2000

minimap2 -x asm20 -c -z400,50 -s 0 -M 0.2 -N 100 -P --hard-mask-level 10k.fa 10_chunk2k.fa > test

awk '{FS=OFS="\t"} match($1, /_([0-9]*)?-/,a){$1="blub"; $3 = $3+a[1]; $4 = $4+a[1]; print $0}' test > test2.paf

../../../R/pafCoordsDotPlotly.R -i test2.paf -o outplot -m 50 -p 10 -q 10
