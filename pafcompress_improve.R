# Ok hang in there just a bit longer. 

devtools::load_all()

#inpaf_link = '/Users/hoeps/PhD/projects/nahrcall/revisions/WGS_mode/aln_lab/t4.paf'
inpaf_link = '/Users/hoeps/PhD/projects/nahrcall/revisions/WGS_mode/aln_lab/chr1/t5_hg38_chr1.paf'
inpaf_link = '/Users/hoeps/PhD/projects/nahrcall/revisions/WGS_mode/aln_lab/NA12878_lab/NA12878_T2T.paf'
inpaf_link = '/Users/hoeps/PhD/projects/nahrcall/revisions/WGS_mode/aln_lab/HG00733_lab/HG00733_T2T_20.paf'
correct_paf(inpaf_link, 't4_fix.paf')

paf = read.table(inpaf_link, fill=T)
paf <- read.table('t4_fix.paf')
paf$V1 = sub("_[0-9]+-[0-9]+$", "", paf$V1)
colnames(paf) <- c(
  "qname",
  "qlen",
  "qstart",
  "qend",
  "strand",
  "tname",
  "tlen",
  "tstart",
  "tend",
  "nmatch",
  "alen",
  "mapq"
)

paf <- transform(
  paf,
  qend = ifelse(strand == "-", qstart, qend),
  qstart = ifelse(strand == "-", qend, qstart)
)


for (row in 1:length(unique(paf$qname))){
  compress_paf_fnct(inpaf_link = NULL, outpaf_link = 'out.paf', inpaf_df = paf, inparam_chunklen = 10000, n_quadrants_per_axis = 10, qname=unique(paf$qname)[row])
  print(row)
}
microbenchmark(compress_paf_fnct(inpaf_link = NULL, outpaf_link = 'out.paf', inpaf_df = paf, inparam_chunklen = 10000, n_quadrants_per_axis = 1000), times=2)


inpaf_link = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/out.paf'
paf = read.table(inpaf_link)
colnames(paf) <- c(
  "qname",
  "qlen",
  "qstart",
  "qend",
  "strand",
  "tname",
  "tlen",
  "tstart",
  "tend",
  "nmatch",
  "alen",
  "mapq"
)
library(ggplot2)






paf_cut = paf[(paf$tname=='chr7'),]
paf_cut <- transform(
  paf_cut,
  qend = ifelse(strand == "-", qstart, qend),
  qstart = ifelse(strand == "-", qend, qstart)
)

ggplot(paf_cut) + geom_segment(aes(x=tstart, xend=tend, y=qstart, yend=qend)) + ggplot2::coord_cartesian(xlim = c(54273650, 54591369), ylim=c(47250000, 47700000))




