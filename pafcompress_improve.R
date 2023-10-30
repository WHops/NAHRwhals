# Ok hang in there just a bit longer. 

devtools::load_all()

#inpaf_link = '/Users/hoeps/PhD/projects/nahrcall/revisions/WGS_mode/aln_lab/t4.paf'
inpaf_link = '/Users/hoeps/PhD/projects/nahrcall/revisions/WGS_mode/aln_lab/chr1/t5_hg38_chr1.paf'

correct_paf(inpaf_link, 't4_fix.paf')

paf = read.table(inpaf_link, fill=T)
paf <- read.table('t4_fix.paf')
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

compress_paf_fnct(inpaf_link = NULL, outpaf_link = 'out.paf', inpaf_df = paf, inparam_chunklen = 10000, n_quadrants_per_axis = 10)
microbenchmark(compress_paf_fnct(inpaf_link = NULL, outpaf_link = 'out.paf', inpaf_df = paf, inparam_chunklen = 10000, n_quadrants_per_axis = 1000), times=2)
