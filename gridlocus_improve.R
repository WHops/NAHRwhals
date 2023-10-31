#Gridlocus improve
library(ggplot2)
library(microbenchmark)
library(profvis)
mb = microbenchmark
devtools::load_all()
paf = read.table('~/PhD/projects/nahrcall/revisions/WGS_mode/bigpaf.paf', header=T)
paf = read.table('res/chr1-144923346-150047635/diff/paf/T2T_NA12_xy.paf', header=T)
paf = read.table('~/Desktop/test.tsv', header=T)
# colnames_paf <- c(
#   "qname",
#   "qlen",
#   "qstart",
#   "qend",
#   "strand",
#   "tname",
#   "tlen",
#   "tstart",
#   "tend",
#   "nmatch",
#   "alen",
#   "mapq"
# )
# colnames(paf)[1:length(paf)] <- colnames_paf
# 
paf <- transform(
  paf,
  qend = ifelse(strand == "-", qstart, qend),
  qstart = ifelse(strand == "-", qend, qstart)
)

paf_df = load_and_prep_paf_for_gridplot_transform('res/chr1-144923346-150047635/diff/paf/T2T_NA12_xy.paf', minlen=10000, compression=10000)


paf_df = paf
gxy = make_xy_grid_fast(paf_df)

xdf = data.frame(x=gxy[[1]])
ydf = data.frame(y=gxy[[2]])

ggplot(paf_df) + geom_segment(aes(x=tstart, xend=tend, y=qstart, yend=qend)) +  geom_hline(data = ydf, aes(yintercept=y)) + geom_vline(data= xdf, aes(xintercept=x))#+  geom_hline(data = yspdf, aes(yintercept=y), color='red') + geom_vline(data= xspdf, aes(xintercept=x), color='red')


#microbenchmark(grid_list = fill_matrix_bressiexchange(paf_df, gxy),times = 2)
#grid_list = fill_matrix_bressiexchange(paf_df, gxy)

#grid_list = data.frame()

#microbenchmark(

#prof = profvis({bressiwrap(paf_df, gxy)}, interval=0.001)
microbenchmark(gl = bressiwrap(paf_df, gxy), times = 2)

#)
gl = bressiwrap(paf_df, gxy)
plot_matrix_ggplot(as.data.frame(gl))

write.table(mat, '~/Desktop/solvertest.tsv', col.names = F, row.names=F, quote=F, sep='\t')
write.table(t(diff(gxy[[2]])), '~/Desktop/solvertest_lens.tsv', col.names = F, row.names=F, quote=F, sep='\t')
write.table(t(diff(gxy[[1]])), '~/Desktop/solvertest_lens.tsv', append=T, col.names = F, row.names=F, quote=F, sep='\t')

julia_cmd = "julia solver.jl -d MAXDEPTH -u MAXDUP -w MAXWIDTH -r MINREPORT compressed_alignment compressed_lengths out_path"
system()













