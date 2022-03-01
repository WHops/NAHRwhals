# Run that stuff.



# runs only when script is run by itself
if (sys.nframe() == 0) {
  
  library(optparse)
  #library(devtools)
  #devtools::load_all()
  library(nahrtoolkit)

  option_list = list(
    optparse::make_option(
      c("-f", "--fasta_x"),
      type = "character",
      default = NULL,
      help = "Fasta file for x axis",
      metavar = "numeric"
    ),
    optparse::make_option(
      c("-g", "--fasta_y"),
      type = "character",
      default = NULL,
      help = "Fasta file for x axis",
      metavar = "character"
    ),
    optparse::make_option(
      c("-p", "--paf"),
      type = "character",
      default = NULL,
      help = "Paffile connecting the two",
      metavar = "character"
    ),
    optparse::make_option(
      c("-o", "--outprefix"),
      type = "character",
      default = "./outputcorr/",
      help = "Run name",
      metavar = "character"
    ),
    optparse::make_option(
      c("-i", "--interval"),
      type = "character",
      default = NULL,
      help = "Interval on fasta_x. Tab-separated.",
      metavar = "character"
    )
  )
  
  opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
  
  debug=F
  if (debug){
    opt = list()
    opt$interval = "chr1-5000000-5100000"
    opt$fasta_x = " ~/PhD/projects/huminvs/genomes/hg38/hg38.fa"
    opt$fasta_y = " ~/Desktop/alns/HG00512_hgsvc_pbsq2-ccs_1000-hifiasm.h2-un.fasta"
    opt$paf = "data/liftover_custom/HG00512_hgsvc_pbsq2-ccs_1000-hifiasm.h2-un_hg38.paf"
    opt$outprefix = 'blub'
  }
  
  interval = strsplit(opt$interval, '-')[[1]]
  wrapper_aln_and_analyse(
    interval[1],
    as.numeric(interval[2]),
    as.numeric(interval[3]),
    opt$fasta_x,
    opt$fasta_y,
    opt$paf,
    runname = opt$outprefix,
    chunklen = 10000,
    compression = 1000
  )
  
  print('done!')
}
