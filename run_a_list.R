# Run that stuff.



# runs only when script is run by itself
if (sys.nframe() == 0) {
  
  library(optparse)
  library(devtools)
  devtools::load_all()
  #library(nahrtoolkit)
  
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
      help = "Interval on fasta_x. Separated with minusees, like 'chr1-10000-20000'",
      metavar = "character"
    ),
    optparse::make_option(
      c("-e", "--example"),
      type = 'logical',
      default = F,
      help = "If T, run an example run to confirm NTK is working.",
      metavar = "character"
    )
  )
  
  opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

  if (opt$example){
    
    print('Running an example locus.')
    opt = list()
    opt$interval = "chr22_hg38_piece-400000-800000"
    opt$fasta_x = "/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/examples/chr22_example_hg38.fa"
    opt$fasta_y = "/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/examples/chr22_example_T2T.fa"
    opt$paf = "/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/examples/chr22_example_T2T_hg38.paf"
    opt$outprefix = 'example'
  }
  
  interval = strsplit(opt$interval, '-')[[1]]
  start = as.numeric(interval[2])
  end = as.numeric(interval[3])

  
  params = list(
    # Interactive debug mode?
    debug = F,
    
    # If T, no dotplot condensation and SV calling is performed.
    plot_only = F,
    
    # Whole_chr just needed for chrY at the moment.
    whole_chr = T,
    
    # Pad the x sequence? Typical: factor 1.5 - 3
    xpad = determine_xpad(start, end),
    
    # Alignment chunklen to use. Default: 1000 - 10000
    chunklen = determine_chunklen_compression(start, end),
    
    # Minlen of alingment to make it into plot (default: 50-2000)
    plot_minlen = determine_plot_minlen(start, end),
    
    # Options for compressed dotplots
    auto_find_compression = T,
    mode = 'precise',
    minlen = NULL,
    compression = NULL,
    max_n_alns = 100,
    n_tests = 20,
    max_size_col_plus_rows = 200,
    baseline_log_minsize_min = log2(1),
    baseline_log_minsize_max = max(log2(20000), log2((end-start) / 10)),
    
    # Mutation search depth
    depth = 1,
    
    # Should alignments and fastas be deleted? (They can get quite large)
    clean_after_yourself = F
    
  )
  
  wrapper_aln_and_analyse(interval[1],
                          start,
                          end,
                          opt$fasta_x,
                          opt$fasta_y,
                          opt$paf,
                          samplename = opt$outprefix,
                          params = params,
                          logfile = paste0('res/',interval[1],'-',interval[2],'-',interval[3],'/calls.tsv'))
  
  
  
  
  print('done!')
}