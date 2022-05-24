library(devtools)
devtools::load_all()

# Specify your inputs!
hg38_fa = "/Users/hoeps/PhD/projects/huminvs/genomes/hg38/hg38.fa"
aln_fa = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/alns/'
conversion_paf = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/liftover_custom/.paf'
interval = list(chr = 'chrY', start =  1 , end=26525500 , sample='T2T')


# Set parameters
params = list(
      # Interactive debug mode?
      debug = F,

      # If T, no dotplot condensation and SV calling is performed.
      plot_only = T,

      # Whole_chr just needed for chrY at the moment. 
      whole_chr = T,

      # Pad the x sequence? Typical: factor 1.5 - 3
      xpad = determine_xpad(interval$start, interval$end),

      # Alignment chunklen to use. Default: 1000 - 10000
      chunklen = determine_chunklen_compression(interval$start, interval$end),
      
      # Minlen of alingment to make it into plot (default: 50-2000)
      plot_minlen = determine_plot_minlen(interval$start, interval$end),
      
      # Options for compressed dotplots
      auto_find_compression = T,
      mode = 'precise',
      minlen = NULL, 
      compression = NULL,
      max_n_alns = 100,
      n_tests = 20,
      max_size_col_plus_rows = 200,
      baseline_log_minsize_min = log2(1),
      baseline_log_minsize_max = max(log2(20000), log2((interval$end-interval$start) / 10)),
      
      # Mutation search depth
      depth = 1,
      
      # Should alignments and fastas be deleted? (They can get quite large) 
      clean_after_yourself = F
    )
    

# Run! 
wrapper_aln_and_analyse(interval$chr,
                        interval$start,
                        interval$end,
                        hg38_fa,
                        aln_fa,
                        conversion_paf,
                        samplename = paste0(interval$sample, '-', interval$category, '-', nrun),
                        params = params,
                        logfile = 'res2/unittest_v4.tsv')
  }
}


