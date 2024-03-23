#!/usr/bin/env Rscript
#
# nahrwhals.R
#
# Created: 2020-07-15
# Author: Wolfram Hops
# Contact: wolfram.hoeps@gmail.com
#
# Description: This script is the main entry point for the nahrwhals package.

#devtools::load_all()


#' NW main
#' @export
#' 
nahrwhals <- function(ref_fa, asm_fa, regions=NULL,
                      genome_y_fa_mmi = "default", anntrack = FALSE,
                      logfile = "res/res.tsv", res_bedfile = 'res/res_all.bed', samplename_x = "Fasta_x",
                      samplename_y = "Fasta_y", compare_full_fastas = FALSE,
                      plot_only = FALSE, self_plots = TRUE,
                      plot_xy_segmented = TRUE, eval_th = 98, depth = 3,
                      chunklen = "default", minlen = "default", compression = "default",
                      max_size_col_plus_rows = 250, max_n_alns = 150,
                      use_paf_library = FALSE, conversionpaf_link = FALSE,
                      xpad = 1, plot_minlen = 350, maxlen_refine = 1e10,
                      n_tests = 10, n_max_testchunks = 5, baseline_log_minsize_min = 8,
                      baseline_log_minsize_max = FALSE, discovery_exact = FALSE,
                      hltrack = FALSE, hllink = FALSE, aln_pad_factor = 1.0,
                      debug = FALSE, clean_after_yourself = FALSE, testrun_std = FALSE,
                      testrun_fullfa = FALSE, threads = 1,
                      noclutterplots = T,   maxdup = 2, minreport = 0.98, init_width = 1000,
                      minimap_cores = 1 ){
  

    library(foreach)
  
    if (is.null(regions)){
        NW_mode = 'whole-genome'
    } else {
        NW_mode = 'genotype'
    }

    if (NW_mode == 'whole-genome'){
      # Step 1: define our windows
      regions_to_genotype_df = scan_for_windows(fasta_x, fasta_y, threads, outfile) #Outfile is optional
    } else if (NW_mode == 'genotype'){ 
      # Step 1: alternative: just load the genotype bedfile
      regions_to_genotype_df = read.table(regions, sep='\t', header=F)
      colnames(regions_to_genotype_df) = c('chr','start', 'end')
    }
    
    # Step 2: run each one
    cl <- parallel::makePSOCKcluster(threads)
    foreach::foreach(idx = 1:nrow(regions_to_genotype_df)) %dopar% {

        row <- regions_to_genotype_df[idx, ]
        seqname_x <- row$chr
        start_x <- row$start
        end_x <- row$end

        (nahrwhals_singlerun(genome_x_fa=ref_fa, genome_y_fa=asm_fa, seqname_x=seqname_x, start_x=start_x, end_x=end_x, 
                            genome_y_fa_mmi=genome_y_fa_mmi, anntrack=anntrack, logfile=logfile, samplename_x=samplename_x, 
                            samplename_y=samplename_y, compare_full_fastas=compare_full_fastas, plot_only=plot_only, 
                            self_plots=self_plots, plot_xy_segmented=plot_xy_segmented, eval_th=eval_th, depth=depth, 
                            chunklen=chunklen, minlen=minlen, compression=compression, max_size_col_plus_rows=max_size_col_plus_rows, 
                            max_n_alns=max_n_alns, use_paf_library=use_paf_library, conversionpaf_link=conversionpaf_link, 
                            xpad=xpad, plot_minlen=plot_minlen, maxlen_refine=maxlen_refine, n_tests=n_tests, 
                            n_max_testchunks=n_max_testchunks, baseline_log_minsize_min=baseline_log_minsize_min, 
                            baseline_log_minsize_max=baseline_log_minsize_max, discovery_exact=discovery_exact, hltrack=hltrack, 
                            hllink=hllink, aln_pad_factor=aln_pad_factor, debug=debug, clean_after_yourself=clean_after_yourself, 
                            testrun_std=testrun_std, testrun_fullfa=testrun_fullfa, noclutterplots=noclutterplots, maxdup=maxdup, 
                            minreport=minreport, init_width=init_width, minimap_cores=minimap_cores)
        )
    }

    # Stop the parallel backend after you're done to free up resources
    parallel::stopCluster(cl)
    
    # Step 4: Run 'regional dominance' analysis
    tsv_to_bed_regional_dominance(logfile, res_bedfile)
        
        

    print('Done!')


}


fasta_x = 'data/athaliana/An-1.chr.all.v2.0.fasta'
fasta_y = 'data/athaliana/C24.chr.all.v2.0.fasta'
#threads <<-
#threads <<- 1 
#Sys.setenv(PATH = paste('/Users/hoeps/.juliaup/bin/', Sys.getenv("PATH"), sep=":"))
#Sys.setenv(PATH = paste("/Users/hoeps/opt/anaconda3/envs/snakemake/bin/", Sys.getenv("PATH"), sep=":"))
#Sys.setenv(PATH = paste("/Users/hoeps/opt/anaconda3/envs/nahrwhalsAPR/bin/", Sys.getenv("PATH"), sep=":"))
#devtools::load_all()
#nahrwhals(fasta_x, fasta_y)



