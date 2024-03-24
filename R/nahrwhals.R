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
nahrwhals <- function(ref_fa, asm_fa, outdir, region=NULL, regionfile=NULL, 
                      genome_y_fa_mmi = "default", anntrack = FALSE, samplename_x = "Fasta_x",
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
                      minimap_cores = 1, region_maxlen = 5000000 ){
  
                        
    # Avoid scientific notations
    options(scipen=999)
    library(doParallel)

    


  
    if (is.null(region) && is.null(regionfile)){
        NW_mode = 'whole-genome'
    } else if (!is.null(region) && is.null(regionfile)){
        NW_mode = 'genotype_1region'
    } else if (is.null(region) && !is.null(regionfile)){
        NW_mode = 'genotype_regionfile'
    } else {
        stop('Please provide either a region or a regionfile, not both.')
    }


    logfile = paste0(outdir, '/nahrwhals_res.tsv')
    res_bedfile = paste0(outdir, '/nahrwhals_res.bed')

    # Little intermezzo to check if the minimap2 index file exists.
    # If no, make it now, before the parallel workers start.
    if (params$genome_y_fa_mmi %in% default_param_values){
        genome_y_fa_mmi = paste0(params$outdir,'/intermediate/mmis/', basename(asm_fa), '.mmi')
    }
    
    mmi_params = c('minimap2_bin' = minimap2_bin, 
                    'genome_y_fa' = genome_y_fa,
                    'genome_y_fa_mmi' = genome_y_fa_mmi)

    create_mmi_if_doesnt_exists(mmi_params)


    #' Creates a minimap2 index (.mmi) file for a genome assembly if it doesn't exist.
#' @param params A list containing parameters for the function.
#' @export
create_mmi_if_doesnt_exists <- function(params) {
  if (file.exists(params$genome_y_fa_mmi)) {
    #print('Found existing minimap2 index ".mmi" file. Skipping re-calculation.')
    return()
  }


  print('No minimap2 index ".mmi" file of the assembly fasta found. Creating one now (takes around 1 minute for a whole genome assembly.)')

  minimap2_indexing_command <- paste0(params$minimap2_bin, " -k 28 -w 255 --idx-no-seq -H -d ", params$genome_y_fa_mmi, " ", params$genome_y_fa)
  run_silent(minimap2_indexing_command)
}



    if (NW_mode == 'whole-genome'){




        window_dir = paste0(outdir, '/whole-genome-windows/', samplename_x, '_', samplename_y)
        dir.create(window_dir, showWarnings = F, recursive = T)

        windows_file = paste0(window_dir,'/nahrwhals_windows.bed')
        regions_to_genotype_df = scan_for_windows(fasta_x, fasta_y, threads, windows_file) #Outfile is optional

        # For debug
        #regions_to_genotype_df = regions_to_genotype_df[1:10,]
        #regions_to_genotype_df$end = regions_to_genotype_df$start + 20000


        print(paste0('Running nahrwhals on ', nrow(regions_to_genotype_df), ' regions.'))
        print(regions_to_genotype_df)

        # Prepare parallel runs
        suppressMessages(suppressWarnings({cl <- parallel::makeCluster(min(threads,nrow(regions_to_genotype_df)), outfile = "")}))
        suppressMessages(suppressWarnings({registerDoParallel(cl)}))

        # Initialize progress bar
        pb <- txtProgressBar(0, nrow(regions_to_genotype_df), style = 3) 

        # Run using parallel workers
        foreach::foreach(idx = 1:nrow(regions_to_genotype_df)) %dopar% {

            setTxtProgressBar(pb, idx)

            suppressMessages(suppressWarnings(devtools::load_all()))

            row <- regions_to_genotype_df[idx, ]
            seqname_x <- row$chr
            start_x <- row$start
            end_x <- row$end
            suppressMessages(suppressWarnings({
            try(nahrwhals_singlerun(genome_x_fa=ref_fa, genome_y_fa=asm_fa, seqname_x=seqname_x, start_x=start_x, end_x=end_x, outdir=outdir,
                              genome_y_fa_mmi=genome_y_fa_mmi, anntrack=anntrack, logfile=logfile, samplename_x=samplename_x, 
                              samplename_y=samplename_y, compare_full_fastas=compare_full_fastas, plot_only=plot_only, 
                              self_plots=self_plots, plot_xy_segmented=plot_xy_segmented, eval_th=eval_th, depth=depth, 
                              chunklen=chunklen, minlen=minlen, compression=compression, max_size_col_plus_rows=max_size_col_plus_rows, 
                              max_n_alns=max_n_alns, use_paf_library=T, conversionpaf_link=paste0(window_dir,'/fullaln.paf'), 
                              xpad=xpad, plot_minlen=plot_minlen, maxlen_refine=maxlen_refine, n_tests=n_tests, 
                              n_max_testchunks=n_max_testchunks, baseline_log_minsize_min=baseline_log_minsize_min, 
                              baseline_log_minsize_max=baseline_log_minsize_max, discovery_exact=discovery_exact, hltrack=hltrack, 
                              hllink=hllink, aln_pad_factor=aln_pad_factor, debug=debug, clean_after_yourself=clean_after_yourself, 
                              testrun_std=testrun_std, testrun_fullfa=testrun_fullfa, noclutterplots=noclutterplots, maxdup=maxdup, 
                              minreport=minreport, init_width=init_width, minimap_cores=minimap_cores, silent=T)

            
            )



            }))


        }   

        stopCluster(cl)
        close(pb)
        # Inform user
        print('')
        print('Finished running nahrwhals on all regions.')

        #tsv_to_bed_regional_dominance(logfile, res_bedfile)

    } else if (NW_mode == 'genotype_regionfile'){ 
        # Step 1: alternative: just load the genotype bedfile
        regions_to_genotype_df = read.table(regionfile, sep='\t', header=F)
        colnames(regions_to_genotype_df) = c('chr','start', 'end')

        # Check if any region is longer than region_maxlen
        if (any(regions_to_genotype_df$end - regions_to_genotype_df$start > region_maxlen)){
            regions_to_genotype_df = split_regions(regions_to_genotype_df, region_maxlen, region_maxlen/2)
            print('Some regions were longer than region_maxlen. They were split into smaller regions.')
        }

        print(paste0('Running nahrwhals on ', nrow(regions_to_genotype_df), ' regions.'))
        print(regions_to_genotype_df)

        # Prepare parallel runs
        suppressMessages(suppressWarnings({cl <- parallel::makeCluster(min(threads,nrow(regions_to_genotype_df)), outfile = "")}))
        suppressMessages(suppressWarnings({registerDoParallel(cl)}))

        # Initialize progress bar
        pb <- txtProgressBar(0, nrow(regions_to_genotype_df), style = 3) 

        # Run using parallel workers
        foreach::foreach(idx = 1:nrow(regions_to_genotype_df)) %dopar% {

            setTxtProgressBar(pb, idx)

            suppressMessages(suppressWarnings(devtools::load_all()))

            row <- regions_to_genotype_df[idx, ]
            seqname_x <- row$chr
            start_x <- row$start
            end_x <- row$end
            suppressMessages(suppressWarnings({
            try(nahrwhals_singlerun(genome_x_fa=ref_fa, genome_y_fa=asm_fa, seqname_x=seqname_x, start_x=start_x, end_x=end_x, outdir=outdir,
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
                              minreport=minreport, init_width=init_width, minimap_cores=minimap_cores, silent=T)

            
            )



            }))


        }   

        stopCluster(cl)
        close(pb)
        # Inform user
        print('')
        print('Finished running nahrwhals on all regions.')

        tsv_to_bed_regional_dominance(logfile, res_bedfile)

    } else if (NW_mode == 'genotype_1region'){
        # Step 1: alternative: just load the genotype bedfile
        # Region is in format chr:start-end. Split it up.
        region = unlist(strsplit(region, ':|-'))
        seqname_x=region[1]
        start_x=as.numeric(region[2])
        end_x=as.numeric(region[3])

        try(nahrwhals_singlerun(genome_x_fa=ref_fa, genome_y_fa=asm_fa, seqname_x=seqname_x, start_x=start_x, end_x=end_x, outdir=outdir,
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
                    minreport=minreport, init_width=init_width, minimap_cores=minimap_cores))
    }
    

    # Step 4: Run 'regional dominance' analysis
    
        
        

}


#fasta_x = 'data/athaliana/An-1.chr.all.v2.0.fasta'
#fasta_y = 'data/athaliana/C24.chr.all.v2.0.fasta'
#threads <<-
#threads <<- 1 
#Sys.setenv(PATH = paste('/Users/hoeps/.juliaup/bin/', Sys.getenv("PATH"), sep=":"))
#Sys.setenv(PATH = paste("/Users/hoeps/opt/anaconda3/envs/snakemake/bin/", Sys.getenv("PATH"), sep=":"))
#Sys.setenv(PATH = paste("/Users/hoeps/opt/anaconda3/envs/nahrwhalsAPR/bin/", Sys.getenv("PATH"), sep=":"))
#devtools::load_all()
#nahrwhals(fasta_x, fasta_y)



