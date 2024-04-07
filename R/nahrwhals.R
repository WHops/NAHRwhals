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





#' Nahrwhals: Comparative Genomics Analysis Tool
#'
#' This function facilitates comparative genomics analysis by comparing a reference genome against an assembly genome. It supports detailed analysis options, including region-specific analysis, and provides extensive visualization features.
#'
#' @param ref_fa Path to the reference genome FASTA file. [required]
#' @param asm_fa Path to the assembly genome FASTA file for comparison. [required]
#' @param outdir Basedirectory to store output files and plots. [required]
#' @param samplename_x Reference genome label for outputs [default: "Fasta_x"]
#' @param samplename_y Assembly genome label for outputs [default: "Fasta_y"]
#' @param region Region to genotype in the format "chr:start-end". Specify this or regionfile, not both. [optional]
#' @param regionfile Region to genotype in the format "chr:start-end". Specify this or region, not both. [optional]
#' @param anntrack Includes annotation tracks (e.g. genes) in dotplot. Specify as 4-column bedfile (col 4: displayed name) [default: FALSE]
#' @param depth depth of the BFS mutation search [default: 3]
#' @param eval_th Evaluation threshold as a percentage. [default: 98]
#' @param chunklen Sequence chunk length alignment, "default" means automatic determination based on sequence length. [default: auto-detection]
#' @param minlen Disregard alignments shorter than this from the segmentation. [default: 350]
#' @param compression Minimum segment length. [default: auto-detection]
#' @param max_size_col_plus_rows Maximum size for the alignment matrix [default: 250]
#' @param max_n_alns Maximum number of alignments for a single query sequence [default: 150]
#' @param self_plots Generates self-comparison plots for the assembly and ref [default: TRUE]
#' @param plot_only Skip BFS/genotyping, just create dotplots [default: FALSE]
#' @param use_paf_library Use an external PAF library for alignments [default: FALSE]
#' @param conversionpaf_link if use_paf_library: provide link to whole-genome paf here [default: FALSE]
#' @param maxdup BFS: max number of duplications per chain [default: 2]
#' @param init_width BFS: follow up only init_width best-scoring nodes per depth [default: 1000] 
#' @param region_maxlen Maximum length of a window that can be analysied. Larger windows will be split into overlapping fragments [default: 5000000]
#' @param testrun_std Test the NAHRwhals installation [default: FALSE]
#' @param threads Number of windows to parallel genotype [default: 1]
#' @param minimap_cores_per_thread Number of cores to give to minimap PER THREAD [default: 1]
#' @param genome_y_fa_mmi Path to pre-indexed assembly genome with minimap2 [default: "default"]
#' @return Performs the specified comparative genomics analysis, generating relevant output files and plots.
#' @export
#'
#' @examples
#' nahrwhals(ref_fa="path/to/ref_genome.fasta", asm_fa="path/to/asm_genome.fasta", outdir="path/to/output")
#'
nahrwhals <- function(ref_fa, asm_fa, outdir='res', region=NULL, regionfile=NULL, anntrack = FALSE, 
                      depth = 3, eval_th = 98, chunklen = "default", minlen = "default", compression = "default",
                      max_size_col_plus_rows = 250, max_n_alns = 150, 
                      self_plots = TRUE, plot_only = FALSE,
                      samplename_x = "Fasta_x", samplename_y = "Fasta_y",
                      use_paf_library = FALSE, conversionpaf_link = FALSE,
                      maxdup = 2, init_width = 1000, region_maxlen = 5000000, testrun_std = FALSE,
                      threads = 1, minimap_cores_per_thread = 1, genome_y_fa_mmi = "default", blacklist = NULL) {


  
                        
    # Avoid scientific notations
    options(scipen=999)
    library(doParallel)

    
    if (is.null(region) && is.null(regionfile) && testrun_std==FALSE){
        print('No region or regionfile provided. Running whole genome discovery mode.')
        NW_mode = 'whole-genome'
    } else if ((!is.null(region) && is.null(regionfile))|| testrun_std==TRUE){
        print('Coordinates provided. Genotyping that region.')
        NW_mode = 'genotype_1region'
    } else if (is.null(region) && !is.null(regionfile)){
        print('Regionfile provided. Genotyping all regions in that file.')
        NW_mode = 'genotype_regionfile'
    } else {
        stop('Please provide either a region or a regionfile, not both.')
    }

    
    # Take some of the more outlandish parameters out of the users hands.
    # TODO: a gardening trip to eradicate some of these parameters
    xpad = 1
    plot_minlen = 350
    maxlen_refine = 1e10
    n_tests = 10
    n_max_testchunks = 5
    baseline_log_minsize_min = 8
    baseline_log_minsize_max = FALSE
    discovery_exact = FALSE
    hltrack = FALSE
    hllink = FALSE
    aln_pad_factor = 1.0
    debug = FALSE
    clean_after_yourself = FALSE
    compare_full_fastas = FALSE
    testrun_fullfa = FALSE
    minreport = 0.98
    noclutterplots = T
    plot_xy_segmented = TRUE
    minimap2_bin = 'minimap2'
    bedtools_bin = 'bedtools'

    logfile = paste0(outdir, '/nahrwhals_res.tsv')
    res_bedfile = paste0(outdir, '/nahrwhals_res.bed')

    # Check and create fasta indices and minimap2 mmi.  
    if(!file.exists(paste0(ref_fa,'.fai'))){
        message("Fasta index ", paste0(ref_fa,'.fai'), " does not exist. Creating it now.")
        run_silent(paste0("samtools faidx ", ref_fa))
    }

    if(!file.exists(paste0(asm_fa,'.fai'))){
        message("Fasta index ", paste0(asm_fa,'.fai'), " does not exist. Creating it now.")
        run_silent(paste0("samtools faidx ", asm_fa))
    }

    if (genome_y_fa_mmi == 'default'){
        genome_y_fa_mmi = paste0(outdir,'/minimap_idxs/', basename(asm_fa), '.mmi')
    }

    create_mmi_if_doesnt_exists(asm_fa, genome_y_fa_mmi)

    # 


    if (NW_mode == 'whole-genome'){

        # If blacklist is not provided, warn the user that this is necessary for whole genomes with centromeres, such as humans.
        if (is.null(blacklist)){
            message('Warning: Whole Genome mode selected, but no blacklist was provided. If you are analysing large complex genomes (e.g. human), as mask for centromeres is required.')
        }

        # Step 0: we also need a ref mmi for whole genome. 
        genome_x_fa_mmi = paste0(outdir,'/minimap_idxs/', basename(ref_fa), '.mmi')
        create_mmi_if_doesnt_exists(ref_fa, genome_x_fa_mmi)

        # Step 1: Scan the whole genome for windows to genotype
        window_dir = paste0(outdir, '/whole-genome-windows/', samplename_x, '_', samplename_y)
        dir.create(window_dir, showWarnings = F, recursive = T)

        windows_file = paste0(window_dir,'/nahrwhals_windows.bed')
        regions_to_genotype_df = scan_for_windows(ref_fa, asm_fa, threads, windows_file, ref_masked_regions=blacklist) #Outfile is optional

        # For debug
        #regions_to_genotype_df = regions_to_genotype_df[1:10,]
        #regions_to_genotype_df$end = regions_to_genotype_df$start + 20000


        #print(paste0('Running nahrwhals on ', nrow(regions_to_genotype_df), ' regions.'))
        #regions_to_genotype_df = regions_to_genotype_df[regions_to_genotype_df$start > 15094630,]
        #regions_to_genotype_df = regions_to_genotype_df[1,]


        #regions_to_genotype_df$start = 18010158
        #regions_to_genotype_df$end = 18010158 + 500000
        print(regions_to_genotype_df)
        # Prepare parallel runs
        browser()
        # Explain to user that we are setting up nthreads parallel workers
        print(paste0('Initializing ', min(threads,nrow(regions_to_genotype_df)), ' parallel workers.'))

        suppressMessages(suppressWarnings({cl <- parallel::makeCluster(min(threads,nrow(regions_to_genotype_df)), outfile = "")}))
        suppressMessages(suppressWarnings({registerDoParallel(cl)}))

        # Initialize progress bar
        #pb <- txtProgressBar(0, nrow(regions_to_genotype_df), style = 3) 
        solver_path = system.file('extdata', 'scripts', 'scripts_nw_main', 'solver.jl', package='nahrwhals')

        # Run using parallel workers
        foreach::foreach(idx = 1:nrow(regions_to_genotype_df)) %dopar% {
            #setTxtProgressBar(pb, idx)
            message(paste0('Running nahrwhals on region ', idx, ' of ', nrow(regions_to_genotype_df)))
            suppressMessages(suppressWarnings(library(nahrwhals)))

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
                              minreport=minreport, init_width=init_width, minimap_cores=minimap_cores_per_thread, julia_solver_path=solver_path, 
                              silent=T)

            )

            print(paste0('Finished analysis of region ', idx, ' of ', nrow(regions_to_genotype_df)))


            }))


        }   

        stopCluster(cl)
        #close(pb)
        # Inform user
        message('Finished running nahrwhals on all regions.')

        tsv_to_bed_regional_dominance(logfile, res_bedfile)

    } else if (NW_mode == 'genotype_regionfile'){ 
        # Step 1: alternative: just load the genotype bedfile

        print(regionfile)
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
        solver_path = system.file('extdata', 'scripts', 'scripts_nw_main', 'solver.jl', package='nahrwhals')

        # Run using parallel workers
        foreach::foreach(idx = 1:nrow(regions_to_genotype_df)) %dopar% {
            #setTxtProgressBar(pb, idx)
            
            message(paste0('Running nahrwhals on region ', idx, ' of ', nrow(regions_to_genotype_df)))
            suppressMessages(suppressWarnings(library(nahrwhals))) #<--- wtf no #TODO

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
                              minreport=minreport, init_width=init_width, minimap_cores=minimap_cores_per_thread, julia_solver_path=solver_path, silent=T)

            
            )


            print(paste0('Finished analysis of region ', idx))

            }))


        }   

        stopCluster(cl)
        close(pb)


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
                    minreport=minreport, init_width=init_width,  minimap_cores=minimap_cores_per_thread))
    }
    

    # Step 4: Run 'regional dominance' analysis
    
        
        

}




