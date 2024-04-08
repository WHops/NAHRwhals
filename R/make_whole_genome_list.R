# #' align_all_vs_all_using_minimap2
# #' @param minimap2_bin Path to the minimap2 binary
# #' @param reference_fasta Path to the reference fasta
# #' @param asm_fasta Path to the assembly fasta
# #' @param out_paf Path to the output paf file
# #' @param threads Number of threads to use
# #' @return Nothing
# #' @export
# align_all_vs_all_using_minimap2 <- function(minimap2_bin, reference_fasta, asm_fasta, out_paf, threads){

#     # Search for a reference_mmi in the same folder as the reference_fasta. If it does not exist, create it.
#     # for example for a file reference.fa, we want to create reference.mmi
#     genome_x_fa_mmi = paste0(outdir,'/intermediate/mmis/', basename(asm_fa), '.mmi')
#     if(!file.exists(genome_x_fa_mmi)){
#         # Inform the user what happened and what we do 
#         #message("Reference index ", reference_mmi, " does not exist. Creating it now.")
#         minimap2_indexing_command <- paste0(minimap2_bin, " -k 28 -w 255 -H -d ", genome_x_fa_mmi, " ", reference_fasta)
#         run_silent(minimap2_indexing_command)
#     } #else {
#         # Inform the user what happened and what we do
#         #message("Found existing reference index ", reference_mmi, ". Skipping re-calculation.")
#     #}

#     # Check if the out_paf exists already. If it does, explain to the user and do not run the code. 
#     if(file.exists(out_paf)){
#         #message("Output file ", out_paf, " already exists. Skipping minimap2.")
#         return()
#     }
#     minimap2_cmd = paste0(minimap2_bin, " -x asm5 -t ",threads, " -c ", reference_mmi, " ", asm_fasta, " > ", out_paf)
#     run_silent(minimap2_cmd)

#     # Inform the user
#     #message("Minimap2 finished successfully. The output is in ", out_paf)
# }

#' align_all_vs_all_using_minimap2
#' @param minimap2_bin Path to the minimap2 binary
#' @param reference_fasta Path to the reference fasta
#' @param asm_fasta Path to the assembly fasta
#' @param out_paf Path to the output paf file
#' @param threads Number of threads to use
#' @return Nothing
#' @export
align_all_vs_all_using_minimap2_shred_merge <- function(awkscript, reference_fasta, asm_fasta, out_paf, threads, outdir){

    minimap2_bin = 'minimap2'
    bedtools_bin = 'bedtools'

    reference_mmi = paste0(outdir,'/../../minimap_idxs/', basename(reference_fasta), '.mmi')

    # Check if the out_paf exists already. If it does, explain to the user and do not run the code. 
    if(file.exists(out_paf)){
        #message("Output file ", out_paf, " already exists. Skipping minimap2.")
        return()
    }

    params=list(bedtools_bin=bedtools_bin, fasta_awk_script=awkscript)
    asm_fa_shredded = paste0(outdir,'/ref_fa_shredded/', basename(asm_fasta), '_10kbp_chunked.fa') #paste0(asm_fasta, "_10kbp_chunked.fa")
    dir.create(dirname(asm_fa_shredded), showWarnings = F, recursive = T)

    alignment_chunked = paste0(out_paf, "_10kbp_chunked.paf")
    alignment_chunked_corrected = paste0(out_paf, "_10kbp_chunked_corrected.paf")

    # if asm_fa_shredded does not exist:
    if (!file.exists(asm_fa_shredded)){
        shred_seq_bedtools_multifasta(asm_fasta, asm_fa_shredded, 10000, params)
    }

    if (!file.exists(alignment_chunked)){
        minimap2_cmd = paste0(minimap2_bin, " -t ",threads, " ", reference_mmi, " ", asm_fa_shredded, " > ", alignment_chunked)
        message("Running minimap2")
        message(minimap2_cmd)
        system(minimap2_cmd)
    }

    if (!file.exists(alignment_chunked_corrected)){
        correct_paf(alignment_chunked, alignment_chunked_corrected)
    }

    
    paf <- read.table(alignment_chunked_corrected)
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
        compress_paf_fnct(inpaf_link = NULL, outpaf_link = out_paf, inpaf_df = paf, inparam_chunklen = 10000, n_quadrants_per_axis = 10, qname=unique(paf$qname)[row])
    }
    # Inform the user
    #message("Minimap2 finished successfully. The output is in ", out_paf)
}


#' extract_test_list_from_paf
#' @param all_vs_all_paf Path to the paf file containing all-vs-all alignments
#' @param out_dir Path to the output directory
#' @param genome_path Path to the genome fasta
#' @param bedtools_bin Path to the bedtools binary
#' @param merge_distance Distance to merge overlapping intervals
#' @return Nothing
#' @export
extract_test_list_from_paf <- function(all_vs_all_paf, out_dir, out_file, genome_path, bedtools_bin, merge_distance, indel_ignore_distance, exclusion_mask='none'){

  
    bashscripts_base = system.file('extdata', 'scripts', 'scripts_whole_genome', package='nahrwhals')
    wg_run_all_script = system.file('extdata', 'scripts', 'scripts_whole_genome', 'wg_run_all.sh', package='nahrwhals')
    bedtools_bin = 'bedtools'

    command = paste0('bash ',wg_run_all_script, ' ', all_vs_all_paf, ' ', out_dir, ' ', out_file, ' ', genome_path, ' ', bedtools_bin, ' ', merge_distance, ' ', indel_ignore_distance, ' ', exclusion_mask, ' ', bashscripts_base)
    system(command)

    # Inform the user
}

  #' Turn fasta.fai into genome file
#' @param fasta_fai Path to the fasta.fai
#' @param out_path Path to the output file
#' @return Nothing
#' @export
make_genome_file <- function(fasta, out_path){

    # If fasta.fai file does not exist, create it
    if(!file.exists(paste0(fasta,'.fai'))){
        message("Fasta index ", paste0(fasta,'.fai'), " does not exist. Creating it now.")
        run_silent(paste0("samtools faidx ", fasta))
    } #else {
    #    message("Found existing fasta index ", paste0(fasta,'.fai'), ". Skipping re-calculation.")
    #}
    # Read the fai file
    fai = read.table(paste0(fasta,'.fai'), sep="\t", header=F, stringsAsFactors=F)
    # Extract the first column
    seqnames = fai[,1]
    # Extract the second column
    lengths = fai[,2]
    # Combine the two into a data.frame
    df = data.frame(seqnames, lengths, stringsAsFactors=F)
    # Write the data.frame to file
    write.table(df, out_path, sep="\t", quote=F, row.names=F, col.names=F)
}

#' @export
wga_write_interval_list <- function(ref_fa, asm_fa, outdir, outfile, merge_distance, indel_ignore_distance, exclusion_mask='none', threads = 1){
    run_silent(paste0("mkdir -p ", outdir , ' 2>&1'))

    allpaf = paste0(outdir, "/fullaln.paf")
    #allpaf = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/out.paf'
    genome_file = paste0(outdir, "/ref.genome")
    # Make all-vs-all alignemnt
    awkscript = system.file('extdata', 'scripts', 'scripts_nw_main', 'awk_on_fasta_gpt4.sh', package='nahrwhals')
    
    align_all_vs_all_using_minimap2_shred_merge(awkscript, ref_fa, asm_fa, allpaf, threads, outdir)
    
    # Write the genome file
    make_genome_file(ref_fa, genome_file)
    # Extract list of breakpoint clusters
    extract_test_list_from_paf(allpaf, outdir, outfile, genome_file, bedtools_bin, merge_distance, indel_ignore_distance, exclusion_mask)

    # return the name of the final output list, which is in outdir/list_final.txt
    # return(outfile)
}

#' @export
inform_about_tests <- function(test_list){
  tests = read.table(test_list, sep='\t')
  colnames(tests) = c('chr', 'start', 'end')

  # Write to the user the number of tests, and their min, max, median and mean length. 
  # Write the lenghts in kilobases, rounded.
  message("Number of tests: ", nrow(tests))
  message("Min length: ", round(min(tests$end - tests$start)/1000, 1), " kb")
  message("Max length: ", round(max(tests$end - tests$start)/1000, 1), " kb")
  message("Median length: ", round(median(tests$end - tests$start)/1000, 1), " kb")
  message("Mean length: ", round(mean(tests$end - tests$start)/1000, 1), " kb")

}

#' @export
make_karyogram <- function(test_list_file, genome_file, second_list = NULL, specified_text='', chr_min=10000000){

    # Quick and hacky
    tests = read.table(test_list_file, sep='\t')
    colnames(tests) = c('chr','start','end')
    bed_data = tests

    genome <- read.table(genome_file, header=F)
    genome = genome[genome$V2 > chr_min,]
    row.names(genome) = genome$V1
    genome <- na.omit(genome[customChromosomeOrder(genome$V1), ])

    # 1. Compute statistics from your BED data
    num_regions <- nrow(tests)
    min_size <- min(tests$end - tests$start + 1)
    max_size <- max(tests$end - tests$start + 1)
    median_size <- median(tests$end - tests$start + 1)

    


    if (!is.null(second_list)){
        tests2 = read.table(second_list, sep='\t')
        colnames(tests2) = c('chr','start','end')
    }
    # 2. Add those statistics and a user-specified text to the karyoplot
    label_text <- paste0(specified_text, ' | ',
                        "Regions:", num_regions, ' | ',
                        " Size:", round(min_size/1000, 1), ' - ', round(max_size/1000, 1), ' kbp ',
                        "(median: ", round(median_size/1000,1), ')')


    endsize = 500000
    # tests, any region smaller than 100 kbp is elongated on both sides so it has a total length of 100000
    tests$len = tests$end - tests$start + 1
    tests$newstart = tests$start
    tests$newend = tests$end
    tests[tests$len < endsize, 'newstart'] = tests[tests$len < endsize, 'start'] - (endsize - tests[tests$len < endsize, 'len'])/2
    tests[tests$len < endsize, 'newend'] = tests[tests$len < endsize, 'end'] + (endsize - tests[tests$len < endsize, 'len'])/2
    tests$start = tests$newstart
    tests$end = tests$newend
    tests$newstart = NULL
    tests$newend = NULL

    tests2$len = tests2$end - tests2$start + 1
    tests2$newstart = tests2$start
    tests2$newend = tests2$end
    tests2[tests2$len < endsize, 'newstart'] = tests2[tests2$len < endsize, 'start'] - (endsize - tests2[tests2$len < endsize, 'len'])/2
    tests2[tests2$len < endsize, 'newend'] = tests2[tests2$len < endsize, 'end'] + (endsize - tests2[tests2$len < endsize, 'len'])/2

    tests2$start = tests2$newstart
    tests2$end = tests2$newend
    tests2$newstart = NULL
    tests2$newend = NULL

    pp <- karyoploteR::getDefaultPlotParams(plot.type=2)
    pp$data1outmargin <- 100
    pp$data2outmargin <- 100
    pp$topmargin <- 450
    kp <- karyoploteR::plotKaryotype(genome=genome, plot.type=2, plot.params = pp)
    karyoploteR::kpAddCytobandsAsLine(kp)

    karyoploteR::kpAddMainTitle(kp, label_text, cex=0.6)
    karyoploteR::kpPlotRegions(kp, data=tests, col="red", data.panel=2, r0=0, r1=1.5)
    # Color of this one should be a good contrast to the first one
    karyoploteR::kpPlotRegions(kp, data=tests2, col='black', data.panel=1, r0=0, r1=1.5)
    # Add a legend also
    # karyoploteR::kpAddLegend(kp, c("NAHRwhals whole-genome", "Yang et al. validated"), fill=c("#FF000088", "#0000FF88"), border=NA, lwd=0, cex=0.6, r0=0.5, r1=1, title="Legend")


}

#' @export
customChromosomeOrder <- function(x) {
  # Extract potential numeric parts of the chromosome names
  numeric_part <- as.numeric(gsub("^chr", "", x))
  
  # Check which are actual numeric chromosome names
  is_numeric <- !is.na(numeric_part)
  
  # Sort numeric chromosome names
  numeric_chromosomes <- x[is_numeric][order(numeric_part[is_numeric])]
  
  # Sort non-numeric chromosome names based on the desired order
  non_numeric_order <- c("X", "Y", "MT", "M")
  non_numeric_chromosomes <- setdiff(x, numeric_chromosomes)
  non_numeric_chromosomes <- non_numeric_chromosomes[order(match(non_numeric_chromosomes, c(paste0("chr", non_numeric_order), non_numeric_order)))]
  
  # Combine the sorted chromosome names in the desired order
  ordered_chromosomes <- c(numeric_chromosomes, non_numeric_chromosomes)
  
  # Return the ordered chromosome names
  return(ordered_chromosomes)
}

