#' align_all_vs_all_using_minimap2
#' @param minimap2_bin Path to the minimap2 binary
#' @param reference_fasta Path to the reference fasta
#' @param asm_fasta Path to the assembly fasta
#' @param out_paf Path to the output paf file
#' @param threads Number of threads to use
#' @return Nothing
#' @export
align_all_vs_all_using_minimap2 <- function(minimap2_bin, reference_fasta, asm_fasta, out_paf, threads){

    # Search for a reference_mmi in the same folder as the reference_fasta. If it does not exist, create it.
    # for example for a file reference.fa, we want to create reference.mmi
    reference_mmi = paste0(dirname(reference_fasta), "/", sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(reference_fasta)), ".mmi")
    if(!file.exists(reference_mmi)){
        # Inform the user what happened and what we do 
        message("Reference index ", reference_mmi, " does not exist. Creating it now.")
        minimap2_indexing_command <- paste0(minimap2_bin, " -k 28 -w 255 -H -d ", reference_mmi, " ", reference_fasta)
        system(minimap2_indexing_command)
    } else {
        # Inform the user what happened and what we do
        message("Found existing reference index ", reference_mmi, ". Skipping re-calculation.")
    }

    # Check if the out_paf exists already. If it does, explain to the user and do not run the code. 
    if(file.exists(out_paf)){
        message("Output file ", out_paf, " already exists. Skipping minimap2.")
        return()
    }
    minimap2_cmd = paste0(minimap2_bin, " -x asm5 -t ",threads, " -c ", reference_mmi, " ", asm_fasta, " > ", out_paf)
    system(minimap2_cmd)

    # Inform the user
    message("Minimap2 finished successfully. The output is in ", out_paf)
}

#' extract_test_list_from_paf
#' @param all_vs_all_paf Path to the paf file containing all-vs-all alignments
#' @param out_dir Path to the output directory
#' @param genome_path Path to the genome fasta
#' @param bedtools_bin Path to the bedtools binary
#' @param merge_distance Distance to merge overlapping intervals
#' @return Nothing
#' @export
extract_test_list_from_paf <- function(all_vs_all_paf, out_dir, genome_path, bedtools_bin, merge_distance){

    command = paste0('bash scripts/wg_run_all.sh ', all_vs_all_paf, ' ', out_dir, ' ', genome_path, ' ', bedtools_bin, ' ', merge_distance)
    system(command)

    # Inform the user
    message("Done. The output is in ", out_dir)
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
        system(paste0("samtools faidx ", fasta))
    } else {
        message("Found existing fasta index ", paste0(fasta,'.fai'), ". Skipping re-calculation.")
    }
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
wga_write_interval_list <- function(ref_fa, asm_fa, outdir, merge_distance, threads,
                                    bedtools_bin = 'bedtools', 
                                    minimap2_bin = '/Users/hoeps/opt/anaconda3/envs/snakemake/envs/nahrwhalsAPR/bin/minimap2'
){
    system(paste0("mkdir -p ", outdir))
    allpaf = paste0(outdir, "/fullaln.paf")
    genome_file = paste0(outdir, "/ref.genome")

    # Make all-vs-all alignemnt
    align_all_vs_all_using_minimap2(minimap2_bin, ref_fa, asm_fa, allpaf, threads)
    # Write the genome file
    make_genome_file(ref_fa, genome_file)
    # Extract list of breakpoint clusters
    print(paste0('Using merge distance: ', merge_distance))
    extract_test_list_from_paf(allpaf, outdir, genome_file, bedtools_bin, merge_distance)

    # return the name of the final output list, which is in outdir/list_final.txt
    return(paste0(outdir, "/list_cut_final.bed"))
}
