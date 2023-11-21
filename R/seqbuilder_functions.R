#' A testfunction to see if nahrwhals is loaded
#' @author Wolfram Hoeps
#' @export
confirm_loaded_nw <- function() {
  print("NAHRwhals loaded and ready to go.")
}

#' Print version number.
#' @author Wolfram Hoeps
#' @export
version_nw <- function() {
  print("This is version 1.2X from Jul 31st, 2023.")
}

#' Perform Chunked Minimap2 Alignment and Generate Dotplot
#'
#' This function performs a chunked minimap2 alignment of two single-sequence FASTA files, generates a PAF file, and creates a dotplot PDF. It utilizes helper functions for sequence chunking, minimap2 alignment, PAF recombination, and dotplot generation.
#'
#' @param params a list of NAHRwhals paramters
#' @param targetfasta A character string containing the path to the 'target' single-sequence FASTA file (e.g., chm13).
#' @param queryfasta A character string containing the path to the 'query' single-sequence FASTA file.
#' @param chunklen A numeric value indicating the length of sequence chunks to split the query into (default = 1000 bp).
#' @param hllink A character string containing the path to an SD annotation file (BED, TSV, or PAF) to include as highlights in the plot (default = FALSE).
#' @param hltype A character string indicating the filetype of hllink. Can be 'NULL', 'bed', 'tsv', or 'paf' (default = FALSE).
#' @param outpaf A character string containing the path to the output PAF file to be written.
#' @param keep_ref A numeric value indicating the number of alignments to keep for plotting (default = 10000).
#' @param minsdlen A numeric value; purpose not entirely clear.
#' @param plot_size A numeric value indicating the plot size in inches (default = 10).
#' @param saveplot A logical value indicating whether or not to save the plot (default = TRUE).
#' @param savepaf A logical value indicating whether or not to save the PAF file (default = TRUE).
#' @param quadrantsize A numeric value indicating the quadrant size (default = 100000).
#' @param hlstart An integer value indicating the start position of the highlight region (default = NULL).
#' @param hlend An integer value indicating the end position of the highlight region (default = NULL).
#' @param targetrange A numeric vector containing the start and end positions of the target sequence to use (default = NULL).
#' @param queryrange A numeric vector containing the start and end positions of the query sequence to use (default = NULL).
#' @param anntrack A logical value indicating whether or not to include an annotation track (default = FALSE).
#' @param x_seqname A character string indicating the name of the reference sequence for the annotation track (default = NULL).
#' @param x_start A numeric value indicating the start position of the reference sequence for the annotation track (default = NULL).
#' @param x_end A numeric value indicating the end position of the reference sequence for the annotation track (default = NULL).
#' @param hltrack A numeric value indicating the track number to use for highlighting (default = NULL).
#' @param onlypafreturn A logical value indicating whether or not to only return the PAF file (default = FALSE).
#' @param aln_type_xx_yy_xy A character string indicating the orientation of the alignment (default = 'xy').
#'
#' @return This function has no return value. Output files are written to the specified paths.
#'
#' @author Wolfram Hoeps
#' @export
make_chunked_minimap_alnment <-
  function(params,
           targetfasta,
           queryfasta,
           outpaf,
           # outplot,
           chunklen = 1000,
           keep_ref = 10000,
           minsdlen = 2000,
           plot_size = 10,
           saveplot = T,
           savepaf = T,
           quadrantsize = 100000,
           hllink = F,
           hltype = F,
           hlstart = NULL,
           hlend = NULL,
           targetrange = NULL,
           queryrange = NULL,
           anntrack = F,
           x_seqname = NULL,
           x_start = NULL,
           x_end = NULL,
           hltrack = NULL,
           onlypafreturn = F,
           aln_type_xx_yy_xy = "xy") {

    # Define intermediate files
    queryfasta_chunk <- paste0(queryfasta, ".chunk.fa")
    outpaf_chunk <- paste0(outpaf, ".chunk")
    outpaf_awk <- paste0(outpaf, ".awked")
    outpaf_filter <- paste0(outpaf, ".filter")
    outpaf_plot <- paste0(outpaf, ".plot")
    outplot <- paste0(outpaf, ".pdf")

    # Run a series of chunking, aligning and merging functions/scripts
    # Single-sequence query fasta gets chopped into pieces.


    shred_seq_bedtools(queryfasta, queryfasta_chunk, chunklen, params)

    # Self explanatory
    run_minimap2(targetfasta, queryfasta_chunk, outpaf_chunk, params)

    # gawk is used to correct the sequence names. This is because I know only there
    # how to use regex...
    correct_paf(outpaf_chunk, outpaf_filter)



    # Check if abandon
    # pafin = read_and_prep_paf(outpaf_filter)
    # if (is_cluttered_paf(pafin) & (onlypafreturn == T | params$noclutterplots == T)){
    #   print('Ok we are in overly cluttery regions; abandon!')
    #   log_collection$cluttered_boundaries <<- T
    #   return(NULL)
    # }
    
    

    # paf of fragmented paf gets put back together.
    compress_paf_fnct(inpaf_link = outpaf_filter, outpaf_link = outpaf, inparam_chunklen = chunklen)

    system(paste0('rm ',outpaf_chunk))
    system(paste0('rm ',outpaf_filter))

    if (onlypafreturn) {
      return(outpaf)
    }

    # Make a dotplot of that final paf (and with sd highlighting).
    miniplot <- pafdotplot_make(
      outpaf,
      outplot,
      keep_ref = keep_ref,
      plot_size = plot_size,
      hllink = hllink,
      hltype = hltype,
      hlstart = hlstart,
      hlend = hlend,
      minsdlen = minsdlen,
      save = saveplot,
      anntrack = anntrack,
      x_start = x_start,
      x_end = x_end,
      x_seqname = x_seqname,
      hltrack = params$hltrack,
      aln_type_xx_yy_xy = aln_type_xx_yy_xy,
      samplename_x = params$samplename_x,
      samplename_y = params$samplename_y,
    )

    if (saveplot == F) {
      print("returning your plot")
      return(miniplot)
    } else {
      ggplot2::ggsave(
        filename = outplot,
        plot = miniplot,
        height = 10,
        width = 10,
        device = "pdf"
      )
      print("Saving")
    }
  }


#' Extracts the name and length of each sequence in a fasta file.
#'
#' This function takes an input fasta file and an output file as arguments,
#' and runs a 'gawk' command to extract the name and length of each sequence
#' in the fasta file. The extracted information is saved to the output file.
#'
#' @param inputfa A character string specifying the path to the input fasta file.
#' @param outputinfo A character string specifying the path to the output file.
#'
#' @return Nothing is returned by the function, but the sequence names and lengths are saved to the output file.
#' @author Wolfram Hoeps
#' @export
get_fasta_info <- function(inputfa, outputinfo) {
  # Run the gawk command to get the fasta file info
  cmd <- paste0(
    "gawk '/^>/{if (l!=\"\") print l; print; l=0; next}{l+=length($0)}END{print l}' ",
    inputfa, " > ", outputinfo
  )
  system(cmd)
}

#' Corrects a PAF file if it's made on scattered input.
#'
#' This function takes an input PAF file and an output PAF file as arguments,
#' and runs an 'awk' command to correct the input PAF file if it's made on
#' scattered input. The corrected PAF file is saved to the output file.
#'
#' @param inputpaf A character string specifying the path to the input PAF file.
#' @param outputpaf A character string specifying the path to the output PAF file.
#'
#' @return Nothing is returned by the function, but the corrected PAF file is saved to the output file.
#'
#' @author Wolfram Hoeps
#' @export
correct_paf <- function(inputpaf, outputpaf) {
  # Run the gawk command to correct the input paf file
  cmd <- paste0(
    "gawk '{FS=OFS=\"\\t\"} match($1, /_([0-9]*)?-/,a){$3 = $3+a[1]; $4 = $4+a[1]; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' ",
    inputpaf, " > ", outputpaf
  )
  system(cmd)
}

#' Chunkify query fasta
#'
#' @description This is a helperfunction calling bedtools to
#' chop a query sequence into chunks.
#'
#' @param infasta A link to a single-seq fasta to be chopped
#' @param outfasta_chunk A link to the output path for a chopped multi-seq fasta.
#' @param chunklen length of sequence chunks in bp
#' @param params a list with all NAHRwhals parameters
#' @return nothing. But output files written.
#'
#' @author Wolfram Hoeps
#' @export
shred_seq_bedtools <- function(infasta,
                               outfasta_chunk,
                               chunklen,
                               params) {
  bedtoolsloc <- params$bedtools_bin
  fasta_awk_script <- params$awkscript_fasta

  # Write a temporary bedfile that will be removed at the end of the function
  bed_tmp_file <- paste0("tmpbed_deleteme_", sprintf("%.0f", runif(1, 1e13, 1e14)), ".bed")
  fasta_awk_tmp_file <- paste0("tmpinfo_deleteme_", sprintf("%.0f", runif(1, 1e13, 1e14)), ".bed")

  # system(paste0(fasta_awk_script, ' ', infasta, ' ', fasta_awk_tmp_file))
  get_fasta_info(infasta, fasta_awk_tmp_file)
  # Collect information about the fasta we want to chop. This info is needed for bedtools getfasta
  # To work well.
  name_len_df <- read.table(fasta_awk_tmp_file)
  contigname <- sub(">", "", (name_len_df)[1, ])
  contiglen <- as.numeric((name_len_df)[2, ])

  bed_df <- data.frame(
    seqnames = contigname,
    start = sprintf("%d", seq(0, contiglen - (contiglen %% chunklen), by = chunklen)),
    end = sprintf("%d", pmin(seq(0, contiglen - (contiglen %% chunklen), by = chunklen) + (chunklen - 1), contiglen))
  )

  write.table(bed_df, file = bed_tmp_file, sep = "\t", quote = F, row.names = F, col.names = F)


  system(paste0("rm ", infasta, ".fai"))

  sedcmd <- "sed -r \'s/(.*):/\\1_/'"
  system(paste0(bedtoolsloc, " getfasta -fi ", infasta, " -bed ", bed_tmp_file, " | ", sedcmd, " > ", outfasta_chunk))

  system(paste0("rm ", bed_tmp_file))
  system(paste0("rm ", fasta_awk_tmp_file))
}


#' Submit a system command to run minimap2
#'
#' @description This is a helperfunction to run minimap2. Also check out:
#' https://github.com/PacificBiosciences/pbmm2/ . Minimap2 parameters:
#'  -k   k-mer size (no larger than 28). (-1)
#' -w   Minimizer window size. (1-)
#' -u   Disable homopolymer-compressed k-mer (compression is active for SUBREAD & UNROLLED presets).
#' -A   Matching score. (-1)
#' -B   Mismatch penalty. (-1)
#' -z   Z-drop score. (-1)
#' -Z   Z-drop inversion score. (-1)
#' -r   Bandwidth used in chaining and DP-based alignment. (-1)
#' -g   Stop chain enlongation if there are no minimizers in N bp. (-1)
#' #'
#' @param targetfasta A link to the 'target' single-sequence fasta (sometimes reference, e.g. chm13.)
#' @param queryfasta A link to the 'query' fasta. Can be single or multi-fasta
#' @param outpaf A Path to the output paffile to be written.
#' @param params a list with all NAHRwhals parameters
#' @param nthreads number of threads to use with Minimap2
#' @return nothing. Only output files written.
#'
#' @author Wolfram Hoeps
#' @export
run_minimap2 <-
  function(targetfasta,
           queryfasta,
           outpaf,
           params,
           nthreads = 4) {
    # system(paste0(minimap2loc," -x asm20 -c -z400,50 -s 0 -M 0.2 -N 100 -P --hard-mask-level ", fastatarget, " ", fastaquery, " > ", outpaf))

    minimap2loc <- params$minimap2_bin

    # Some self-defined parameters
    invisible(system(
      paste0(
        minimap2loc,
        " -x asm20 -P -c -s 0 -M 0.2 -t ",
        nthreads,
        " ",
        targetfasta,
        " ",
        queryfasta,
        " > ",
        outpaf
      )
    ))
    # Check if that was successful.
    stopifnot(
      "Alignment error: Minimap2 has not reported any significant alignment.
              Check if your input sequence is sufficiently long." =
        file.size(outpaf) != 0
    )
  }



#' Helperfunction to save a fasta file.
#'
#' @description a simple function that takes a data.frame that has a column name
#' and seq and writes a fasta file from it. Taken from
#' https://bootstrappers.umassmed.edu/guides/main/r_writeFasta.html.
#'
#' @param data A data frame with column 'name' and 'seq'
#' @param filename An output filename.
#'
#' @return nothing. Only output files written.
#'
#' @author Nicholas Hathaway
#' @export
writeFasta <- function(data, filename) {
  fastaLines <- c()
  for (rowNum in 1:nrow(data)) {
    fastaLines <- c(fastaLines, as.character(paste(">", data[rowNum, "name"], sep = "")))
    fastaLines <- c(fastaLines, as.character(data[rowNum, "seq"]))
  }
  fileConn <- file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}


#' https://rdrr.io/cran/insect/src/R/complement.R
#' Reverse complement DNA in character string format.
#'
#' This function reverse complements a DNA sequence or vector of DNA
#'   sequences that are stored as character strings.
#'
#' @param z a vector of DNA sequences in upper case character string format.
#' @return a vector of DNA sequences as upper case character strings.
#' @details This function accepts only DNA sequences in concatenated character
#'   string format.
#' @author Shaun Wilkinson
################################################################################
rc <- function(z) {
  rc1 <- function(zz) {
    s <- strsplit(zz, split = "")[[1]]
    s <- rev(s)
    dchars <- strsplit("ACGTMRWSYKVHDBNI", split = "")[[1]]
    comps <- strsplit("TGCAKYWSRMBDHVNI", split = "")[[1]]
    s <- s[s %in% dchars] # remove spaces etc
    s <- dchars[match(s, comps)]
    s <- paste0(s, collapse = "")
    return(s)
  }
  z <- toupper(z)
  tmpnames <- names(z)
  res <- unname(sapply(z, rc1))
  if (!is.null(attr(z, "quality"))) {
    strev <-
      function(x) {
        sapply(lapply(lapply(unname(x), charToRaw), rev), rawToChar)
      }
    attr(res, "quality") <-
      unname(sapply(attr(z, "quality"), strev))
  }
  names(res) <- tmpnames
  return(res)
}
################################################################################




#' colMax function
#'
#' Computes the maximum value for each column of a matrix or data frame.
#'
#' @param data a matrix or data frame.
#' @return a numeric vector containing the maximum value for each column of \code{data}.
#' @export
colMax <- function(data) {
  sapply(data, max, na.rm = TRUE)
}

#' shred_seq_bedtools_multifasta
#' @description This is a helperfunction calling bedtools to
#' chop a query sequence into chunks.
#'
#' @param infasta A link to a multi-seq fasta to be chopped
#' @param outfasta_chunk A link to the output path for a chopped multi-seq fasta.
#' @param chunklen length of sequence chunks in bp
#' @param params a list with all NAHRwhals parameters
#' @return nothing. But output files written.
#'
#' @author Wolfram Hoeps
shred_seq_bedtools_multifasta <- function(infasta,
                                outfasta_chunk,
                                chunklen,
                                params) {
  bedtoolsloc <- params$bedtools_bin
  fasta_awk_script <- params$awkscript_fasta
  # Write a temporary bedfile that will be removed at the end of the function
  bed_tmp_file <- paste0("tmpbed_deleteme_", sprintf("%.0f", runif(1, 1e13, 1e14)), ".bed")
  fasta_awk_tmp_file <- paste0("tmpinfo_deleteme_", sprintf("%.0f", runif(1, 1e13, 1e14)), ".bed")

  get_fasta_info_multi(infasta, fasta_awk_tmp_file)

  # Collect information about the fasta we want to chop. This info is needed for bedtools getfasta to work well.
  name_len_df <- read.table(fasta_awk_tmp_file, stringsAsFactors = FALSE)
  
  bed_df_list <- lapply(1:nrow(name_len_df), function(i) {
    contigname <- sub(">", "", name_len_df[i, 1])
    contiglen <- as.numeric(name_len_df[i, 2])
    data.frame(
      seqnames = contigname,
      start = sprintf("%d", seq(0, contiglen - (contiglen %% chunklen), by = chunklen)),
      end = sprintf("%d", pmin(seq(0, contiglen - (contiglen %% chunklen), by = chunklen) + (chunklen - 1), contiglen))
    )
  })

  bed_df <- do.call(rbind, bed_df_list)

  write.table(bed_df, file = bed_tmp_file, sep = "\t", quote = F, row.names = F, col.names = F)

  system(paste0("rm ", infasta, ".fai"))

  sedcmd <- "sed -r \'s/(.*):/\\1_/'"
  system(paste0(bedtoolsloc, " getfasta -fi ", infasta, " -bed ", bed_tmp_file, " | ", sedcmd, " > ", outfasta_chunk))

  system(paste0("rm ", bed_tmp_file))
  system(paste0("rm ", fasta_awk_tmp_file))
}

get_fasta_info_multi <- function(inputfa, outputinfo) {
  # Run the gawk command to get the fasta file info
  cmd <- paste0(
    "gawk '/^>/{if (l!=\"\") {print seqname, l} seqname=$0; l=0; next}{l+=length($0)}END{print seqname, l}' ",
    inputfa, " > ", outputinfo
  )
  system(cmd)
}

# in_fa = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/chunk_lab/NA12878_giab_pbsq2-ccs_1000-hifiasm.h1-un.fasta'
# out_fa = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/data/chunk_lab/NA12878_giab_pbsq2-ccs_1000-hifiasm.h1-un_chunked.fasta'
# params = list(bedtools_bin = 'bedtools', awkscript_fasta = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/scripts/awk_on_fasta_gpt4.sh')
# shred_seq_bedtools_multifasta(in_fa, out_fa, 10000, params)
