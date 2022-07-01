#!/usr/bin/env Rscript


#' @export
add_file_highlight <- function(gp, opt){
  # if sdlink links to sd file (e.g. because it is a groundtruth):
  # Format to bed (... why exactly?)
  # sd_simple = sd_to_bed(opt$sdlink)
  
  # Load in the annotation file. Can be different formats.
  if (opt$hltype == 'bed') {
    print('Taking a bedfile as highlighting input')
    sd_simple = read.table(opt$hllink, sep = '\t')
    colnames(sd_simple) =  c(
      'chrom',
      'chromStart',
      'chromEnd',
      'uid',
      'otherChrom',
      'otherStart',
      'otherEnd',
      'strand',
      'fracMatch'
    )
  } else if (opt$hltype == 'sd') {
    print('Taking an sdfile as highlighting input')
    sd_simple = sd_to_bed(opt$hllink, outbedfile = NULL)
  } else if (opt$hltype == 'paf') {
    print('Taking a paffile as highlighting input')
    print(opt$hltype)
    sd_simple = paf_write_bed(opt$hllink, outsdbed_link = NULL)
    
    # [W] playing around with colnames a bit here. We want to plot
    # target on X, and query on Y. The paf format has
    # first query, and then target.
    colnames(sd_simple) = c(
      'otherName',
      'otherStart',
      'otherEnd',
      'uid',
      'chromName',
      'chromStart',
      'chromEnd',
      'strand',
      'id'
    )
  }
  sd_simple = sd_simple[((sd_simple$chromEnd - sd_simple$chromStart) > opt$minsdlen),]
  
  print(paste0('Printing a final thing with ', dim(sd_simple)[1], ' entries.'))
  # Assumes sd_simple is in BED format - each SD only once. (?)
  gp = gp + ggplot2::geom_rect(
    data = sd_simple,
    ggplot2::aes(
      xmin = chromStart,
      xmax = chromEnd,
      ymin = otherStart,
      ymax = otherEnd
    ),
    alpha = 0.25,
    fill = 'orange'
  ) +
    ggplot2::geom_rect(
      data = sd_simple,
      ggplot2::aes(
        xmin = chromStart,
        xmax = chromEnd,
        ymin = chromStart,
        ymax = otherEnd
      ),
      alpha = 0.25,
      fill = 'grey'
    ) +
    ggplot2::geom_rect(
      data = sd_simple,
      ggplot2::aes(
        xmin = chromStart,
        xmax = otherEnd,
        ymin = otherStart,
        ymax = otherEnd
      ),
      alpha = 0.25,
      fill = 'grey'
    )
  
}







## Make Dot Plot with Percent Divergence on color scale

#' A core function for plotting a paf.
#'
#' @description This is extracted and an edited version of
#' the script "pafCoordsDotPlotly.R" from the dotplotly package
#' (https://github.com/tpoorten/dotPlotly) by Tom Poorten.
#' Should be run through the 'pafdotplot_make' wrapper function
#' which helps in selecting input parameters.
#'
#' @param opt a list of parameters.
#' @return Nothing, but writes a pdf with a dotplot.
#'
#' @author Tom Poorten, edited by Wolfram Höps
#' @export
dotplotly_dotplot_return_aln <- function(opt) {
  if (opt$v) {
    cat(paste0("PARAMETERS:\ninput (-i): ", opt$input_filename, "\n"))
    cat(paste0("output (-o): ", opt$output_filename, "\n"))
    cat(paste0(
      "minimum query aggregate alignment length (-q): ",
      opt$min_query_aln,
      "\n"
    ))
    cat(paste0("minimum alignment length (-m): ", opt$min_align, "\n"))
    cat(paste0("plot size (-p): ", opt$plot_size, "\n"))
    cat(paste0("show horizontal lines (-l): ", opt$h_lines, "\n"))
    cat(paste0(
      "number of reference chromosomes to keep (-k): ",
      opt$keep_ref,
      "\n"
    ))
    cat(paste0("show % identity (-s): ", opt$similarity, "\n"))
    cat(paste0(
      "show % identity for on-target alignments only (-t): ",
      opt$similarity,
      "\n"
    ))
    cat(paste0("produce interactive plot (-x): ", opt$interactive, "\n"))
    cat(paste0("reference IDs to keep (-r): ", opt$refIDs, "\n"))
    
  }
  #opt$output_filename = unlist(strsplit(opt$output_filename, "/"))[length(unlist(strsplit(opt$output_filename, "/")))]
  
  debug = F
  if (debug) {
    opt = list()
    opt$min_query_aln = 0
    opt$min_align = 0
    opt$on_target = T
    opt$input_filename = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/res/paf/invs_mut_self.paf'
  }
  # read in alignments
  alignments = read.table(
    opt$input_filename,
    stringsAsFactors = F,
    fill = T,
    row.names = NULL,
    comment.char = ''
  )
  # read in originvbed
  #originv = read.table(opt$originvbed, stringsAsFactors = F)
  #colnames(originv) = c('chrom','start','end')
  # set column names
  # PAF IS ZERO-BASED - CHECK HOW CODE WORKS
  colnames(alignments)[1:12] = c(
    "queryID",
    "queryLen",
    "queryStart",
    "queryEnd",
    "strand",
    "refID",
    "refLen",
    "refStart",
    "refEnd",
    "numResidueMatches",
    "lenAln",
    "mapQ"
  )
  
  # Fixes for PAF
  # Some measure of similarity - need to check on this
  alignments$percentID = alignments$numResidueMatches / alignments$lenAln
  queryStartTemp = alignments$queryStart
  # Flip starts, ends for negative strand alignments
  alignments$queryStart[which(alignments$strand == "-")] = alignments$queryEnd[which(alignments$strand == "-")]
  alignments$queryEnd[which(alignments$strand == "-")] = queryStartTemp[which(alignments$strand == "-")]
  rm(queryStartTemp)
  cat(paste0("\nNumber of alignments: ", nrow(alignments), "\n"))
  cat(paste0("Number of query sequences: ", length(unique(
    alignments$queryID
  )), "\n"))
  
  # sort by ref chromosome sizes, keep top X chromosomes OR keep specified IDs
  if (is.null(opt$refIDs)) {
    chromMax = tapply(alignments$refEnd, alignments$refID, max)
    if (is.null(opt$keep_ref)) {
      opt$keep_ref = length(chromMax)
    }
    refIDsToKeepOrdered = names(sort(chromMax, decreasing = T)[1:opt$keep_ref])
    alignments = alignments[which(alignments$refID %in% refIDsToKeepOrdered), ]
    
  } else {
    refIDsToKeepOrdered = unlist(strsplit(opt$refIDs, ","))
    alignments = alignments[which(alignments$refID %in% refIDsToKeepOrdered), ]
  }
  
  # filter queries by alignment length, for now include overlapping intervals
  queryLenAgg = tapply(alignments$lenAln, alignments$queryID, sum)
  alignments = alignments[which(alignments$queryID %in% names(queryLenAgg)[which(queryLenAgg > opt$min_query_aln)]), ]
  
  # filter alignment by length
  alignments = alignments[which(alignments$lenAln > opt$min_align), ]
  
  # re-filter queries by alignment length, for now include overlapping intervals
  queryLenAgg = tapply(alignments$lenAln, alignments$queryID, sum)
  alignments = alignments[which(alignments$queryID %in% names(queryLenAgg)[which(queryLenAgg > opt$min_query_aln)]), ]
  
  cat(paste0(
    "\nAfter filtering... Number of alignments: ",
    nrow(alignments),
    "\n"
  ))
  cat(paste0(
    "After filtering... Number of query sequences: ",
    length(unique(alignments$queryID)),
    "\n\n"
  ))
  
  # sort df on ref
  alignments$refID = factor(alignments$refID, levels = refIDsToKeepOrdered) # set order of refID
  alignments = alignments[with(alignments, order(refID, refStart)), ]
  chromMax = tapply(alignments$refEnd, alignments$refID, max)
  
  # make new ref alignments for dot plot
  if (length(levels(alignments$refID)) > 1) {
    alignments$refStart2 = alignments$refStart + sapply(as.character(alignments$refID), function(x)
      ifelse(x == names((chromMax))[1], 0, cumsum(as.numeric(chromMax))[match(x, names(chromMax)) - 1]))
    alignments$refEnd2 = alignments$refEnd +     sapply(as.character(alignments$refID), function(x)
      ifelse(x == names((chromMax))[1], 0, cumsum(as.numeric(chromMax))[match(x, names(chromMax)) - 1]))
  } else {
    alignments$refStart2 = alignments$refStart
    alignments$refEnd2 = alignments$refEnd
    
  }
  
  ## queryID sorting step 1/2
  # sort levels of factor 'queryID' based on longest alignment
  alignments$queryID = factor(alignments$queryID, levels = unique(as.character(alignments$queryID)))
  queryMaxAlnIndex = tapply(alignments$lenAln,
                            alignments$queryID,
                            which.max,
                            simplify = F)
  alignments$queryID = factor(alignments$queryID, levels = unique(as.character(alignments$queryID))[order(mapply(
    function(x, i)
      alignments$refStart2[which(i == alignments$queryID)][x],
    queryMaxAlnIndex,
    names(queryMaxAlnIndex)
  ))])
  
  ## queryID sorting step 2/2
  ## sort levels of factor 'queryID' based on longest aggregrate alignmentst to refID's
  # per query ID, get aggregrate alignment length to each refID
  queryLenAggPerRef = sapply((levels(alignments$queryID)), function(x)
    tapply(alignments$lenAln[which(alignments$queryID == x)], alignments$refID[which(alignments$queryID == x)], sum))
  if (length(levels(alignments$refID)) > 1) {
    queryID_Ref = apply(queryLenAggPerRef, 2, function(x)
      rownames(queryLenAggPerRef)[which.max(x)])
  } else {
    queryID_Ref = sapply(queryLenAggPerRef, function(x)
      names(queryLenAggPerRef)[which.max(x)])
  }
  # set order for queryID
  alignments$queryID = factor(alignments$queryID, levels = (levels(alignments$queryID))[order(match(queryID_Ref, levels(alignments$refID)))])
  
  # get mean percent ID per contig
  #   calc percent ID based on on-target alignments only
  if (opt$on_target & length(levels(alignments$refID)) > 1) {
    alignments$queryTarget = queryID_Ref[match(as.character(alignments$queryID), names(queryID_Ref))]
    alignmentsOnTarget = alignments[which(as.character(alignments$refID) == alignments$queryTarget), ]
    scaffoldIDmean = tapply(alignmentsOnTarget$percentID,
                            alignmentsOnTarget$queryID,
                            mean)
    alignments$percentIDmean = as.numeric(scaffoldIDmean[match(as.character(alignments$queryID), names(scaffoldIDmean))])
    alignments$percentIDmean[which(as.character(alignments$refID) != alignments$queryTarget)] = NA
  } else{
    scaffoldIDmean = tapply(alignments$percentID, alignments$queryID, mean)
    alignments$percentIDmean = as.numeric(scaffoldIDmean[match(as.character(alignments$queryID), names(scaffoldIDmean))])
  }
  
  alignments = alignments[(abs(alignments$refEnd - alignments$refStart) > opt$minsdlen),]

  # and a filter for artifacts... :) 
  alignments = alignments[  between((abs(alignments$refEnd - alignments$refStart) / 
                            abs(alignments$queryEnd - alignments$queryStart)), 0.5, 2),  ]
  
  
  return(alignments)
}



#' @export
plot_alignments <-function(alignments, opt){
  ### Plot ###
  
  
  yTickMarks = tapply(alignments$queryEnd, alignments$queryID, max)

  gp = 
    ggplot2::ggplot(alignments) +
      ggplot2::geom_segment(
        ggplot2::aes(
          y = queryStart,
          yend = queryEnd,
          x = refStart,
          xend = refEnd,
          color = percentID,
          #sign(queryEnd-queryStart),
          text = sprintf(
            'Query ID: %s<br>Query Start Pos: %s<br>Query End Pos: %s<br>Target ID: %s<br>Target Start Pos: %s<br>Target End Pos: %s<br>Length: %s kb',
            queryID,
            queryStart,
            queryEnd,
            refID,
            refStart,
            refEnd,
            round(lenAln / 1000, 1)
          )
        )
      ) +
      ggplot2::labs(color = "Percent ID",
                    title = opt$input_filename) +
      ggplot2::xlab(as.character(alignments$refID[1])) +
      ggplot2::ylab(as.character(alignments$queryID[1])) +
      ggplot2::coord_fixed() + 
      ggplot2::scale_x_continuous(labels = scales::comma) +
      ggplot2::scale_y_continuous(labels = scales::comma) +
      # ggplot2::coord_fixed(
      #   ratio = 1,
      #   xlim = NULL,
      #   ylim = NULL,
      #   expand = TRUE,
      #   clip = "on"
      # ) +
      ggplot2::labs(title=NULL) + 
      ggplot2::theme_bw() + 
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())
    
  if ((!is.null(opt$x_start)) & (!is.null(opt$x_end))){
    gp = gp + 
      ggplot2::scale_x_continuous(labels = scales::comma, limits = c(0, opt$x_end-opt$x_start)) +
      ggplot2::coord_fixed()
    
  }
  
  # If a link is given for highlighting
  if (opt$hllink != F) {
    
    gp = add_file_highlight(gp, opt)
    
  }
  


  # If a (colored) highlight track is given...
  if (!is.null(opt$hltrack)){
    gp = highlight_region(gp, opt)

  }

  # If a (colored) block track is given...
  if (!is.null(opt$anntrack)){

    
    gp = add_ann_blocks(gp, opt)
  }
  


  return(gp)
  
}



#' @export
add_ann_blocks <- function(gp, opt){

    max_ann_blocks = 100
    library(patchwork)
    
    annotation = read.table(opt$anntrack, header=F, comment.char="")
    a_int = annotation[annotation$V1 == opt$x_seqname,]
    
    
    a_int2 = a_int[(a_int$V2 > opt$x_start &
                      a_int$V2 < opt$x_end) | # Start is Embedded
                     (a_int$V3 > opt$x_start &
                        a_int$V3 < opt$x_end) | # End is embedded
                     (a_int$V2 < opt$x_start &
                        a_int$V3 > opt$x_end)  ,] # Annotation spans whole region
    
    cond1 = dim(a_int2)[1] == 0
    cond2 = dim(a_int2)[1] > max_ann_blocks
    
    if (cond1 | cond2){
      return(gp)
    }
    
    a_int2$plotstart = (a_int2$V2 - opt$x_start)
    a_int2$plotend = (a_int2$V3 - opt$x_start)
    a_int2$plotstart[a_int2$plotstart <= 0] = 1
    a_int2$plotend[a_int2$plotend > (opt$x_end - opt$x_start)] = (opt$x_end - opt$x_start) - 1
    
    if (is.null(a_int2$V5)) {
      a_int2$V5 = rainbow(dim(a_int2)[1])
    }
    
    if (sum(!areColors(a_int2$V5)) > 0) {
      a_int2[!areColors(a_int2$V5), ]$V5 = rainbow(length(a_int2[!areColors(a_int2$V5), ]$V5))
    }
    
    h2 = ggplot2::ggplot(a_int2) + 
         ggplot2::geom_rect(
         ggplot2::aes(xmin= plotstart,
                     xmax= plotend,
                     ymin = 0,
                     ymax = 1,#0 + runif(dim(a_int2)[1], 0.75, 1.25),
                     fill = V5),
        alpha = 0.5
      ) + 
      ggplot2::scale_fill_identity(guide="legend",breaks=a_int2$V5) + 

      ggplot2::xlim(c(0,opt$x_end - opt$x_start)) + 
      ggplot2::ylim(c(0,2)) +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position='none',
                     axis.ticks = ggplot2::element_blank(),
                     panel.spacing.x = ggplot2::unit(1, "mm"),
                     axis.title.x = ggplot2::element_blank(),
                     strip.background.x = ggplot2::element_blank(),
                     strip.text.x = ggplot2::element_blank())
    
  if (dim(a_int2)[1] < 20){
    h2 = h2 + ggrepel::geom_text_repel(
      ggplot2::aes(x= ((V3 - opt$x_start) + (V2 - opt$x_start)) / 2,y = 1, label=V4),
      color=a_int2$V5,
      max.overlaps=100
    ) 
  }
  labels = sub(".*/", "", c(opt$anntrack, 'plot'))
  gp_out = 
    (h2 + ggplot2::coord_fixed(ratio = (opt$x_start - opt$x_end) / 50)) + (gp + ggplot2::coord_fixed()) + 
    patchwork::plot_layout(ncol=1)
  

  
  return(gp_out)
  
  # Filter annotation to interesting loci
  # annotation_int = 
  
}

#' @export
highlight_region <- function(gp, opt){
  
  #
  stopifnot(class(opt$hltrack) == 'character')
  annotation = read.table(opt$hltrack, header=F, comment.char="")
  a_int = annotation[annotation$V1 == opt$x_seqname, ]
  
  a_int2 = a_int[(a_int$V2 > opt$x_start &
                    a_int$V2 < opt$x_end) | # Start is Embedded
                 (a_int$V3 > opt$x_start &
                    a_int$V3 < opt$x_end) | # End is embedded
                 (a_int$V2 < opt$x_start &
                    a_int$V3 > opt$x_end)  ,] # Annotation spans whole region
  
  if (dim(a_int2)[1] == 0){
    return(gp)
  }
  if (is.null(a_int2$V5)) {
    a_int2$V5 = rainbow(dim(a_int2)[1])
  }
  
  if (sum(!areColors(a_int2$V5)) > 0) {
    a_int2[!areColors(a_int2$V5), ]$V5 = rainbow(length(a_int2[!areColors(a_int2$V5), ]$V5))
  }
  
  a_int2$plotstart = (a_int2$V2 - opt$x_start)
  a_int2$plotend = (a_int2$V3 - opt$x_start)
  a_int2$plotstart[a_int2$plotstart <= 0] = 1
  a_int2$plotend[a_int2$plotend > (opt$x_end - opt$x_start)] = (opt$x_end - opt$x_start) - 1
  
  
  gp = gp + ggplot2::geom_rect(
    data = a_int2,
    ggplot2::aes(
      xmin = plotstart,
      xmax = plotend,
      ymin = -0.5,
      ymax = Inf,
      fill = V5
    ),
    alpha = 0.5
  ) + ggplot2::scale_fill_identity(guide = "legend", breaks = a_int2$V5) +
      ggplot2::theme(legend.position = 'none')
  
  
  return(gp)
  
  
}

#' @export
areColors <- function(x) {
  x[is.na(x)] = 'none'
  sapply(x, function(X) {
    tryCatch(is.matrix(col2rgb(X)), 
             error = function(e) FALSE)
  })
}




#' Wrapperfunction for turning a paf (and an SD highlight annotation) into a dotplot pdf.
#'
#' @description T
#'
#' @param inpaf_link [character/link] a link to the paf to be plotted
#' @param outplot_link [character/link] a link to where the output pdf will be saved
#' @param min_align not sure. set to 11, keep it there.
#' @param min_query_aln not sure. set to 11, keep it there.
#' @param keep_ref number of sequences to keep or something. Dont touch.
#' @param similarity [T/F] indicate sequence similarity by color code (default: T)
#' @param h_lines [T/F] plot horizontal lines (default: T)
#' @param interactive [T/F] unsure. Dotplotly? Untested, keep at F.
#' @param plot_size [numeric] plot size in inch.
#' @param on_target [T/F] dont remember. default: T
#' @param v [T/F] unsure. default: F
#' @param hllink [character/link] link to an SD annotation file (bed, tsv or paf) to include as
#' highlights in the plot.
#' @param hltype [character] filetype of hllink. Can be 'NULL', 'bed', 'tsv', 'paf'.
#' @return Nothing, but writes a pdf with a dotplot.
#'
#' @author Wolfram Höps
#' @export
pafdotplot_make <-
  function(inpaf_link,
           outplot_link,
           min_align = 11,
           min_query_aln = 11,
           keep_ref = 10000,
           similarity = T,
           h_lines = T,
           interactive = F,
           plot_size = 10,
           on_target = T,
           v = F,
           hllink = F,
           hltype = F,
           hlstart = NULL,
           hlend = NULL,
           save = T,
           minsdlen = 5000,
           anntrack = NULL,
           x_seqname = NULL,
           x_start = NULL,
           x_end = NULL,
           hltrack = NULL) {
    
    # I deserve death penalty for this :)
    opt = list(
      input_filename = inpaf_link,
      output_filename = outplot_link,
      min_align = min_align,
      min_query_aln = min_query_aln,
      keep_ref = keep_ref,
      similarity = similarity,
      h_lines = h_lines,
      interactive = interactive,
      plot_size = plot_size,
      on_target = on_target,
      v = v,
      hllink = hllink,
      hltype = hltype,
      hlstart = hlstart,
      hlend = hlend,
      save = save,
      minsdlen = minsdlen,
      x_seqname = x_seqname,
      x_start = x_start,
      x_end = x_end,
      anntrack = anntrack,
      hltrack = hltrack
    )
    
    alignments = dotplotly_dotplot_return_aln(opt)
    plot = plot_alignments(alignments, opt)
    return(plot)
  }


# # runs only when script is run by itself
# if (sys.nframe() == 0) {
#   debugmode = F
#   if (!debugmode) {
#     option_list <- list(
#       optparse::make_option(
#         c("-i", "--input"),
#         type = "character",
#         default = NULL,
#         help = "coords file from mummer program 'show.coords' [default %default]",
#         dest = "input_filename"
#       ),
#       optparse::make_option(
#         c("-o", "--output"),
#         type = "character",
#         default = "out",
#         help = "output filename prefix [default %default]",
#         dest = "output_filename"
#       ),
#       optparse::make_option(
#         c("-v", "--verbose"),
#         action = "store_true",
#         default = TRUE,
#         help = "Print out all parameter settings [default]"
#       ),
#       optparse::make_option(
#         c("-q", "--min-query-length"),
#         type = "numeric",
#         default = 400000,
#         help = "filter queries with total alignments less than cutoff X bp [default %default]",
#         dest = "min_query_aln"
#       ),
#       optparse::make_option(
#         c("-m", "--min-alignment-length"),
#         type = "numeric",
#         default = 10000,
#         help = "filter alignments less than cutoff X bp [default %default]",
#         dest = "min_align"
#       ),
#       optparse::make_option(
#         c("-p", "--plot-size"),
#         type = "numeric",
#         default = 15,
#         help = "plot size X by X inches [default %default]",
#         dest = "plot_size"
#       ),
#       optparse::make_option(
#         c("-l", "--show-horizontal-lines"),
#         action = "store_true",
#         default = FALSE,
#         help = "turn on horizontal lines on plot for separating scaffolds  [default %default]",
#         dest = "h_lines"
#       ),
#       optparse::make_option(
#         c("-k", "--number-ref-chromosomes"),
#         type = "numeric",
#         default = NULL,
#         help = "number of sorted reference chromosomes to keep [default all chromosmes]",
#         dest = "keep_ref"
#       ),
#       optparse::make_option(
#         c("-s", "--identity"),
#         action = "store_true",
#         default = FALSE,
#         help = "turn on color alignments by % identity [default %default]",
#         dest = "similarity"
#       ),
#       optparse::make_option(
#         c("-t", "--identity-on-target"),
#         action = "store_true",
#         default = FALSE,
#         help = "turn on calculation of % identity for on-target alignments only [default %default]",
#         dest = "on_target"
#       ),
#       optparse::make_option(
#         c("-x", "--interactive-plot-off"),
#         action = "store_false",
#         default = TRUE,
#         help = "turn off production of interactive plotly [default %default]",
#         dest = "interactive"
#       ),
#       optparse::make_option(
#         c("-r", "--reference-ids"),
#         type = "character",
#         default = NULL,
#         help = "comma-separated list of reference IDs to keep [default %default]",
#         dest = "refIDs"
#       ),
#       optparse::make_option(
#         c("-j", "--originvbed"),
#         type = "character",
#         default = NULL,
#         help = "inversion length in bp",
#         dest = "originvbed"
#       )
#     )
#     
#     options(error = traceback)
#     
#     parser <-
#       optparse::OptionParser(usage = "%prog -i alignments.coords -o out [options]", option_list =
#                                option_list)
#     opt = optparse::parse_args(parser)
#     
#   } else if (debugmode == T) {
#     rm(list = ls())
#     opt = list(
#       input_filename = "~/s/g/korbel/hoeps/lab/assemblies/new_invextracts/chr1-26641304-INV-5445_backup/paf/chr1-26641304-INV-5445_HG00096_h1.paf",
#       output_filename = "~/s/g/korbel/hoeps/lab/assemblies/new_invextracts/chr10-55007429-55013435/dotplots/blub",
#       originvbed = "~/s/g/korbel/hoeps/lab/assemblies/new_invextracts/beds/chr1-26641304-INV-5445.bed",
#       min_align = 11,
#       min_query_aln = 11,
#       keep_ref = 23,
#       similarity = T,
#       h_lines = T,
#       interactive = F,
#       plot_size = 15,
#       on_target = T,
#       v = FALSE
#     )
#   }
#   
#   alignments = dotplotly_dotplot_return_aln(opt)
#   plot = plot_alignments(alignments)
#   
# }



#   gp = gp + ggplot2::geom_rect(
#     ggplot2::aes(
#       xmin = opt$hlstart,
#       xmax = opt$hlend,
#       ymin = -0.5,
#       ymax = Inf
#     ),
#     alpha = 0.005,
#     fill = 'orange'