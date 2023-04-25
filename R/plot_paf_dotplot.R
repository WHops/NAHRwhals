#' Add file-based highlights to a ggplot object
#'
#' @param gp A ggplot object to add file-based highlights to
#' @param opt A list containing necessary options for the highlighting
#' @return A ggplot object with file-based highlights added
#' @export
add_file_highlight <- function(gp, opt) {
  # Load in the annotation file based on the specified format
  if (opt$hltype == 'bed') {
    cat('Taking a bedfile as highlighting input\n')
    highlight_data = read.table(opt$hllink, sep = '\t')
    colnames(highlight_data) =  c(
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
    cat('Taking an sdfile as highlighting input\n')
    highlight_data = sd_to_bed(opt$hllink, outbedfile = NULL)
  } else if (opt$hltype == 'paf') {
    cat('Taking a paffile as highlighting input\n')
    highlight_data = paf_write_bed(opt$hllink, outsdbed_link = NULL)
    
    # Modify column names to plot target on X-axis and query on Y-axis
    colnames(highlight_data) = c(
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
  
  # Filter out rows with length less than the specified minimum length
  highlight_data = highlight_data[((highlight_data$chromEnd - highlight_data$chromStart) > opt$minsdlen),]
  
  cat(paste0('Printing a final thing with ', nrow(highlight_data), ' entries.\n'))
  
  # Add geom_rect layers to the ggplot object for the file-based highlights
  gp = gp + ggplot2::geom_rect(
    data = highlight_data,
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
      data = highlight_data,
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
      data = highlight_data,
      ggplot2::aes(
        xmin = chromStart,
        xmax = otherEnd,
        ymin = otherStart,
        ymax = otherEnd
      ),
      alpha = 0.25,
      fill = 'grey'
    )
  
  return(gp)
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
#' @author Tom Poorten (https://github.com/tpoorten/dotPlotly), edited by Wolfram Höps.
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



#' Plot alignments
#'
#' @param alignments A data frame containing alignment information
#' @param opt A list of options for the plot
#' @export
plot_alignments <- function(alignments, opt) {
  # Print x_end and x_start
  print(opt$x_end)
  print(opt$x_start)

  # Calculate yTickMarks
  yTickMarks = tapply(alignments$queryEnd, alignments$queryID, max)

  # Make y name
  stringlist = (lapply(stringr::str_split(as.character(alignments[, 'queryID']), '_'), head, -1))
  all_y_names = unlist(lapply(stringlist, paste, collapse = ''))
  common_y_name = names(which.max(table(all_y_names)))

  # Set samplename_x and samplename_y based on aln_type_xx_yy_xy
  if (opt$aln_type_xx_yy_xy == 'xx') {
    samplename_x = params$samplename_x
    samplename_y = params$samplename_x
  } else if (opt$aln_type_xx_yy_xy == 'yy') {
    samplename_x = params$samplename_y
    samplename_y = params$samplename_y
  } else if (opt$aln_type_xx_yy_xy == 'xy') {
    samplename_x = params$samplename_x
    samplename_y = params$samplename_y
  }

  # Create the ggplot object
  gp = ggplot2::ggplot(alignments) +
    ggplot2::geom_segment(
      ggplot2::aes(
        y = queryStart,
        yend = queryEnd,
        x = refStart,
        xend = refEnd,
        color = percentID,
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
    ggplot2::ylab(paste0('[', samplename_y, '] ', common_y_name)) +
    ggplot2::xlab(paste0('[', samplename_x, '] ',
                         paste0(translate_t2t_chr_to_readable(stringr::str_split(as.character(alignments$refID[1]), ':')[[1]][1]), ':', stringr::str_split(as.character(alignments$refID[1]), ':')[[1]][2]))) +
    ggplot2::coord_fixed() +
    ggplot2::scale_x_continuous(labels = scales::comma) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::labs(title = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())

  # Set x-axis scale and limits if x_start and x_end are not null
  if ((!is.null(opt$x_start)) & (!is.null(opt$x_end))) {
    gp = gp + ggplot2::scale_x_continuous(labels = scales::comma, limits = c(0, opt$x_end - opt$x_start)) +
      ggplot2::coord_fixed()
  }

  # If a link is given for highlighting
  if (opt$hllink != FALSE) {
    gp = add_file_highlight(gp, opt)
  }

  # If a (colored) highlight track is given
  if (opt$hltrack != FALSE) {
    gp = highlight_region(gp, opt)
  }

  # If a (colored) block track is given
  if (opt$anntrack != FALSE) {
    gp = add_ann_blocks(gp, opt)
  }

  return(gp)
}


#' Add annotation blocks to a given ggplot object
#'
#' @param gp A ggplot object to add annotation blocks to
#' @param opt A list containing necessary options for the annotation
#' @return A ggplot object with annotation blocks added
#' @export
add_annotation_blocks <- function(gp, opt) {
  max_annotation_blocks = 1000
  
  # Read annotation data
  annotation = read.table(opt$anntrack, header = FALSE, comment.char = "")
  annotation_interval = annotation[annotation$V1 == opt$x_seqname,]
  
  # Filter annotation data based on given options
  filtered_annotation = annotation_interval[
    (annotation_interval$V2 > opt$x_start & annotation_interval$V2 < opt$x_end) | # Start is embedded
      (annotation_interval$V3 > opt$x_start & annotation_interval$V3 < opt$x_end) | # End is embedded
      (annotation_interval$V2 < opt$x_start & annotation_interval$V3 > opt$x_end), # Annotation spans the whole region
    ]
  
  no_annotations = nrow(filtered_annotation) == 0
  too_many_annotations = nrow(filtered_annotation) > max_annotation_blocks
  
  if (no_annotations | too_many_annotations) {
    return(gp)
  }
  
  # Calculate start and end positions for the plot
  filtered_annotation$plot_start = (filtered_annotation$V2 - opt$x_start)
  filtered_annotation$plot_end = (filtered_annotation$V3 - opt$x_start)
  filtered_annotation$plot_start[filtered_annotation$plot_start <= 0] = 1
  filtered_annotation$plot_end[filtered_annotation$plot_end > (opt$x_end - opt$x_start)] = (opt$x_end - opt$x_start) - 1
  
  if (is.null(filtered_annotation$V5)) {
    filtered_annotation$V5 = rainbow(nrow(filtered_annotation))
  }
  
  if (sum(!areColors(filtered_annotation$V5)) > 0) {
    filtered_annotation[!areColors(filtered_annotation$V5), ]$V5 = rainbow(length(filtered_annotation[!areColors(filtered_annotation$V5), ]$V5))
  }
  
  # Create the annotation plot
  annotation_plot = ggplot2::ggplot(filtered_annotation) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = plot_start,
                   xmax = plot_end,
                   ymin = 0,
                   ymax = 1,
                   fill = V5),
      alpha = 0.5
    ) +
    ggplot2::scale_fill_identity(guide = "legend", breaks = filtered_annotation$V5) +
    ggplot2::xlim(c(0, opt$x_end - opt$x_start)) +
    ggplot2::ylim(c(0, 5)) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = 'none',
                   axis.ticks = ggplot2::element_blank(),
                   panel.spacing.x = ggplot2::unit(1, "mm"),
                   axis.title.x = ggplot2::element_blank(),
                   strip.background.x = ggplot2::element_blank(),
                   strip.text.x = ggplot2::element_blank())
  
  # Add text labels to the annotation plot
  annotation_plot = annotation_plot + ggrepel::geom_text_repel(
    ggplot2::aes(x = ((V3 - opt$x_start) + (V2 - opt$x_start)) / 2, y = 2, label = V4),
    color = filtered_annotation$V5,
    max.overlaps = 10
  )
  
  labels = sub(".*/", "", c(opt$anntrack, 'plot'))
  
  # Combine the annotation plot with the input ggplot object (gp) using patchwork
  gp_out = (annotation_plot + ggplot2::coord_fixed(ratio = (opt$x_start - opt$x_end) / 50)) +
    (gp + ggplot2::coord_fixed()) +
    patchwork::plot_layout(ncol = 1)
  
  return(gp_out)
}


#' Highlight a specific region in a given ggplot object
#'
#' @param gp A ggplot object to highlight a specific region in
#' @param opt A list containing necessary options for the highlighting
#' @return A ggplot object with the specified region highlighted
#' @export
highlight_region <- function(gp, opt) {
  # Check if the hltrack option is of the required class
  stopifnot(class(opt$hltrack) == 'character')
  
  # Read annotation data
  annotation = read.table(opt$hltrack, header = FALSE, comment.char = "")
  annotation_interval = annotation[annotation$V1 == opt$x_seqname, ]
  
  # Filter annotation data based on given options
  filtered_annotation = annotation_interval[
    (annotation_interval$V2 > opt$x_start & annotation_interval$V2 < opt$x_end) | # Start is embedded
      (annotation_interval$V3 > opt$x_start & annotation_interval$V3 < opt$x_end) | # End is embedded
      (annotation_interval$V2 < opt$x_start & annotation_interval$V3 > opt$x_end), # Annotation spans the whole region
    ]
  
  # Check if there are no annotations
  if (nrow(filtered_annotation) == 0) {
    return(gp)
  }
  
  if (is.null(filtered_annotation$V5)) {
    filtered_annotation$V5 = rainbow(nrow(filtered_annotation))
  }
  
  if (sum(!areColors(filtered_annotation$V5)) > 0) {
    filtered_annotation[!areColors(filtered_annotation$V5), ]$V5 = rainbow(length(filtered_annotation[!areColors(filtered_annotation$V5), ]$V5))
  }
  
  # Calculate start and end positions for the plot
  filtered_annotation$plot_start = (filtered_annotation$V2 - opt$x_start)
  filtered_annotation$plot_end = (filtered_annotation$V3 - opt$x_start)
  filtered_annotation$plot_start[filtered_annotation$plot_start <= 0] = 1
  filtered_annotation$plot_end[filtered_annotation$plot_end > (opt$x_end - opt$x_start)] = (opt$x_end - opt$x_start) - 1
  
  # Add a geom_rect layer to highlight the specified region in the ggplot object
  gp = gp + ggplot2::geom_rect(
    data = filtered_annotation,
    ggplot2::aes(
      xmin = plot_start,
      xmax = plot_end,
      ymin = -0.5,
      ymax = Inf,
      fill = V5
    ),
    alpha = 0.5
  ) + ggplot2::scale_fill_identity(guide = "legend", breaks = filtered_annotation$V5) +
    ggplot2::theme(legend.position = 'none')
  
  return(gp)
}

#' Check if given character vector elements are valid color representations
#'
#' This function checks if the given input is a valid color name, hex code or RGB value.
#'
#' @param x A character vector of color names, hex codes or RGB values.
#' @return A logical vector of the same length as `x`, indicating whether each element is a valid color representation.
#'
#' @details The function replaces missing values with 'none' and tries to convert each element to an RGB matrix using `col2rgb`.
#' If successful, the element is considered a valid color representation.
#'
#' @examples
#' areColors(c("red", "#FF0000", "rgb(255, 0, 0)", "blue", "invalid"))
#' areColors(c("green", "#00FF00", "rgb(0, 255, 0)", NA))
#'
#' @export
areColors <- function(x) {
  x[is.na(x)] <- 'none'
  sapply(x, function(X) {
    tryCatch(is.matrix(col2rgb(X)), error = function(e) FALSE)
  })
}



#' @title Generate a Dotplot from a PAF File with Optional SD Highlight Annotation
#' @description This function generates a dotplot from a PAF file and saves it as a PDF.
#' It can also include an optional SD highlight annotation file for enhanced visualization.
#'
#' @param inpaf_link Character/link: Path to the PAF file to be plotted.
#' @param outplot_link Character/link: Path to the output dotplot PDF to be saved.
#' @param min_align Numeric: Minimum alignment length (default: 11).
#' @param min_query_aln Numeric: Minimum query alignment length (default: 11).
#' @param keep_ref Numeric: Number of sequences to keep for plotting.
#' @param similarity Logical: Indicate sequence similarity by color code (default: TRUE).
#' @param h_lines Logical: Plot horizontal lines (default: TRUE).
#' @param interactive Logical: Use interactive dotplot (Dotplotly)? Untested, keep at FALSE.
#' @param plot_size Numeric: Plot size in inches.
#' @param on_target Logical: Unknown parameter, default: TRUE.
#' @param v Logical: Unknown parameter, default: FALSE.
#' @param hllink Character/link: Path to an SD annotation file (BED, TSV, or PAF) to include as highlights in the plot.
#' @param hltype Character: Filetype of hllink. Can be 'NULL', 'bed', 'tsv', 'paf'.
#' @param hlstart Numeric: Optional start position for highlights.
#' @param hlend Numeric: Optional end position for highlights.
#' @param save Logical: Save the dotplot to the specified output path (default: TRUE).
#' @param minsdlen Numeric: Minimum length for SD annotations (default: 5000).
#' @param anntrack Logical: Include annotation tracks (default: FALSE).
#' @param x_seqname Character: Sequence name for X-axis.
#' @param x_start Numeric: Start position for X-axis.
#' @param x_end Numeric: End position for X-axis.
#' @param hltrack Logical: Include highlight track (default: NULL).
#' @param aln_type_xx_yy_xy Character: Alignment type (default: 'xy').
#' @return The generated dotplot. If 'save' is TRUE, the dotplot is saved as a PDF.
#' @author Wolfram Höps
#' @export
pafdotplot_make <- function(inpaf_link,
                            outplot_link,
                            min_align = 11,
                            min_query_aln = 11,
                            keep_ref = 10000,
                            similarity = TRUE,
                            h_lines = TRUE,
                            interactive = FALSE,
                            plot_size = 10,
                            on_target = TRUE,
                            v = FALSE,
                            hllink = FALSE,
                            hltype = FALSE,
                            hlstart = NULL,
                            hlend = NULL,
                            save = TRUE,
                            minsdlen = 5000,
                            anntrack = FALSE,
                            x_seqname = NULL,
                            x_start = NULL,
                            x_end = NULL,
                            hltrack = NULL,
                            aln_type_xx_yy_xy = 'xy') {
  
  options <- list(
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
    hltrack = hltrack,
    aln_type_xx_yy_xy = aln_type_xx_yy_xy
  )

  # Obtain alignments using the provided options
  alignments = dotplotly_dotplot_return_aln(opt)
  # Generate the dotplot from the alignments
  plot = plot_alignments(alignments, opt)

  return(plot)
                            }

#' Plot a pairwise alignment format (PAF) data frame.
#'
#' This function creates a plot of a pairwise alignment format (PAF) data frame, showing the alignment between two sequences.
#'
#' @param paf_df a data frame containing PAF data. Required columns are: "tname", "tstart", "tend", "qname", "qstart", "qend", and "strand".
#' @param gridlines.x numeric vector of positions along the x-axis where vertical gridlines should be drawn. Default is NULL.
#' @param gridlines.y numeric vector of positions along the y-axis where horizontal gridlines should be drawn. Default is NULL.
#' @param wiggle numeric value indicating the amount of random jitter to add to the gridlines. Default is 0.
#' @param plot_gridlines logical value indicating whether or not to plot the gridlines. Default is TRUE.
#' @param linewidth numeric value indicating the width of the plotted gridlines. Default is 1.
#'
#' @return a ggplot object containing the alignment plot.
#'
#' @author Wolfram Hoeps
#' @export
plot_paf <- function(paf_df,
                     gridlines.x = NULL,
                     gridlines.y = NULL,
                     wiggle = 0,
                     plot_gridlines = T, 
                     linewidth = 1) {
  plot = ggplot2::ggplot() +
    ggplot2::geom_segment(data = paf_df,
                          ggplot2::aes(
                            x = tstart,
                            xend = tend,
                            y = qstart,
                            yend = qend
                          )) +
    ggplot2::coord_fixed(
      ratio = 1,
      xlim = NULL,
      ylim = NULL,
      expand = TRUE,
      clip = "on"
    ) +
    ggplot2::theme_bw() + 
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())
  
  if (plot_gridlines) {
    plot = plot +
      ggplot2::geom_hline(ggplot2::aes(yintercept = gridlines.y + runif(
        length(gridlines.y),-wiggle, wiggle
      ))
      , color =
        'grey', size = linewidth) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = gridlines.x + runif(
        length(gridlines.x), wiggle, wiggle
      ))
      , color =
        'grey', size = linewidth) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())
  }
  return(plot)
}
