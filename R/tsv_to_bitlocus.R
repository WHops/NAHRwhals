#' bounce point. Algorithm for making the grid. 
#' Description TBD. 
#' @author Wolfram Hoeps
#' @export
bounce_point <- function(vectors, point){
  

  newpoints = data.frame(matrix(ncol = 2, nrow = 0))
  newpoints[1,] = point
  colnames(newpoints) = c('x', 'y')
  
  # vertical overlaps
  y_overlap_vecs = vectors[((vectors$tstart < as.numeric(point[1])) &
                              (vectors$tend > as.numeric(point[1]))), ]
  
  # horizontal overlaps

  
  x_overlap_vecs = vectors[((vectors$qlow < as.numeric(point[2])) &
                              (vectors$qhigh > as.numeric(point[2]))), ]
  
  
  
  # If there is overlap, calculate overlap of point with every overlapper
  if (dim(y_overlap_vecs)[1] > 0) {
    for (i in 1:dim(y_overlap_vecs)[1]) {
      vec = y_overlap_vecs[i, ]
      newpoints = rbind(newpoints, c(point[1], ((vec$slope) * as.numeric(point[1])) + vec$y_intercept))
      #y_i = (point[1] + vec$x_intercept) / vec$slope_inv
    }
  }
  
  if (dim(x_overlap_vecs)[1] > 0) {
    for (i in 1:dim(x_overlap_vecs)[1]) {
      vec = x_overlap_vecs[i, ]
      newpoints = rbind(newpoints, 
                        c( ((point[2] - vec$y_intercept) / vec$slope), point[2])
      )
    }
  }
  return(newpoints)  
}

#' wrapper_paf_to_bitlocus
#' @author Wolfram Hoeps
#' @export
wrapper_paf_to_bitlocus <-
  function(inpaf,
           realplot = T,
           bitlocusplot = T,
           saveplot = F,
           minlen = 1000,
           compression = 1000) {
    # Read paf
    paf = read.table(inpaf)
    colnames(paf) = c(
      'qname',
      'qlen',
      'qstart',
      'qend',
      'strand',
      'tname',
      'tlen',
      'tstart',
      'tend',
      'nmatch',
      'alen',
      'mapq'
    )
    
    # A) PAF prep
    
    # Filter alignments by length
    paf = paf[paf$alen > minlen, ]
    
    # Round start/end by compression factor
    paf[,c('qstart','qend','tstart','tend')] = round(data.frame(paf[,c('qstart','qend','tstart','tend')] / compression),0) * compression
    
    # If any entry is not slope 1 yet (which is unlikely after compression), then make it so. 
    paf = enforce_slope_one(paf)
    
    # In paf, start is always smaller than end. For our slope etc calculation, it will be easier
    # to change the order of start and end, if the orientation of the alignment is negative. 
    paf = transform(
      paf,
      qend = ifelse(strand == '-', qstart, qend),
      qstart = ifelse(strand == '-', qend, qstart)
    )
    
    # Add slope information to the paf. 
    paf = cbind(paf, add_slope_intercept_info(paf))
    
    # The Transform from earlier shoots back at us here. We now have to 
    # make a new columns for qlow and qhigh - basically what was start and
    # end originally. 
    paf$qlow = -(apply(-paf[,c('qstart','qend')], 1, function(x) max(x)))
    paf$qhigh = apply(paf[,c('qstart','qend')], 1, function(x) max(x))
    
    # Make a grid, using the bounce algorithm. 
    print('making grid')
    gxy = make_xy_grid(paf, n_additional_bounces = 5)
    print('grid made')
    gridlines.x = gxy[[1]]
    gridlines.y = gxy[[2]]
    
    # Run bressi to fill the grid
    grid_list = data.frame()
    for (i in 1:dim(paf)[1]) {
      grid_list = rbind(grid_list, as.data.frame(
        bresenham(
          x = as.numeric(paf[i, c('tstart', 'tend')]) ,
          y = as.numeric(paf[i, c('qstart', 'qend')]) ,
          gridlines.x,
          gridlines.y,
          debug = F
        )
      ))
    }
    
    # Sort by x
    grid_list = grid_list[order(grid_list$x),]
    
    # EXPERIMENTAL PART #
    # Remove duplicates TRIPLE. First, remove obvious duplicates (all 3 same)
    # Then, identify those with same x/y but different z. The get 0.
    # Finally, the second step has produced new duplicates, so be mean again. 
    
    # Remove duplicates. This happens if SDs overlap (which can be a result
    # of the compression step)
    grid_list = grid_list[!duplicated(grid_list),]
    
    # If a value appears both positive and negative, we have to set it to 0, 
    # rather than keeping randomly one of the two.
    grid_list[duplicated(grid_list[,c('x','y')]),'z'] = 0
    
    # Remove duplicates again. 
    grid_list = grid_list[!duplicated(grid_list[,c('x','y')], fromLast=T),]
    
    # EXPERIMENTAL PART OVER
    
    # Convert the grid to a matrix. This is the data we will eventually
    # be working with. 
    # Previous comment: "This can come with data loss". This refers to empty
    # gridlines/columns.
    ## [Feb 15, commenting this out because unused.]
    ##grid_matrix = reshape2::dcast(grid_list, y ~ x, fill=0); grid_matrix$y = NULL
    
    #df_mat2 = df_mat[(10**(apply(df_mat, 1, function(x) max(abs(x)))) == 0),]
    #grid_matrix = df_mat2[, (10**(apply(df_mat2, 2, function(x) max(abs(x)))) == 0)]
    
    # [Feb 15, commenting this out because unused.]
    ##colnames(grid_matrix)  = as.character(1:dim(grid_matrix)[2])
    ##row.names(grid_matrix) = as.character(1:dim(grid_matrix)[1])
    
    #grid_list = reshape2::melt(as.matrix(grid_matrix), value.name = 'value')
    #colnames(grid_list) = c('y','x','z')
    # grid_list = grid_list[grid_list$z != 0,]
    factor = 0
    if (realplot) {
      plot = ggplot2::ggplot() + ggplot2::geom_segment(data = paf,
                                                       ggplot2::aes(
                                                         x = tstart,
                                                         xend = tend,
                                                         y = qstart,
                                                         yend = qend
                                                       )) +
        ggplot2::geom_hline(ggplot2::aes(yintercept = gridlines.y + runif(length(gridlines.y),-4000*factor,4000*factor))
                            , color =
                              'grey') +
        ggplot2::geom_vline(ggplot2::aes(xintercept = gridlines.x + runif(length(gridlines.x),-4000*factor,4000*factor))
                            , color =
                              'grey') +
        ggplot2::coord_fixed(
          ratio = 1,
          xlim = NULL,
          ylim = NULL,
          expand = TRUE,
          clip = "on"
        ) +
        ggplot2::theme_bw()
      print(plot)
    }
    
    if (bitlocusplot) {
      
      p = plot_matrix_ggplot(grid_list)
      
      if (saveplot == F){
        print(p)
      } else {
        ggplot2::ggsave(p, file=saveplot, height = 20, width = 20, units = 'cm', device='pdf')
      }
    }
    
    return(list(gridlines.x, gridlines.y, grid_list))
  }



#' enforce_slope_one
#' @author Wolfram Hoeps
#' @export
enforce_slope_one <- function(df){
  
  df$qlen_aln = abs(df$qend - df$qstart)
  df$tlen_aln = abs(df$tend - df$tstart)
  
  df$squaresize = -(apply(-df[,c('qlen_aln','tlen_aln')], 1, function(x) max(x)))
  
  df$qend = df$qstart + df$squaresize
  df$tend = df$tstart + df$squaresize
  
  df[,c('qlen_aln','tlen_aln','squaresize')] = NULL
  
  return(df)
  
}


#' make_xy_grid
#' @author Wolfram Hoeps
#' @export
make_xy_grid <- function(paf, n_additional_bounces=10){
  
  # Points is 'gridpoints'
  points_overall = data.frame()
  
  # Bounce points a couple of times. For the bouncing algorithm, 
  # see separate description. 
  for (i in 1:dim(paf)[1]){
    
    # Bounce start once
    p1_bounced = bounce_point(paf, c(paf[i,]$tstart, paf[i,]$qlow))
    # Bounce end once
    p2_bounced = bounce_point(paf, c(paf[i,]$tend, paf[i,]$qhigh))
    points_overall = rbind(rbind(points_overall, p1_bounced), p2_bounced)
    
    # Bounce bounced #1
    for (j in 1:dim(p1_bounced)[1]){
      points_overall = rbind(points_overall, bounce_point(paf, as.numeric(p1_bounced[j,]))) 
    }
    # Bounce bounced #2
    for (k in 1:dim(p2_bounced)[1]){
      points_overall = rbind(points_overall, bounce_point(paf, as.numeric(p2_bounced[k,]))) 
    }
    
  }
  
  points_overall = unique((points_overall))
  
  # # bounce some more
  for (i in (1:n_additional_bounces)){

    
    points_overall = rbind(points_overall, bounce_point(paf, as.numeric(points_overall[i,])))
    points_overall = unique((points_overall))
    
  }
  
  gridlines.y = unique(round(sort(c(0,points_overall$y))))
  gridlines.x = unique(round(sort(c(0,points_overall$x))))
  
  return(list(gridlines.x, gridlines.y))
}



#ggplot2::ggplot(as.data.frame(grid)) + ggplot2::geom_tile()

