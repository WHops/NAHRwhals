#library(nahrtoolkit)

# Whoeps, 8th Jan 2021


#' find_x_intersection

#' @author Wolfram Hoeps
#' @rdname alignment
#' @export
find_x_intersection <- function(vectors, point) {
  # Find vectors that overlap with the point along y axis
  overlap_vecs = vectors[(vectors$qstart < point[2]) &
                           (vectors$qend > point[2]), ]
  if (dim(overlap_vecs)[1] == 0) {
    # Early return if there is no overlap
    return(point[1])
  }
  
  # If there is overlap, calculate overlap of point with every overlapping line
  x_overlap = c(point[1])
  for (i in 1:dim(overlap_vecs)[1]) {
    vec = overlap_vecs[i, ]
    x_i = (point[2] - vec$y_intercept) / vec$slope
    x_overlap = c(x_overlap, x_i)
  }
  
  return(x_overlap)
}





#' find_y_intersection

#' @author Wolfram Hoeps
#' @rdname alignment
#' @export
find_y_intersection <- function(vectors, point) {
  
  # DEBUG 
  #vectors = paf
  #point = paf[4, c('tstart', 'qstart')]
  # Find vectors that overlap with the point along y axis
  overlap_vecs = vectors[((vectors$tstart < as.numeric(point[1])) &
                           (vectors$tend > as.numeric(point[1]))), ]
  if (dim(overlap_vecs)[1] == 0) {
    # Early return if there is no overlap
    return(point[2])
  }
  
  # If there is overlap, calculate overlap of point with every overlapper
  y_overlap = c(point[2])
  for (i in 1:dim(overlap_vecs)[1]) {
    vec = overlap_vecs[i, ]
    y_i = ((vec$slope) * as.numeric(point[1])) + vec$y_intercept
    #y_i = (point[1] + vec$x_intercept) / vec$slope_inv
    y_overlap = c(y_overlap, y_i)
  }
  
  return(y_overlap)
}

#' add_slope_intercept_info

#' @author Wolfram Hoeps
#' @rdname alignment
#' @export
add_slope_intercept_info <- function(vector_df,
                                     xstart = 'tstart',
                                     xend = 'tend',
                                     ystart = 'qstart',
                                     yend = 'qend') {
  dx = vector_df[xend] - vector_df[xstart]
  dy = vector_df[yend] - vector_df[ystart]
  
  # Slope
  m = dy / dx
  
  # m_inv explicitly
  m_inv = 1 / m
  
  # y = mx + c    -->    c = y - mx
  y_intercept = vector_df[ystart] - (m * vector_df[xstart])
  
  # x intercept
  x_intercept = y_intercept / m
  
  slope_df = data.frame(m, m_inv, y_intercept, x_intercept)
  colnames(slope_df) = c('slope', 'slope_inv', 'y_intercept', 'x_intercept')
  
  return(slope_df)
}





#' get_gridlines.x

#' @author Wolfram Hoeps
#' @rdname alignment
#' @export
get_gridlines.x <- function(paf, gp = 10) {
  # Let's do it slow and bad for now. This doesn't seem to be very time critical anyway
  # because it only has to be run once per locus.
  gridlines.x = c()
  for (i in 1:dim(paf)[1]) {
    gridlines.x = c(gridlines.x, find_x_intersection(paf, as.numeric(paf[i, c('tstart', 'qstart')])))
    gridlines.x = c(gridlines.x, find_x_intersection(paf, as.numeric(paf[i, c('tend', 'qend')])))
  }
  gridlines.x = sort(unique(gridlines.x[gridlines.x >= 0]))
  gridlines.x = gridlines.x[diff(gridlines.x) > gp]
  
  return(as.integer(gridlines.x))
}


#' get_gridlines.y

#' @author Wolfram Hoeps
#' @rdname alignment
#' @export
get_gridlines.y <- function(paf, gp = 10) {
  gridlines.y = c()
  for (i in 1:dim(paf)[1]) {
    gridlines.y = c(gridlines.y, 
                    find_y_intersection(paf, as.numeric(paf[i, c('tstart', 'qstart')])))
    gridlines.y = c(gridlines.y, 
                    find_y_intersection(paf, as.numeric(paf[i, c('tend', 'qend')])))
  }
  
  gridlines.y = sort(unique(gridlines.y[gridlines.y >= 0]))
  gridlines.y = gridlines.y[c(T, diff(gridlines.y) > gp)]
  
  return(as.integer(gridlines.y))
}


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
#' @rdname alignment
#' @export
wrapper_paf_to_bitlocus <-
  function(inpaf,
           realplot = T,
           bitlocusplot = T,
           minlen = 1000,
           gp = 10,
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
    

    paf = paf[paf$alen > minlen, ]
    
    paf[,c('qstart','qend','tstart','tend')] = round(data.frame(paf[,c('qstart','qend','tstart','tend')] / compression),0) * compression
    #paf = paf[paf$qend > 100000,]
    #  print(dim(paf))
    #}
    
    paf = enforce_slope_one(paf)
    
    
    #paf[5,]$tstart = paf[5,]$tstart + 10000
    #paf[5,]$tend = paf[5,]$tend + 10000
    #paf = paf[c(1,2,5),]
    #paf[2,'tstart'] = paf[2,'tstart'] - 10000
    #paf[2,'qstart'] = paf[2,'qstart'] - 10000
    
    #paf = paf[c(1:3, 5),]
    # For negative alignments, interchange start and end. So that
    #
    
    
    paf = transform(
      paf,
      qend = ifelse(strand == '-', qstart, qend),
      qstart = ifelse(strand == '-', qend, qstart)
    )
    
    paf = cbind(paf, add_slope_intercept_info(paf))
    
    #gp = 1
    
    paf$qlow = -(apply(-paf[,c('qstart','qend')], 1, function(x) max(x)))
    paf$qhigh = apply(paf[,c('qstart','qend')], 1, function(x) max(x))
    
    points_overall = data.frame()
    
    for (i in 1:dim(paf)[1]){
      p1_bounced = bounce_point(paf, c(paf[i,]$tstart, paf[i,]$qlow))
      p2_bounced = bounce_point(paf, c(paf[i,]$tend, paf[i,]$qhigh))
      points_overall = rbind(rbind(points_overall, p1_bounced), p2_bounced)

        for (j in 1:dim(p1_bounced)[1]){
          points_overall = rbind(points_overall, bounce_point(paf, as.numeric(p1_bounced[j,]))) 
        }
        for (k in 1:dim(p2_bounced)[1]){
          points_overall = rbind(points_overall, bounce_point(paf, as.numeric(p2_bounced[k,]))) 
        }

    }

    points_overall = unique((points_overall))
     
    # # bounce some more
    # for (i in (1:10)){
    #   
    #   #print(gridlines.y)
    #  # print(gridlines.x)
    #   
    #   points_overall = rbind(points_overall, bounce_point(paf, as.numeric(points_overall[i,]))) 
    #   points_overall = unique((points_overall))
    #   #print(points_overall)
    #   #print(dim(points_overall))
    #   
    #   
    #   gridlines.y = unique(round(sort(c(0,points_overall$y))))
    #   gridlines.x = unique(round(sort(c(0,points_overall$x))))
    #   
    # 
    # }
    # 
    gridlines.y = unique(round(sort(c(0,points_overall$y))))
    gridlines.x = unique(round(sort(c(0,points_overall$x))))

    print(gridlines.y)
    print(gridlines.x)
    #gridlines.y == gridlines.x
    # paf = transform(
    #   paf,
    #   qend = ifelse(strand == '-', qstart, qend),
    #   qstart = ifelse(strand == '-', qend, qstart)
    # )
    
    df = data.frame()
    for (i in 1:dim(paf)[1]) {
      df = rbind(df, as.data.frame(
        bresenham(
          x = as.numeric(paf[i, c('tstart', 'tend')]) ,
          y = as.numeric(paf[i, c('qstart', 'qend')]) ,
          gridlines.x,
          gridlines.y,
          debug = F
        )
      ))
    }
    
    #gp=100
    # This here can come with loss of data. ..
    df_mat = reshape2::dcast(df, y ~ x, fill=0); df_mat$y = NULL
    df_mat2 = df_mat[(10**(apply(df_mat, 1, function(x) max(abs(x)))) > gp) |
                     (10**(apply(df_mat, 1, function(x) max(abs(x)))) == 0),]
    
    df_mat3 = df_mat2[,(10**(apply(df_mat2, 2, function(x) max(abs(x)))) > gp) |
                       (10**(apply(df_mat2, 2, function(x) max(abs(x)))) == 0)]
    
    colnames(df_mat3) = as.character(1:dim(df_mat3)[2])
    row.names(df_mat3) = as.character(1:dim(df_mat3)[1])
    
    df_back = reshape2::melt(as.matrix(df_mat3), value.name = 'value')
    colnames(df_back) = c('y','x','z')
    df_back = df_back[df_back$z != 0,]
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
      p = ggplot2::ggplot(df_back) + ggplot2::geom_tile(ggplot2::aes(
        x = x,
        y = y,
        fill = sign(z) * log10(abs(z))
      )) +
        ggplot2::scale_fill_gradient2(low = 'red',
                                      mid = 'black',
                                      high = 'blue') +
        ggplot2::coord_fixed(
          ratio = 1,
          xlim = NULL,
          ylim = NULL,
          expand = TRUE,
          clip = "on"
        ) +
        ggplot2::theme_bw()
      print(p)
    }
    
    return(df_back)
  }

#' get_aln_overlap_in_sector

#' @author Wolfram Hoeps
#' @rdname alignment
#' @export
get_aln_overlap_in_sector <-
  function(paf,
           x_start,
           y_start,
           x_mid,
           y_mid,
           x_end,
           y_end,
           gp = 10) {
    # We say that a line traverses a gridpart if the mid of the grid
    # is on the line. We only have to make sure (later) that the
    # line actually extends that long (remember we have vectors, not lines.)
    hits = which(abs(((x_mid * paf$slope) + paf$y_intercept) - y_mid) <= gp)
    
    #x_mid
    #everyones_y = ((x_mid * paf$slope) + paf$y_intercept)
    
    
    # If anything is found, continue
    if (length(hits) > 0) {
      hitpaf = paf[hits, ]
      if (y_start > 2500) {
        browser()
      }
      paf_ystart_at_x = ((x_start * hitpaf$slope) + hitpaf$y_intercept)
      paf_yend_at_x = ((x_end * hitpaf$slope) + hitpaf$y_intercept)
      
      hitpaf = hitpaf[(abs(paf_ystart_at_x - y_start) < gp) &
                        (abs(paf_yend_at_x - y_end) < gp),]
      
      # # Criterion plus: 'Positive alignment'
      # crit_plus = hitpaf$slope > 0
      #
      # # Criterion A: 'Line Transitioning gridppoint':
      # # The gridpoint lies 'strictly' between start and end of the line.
      # crit_A = ((hitpaf$tstart - gp < x_start) &
      #             (hitpaf$tend  + gp > x_start))
      # # Criterion B: 'Line starting at gridpoint':
      # # This line starts in this quadrant.
      # crit_B = abs(hitpaf$tstart - x_start) <= gp
      #
      # # Criterion C: 'Line starting at gridpoint but inverted
      # crit_C = abs(hitpaf$tend - x_start) <= gp
      
      #hitpaf = hitpaf[(crit_A | (crit_plus & crit_B) | ((!crit_plus) & crit_C)),]
      
      if (dim(hitpaf)[1] > 0) {
        if (dim(hitpaf)[1] != 1) {
          browser()
        }
        stopifnot(dim(hitpaf)[1] == 1)
        return(abs(x_start - x_mid) * 2 * sign(hitpaf[1, 'slope']))
      } else {
        return(0)
      }
    } else {
      return(0)
    }
  }

#grid = wrapper_paf_to_bitlocus(samplepaf_link, gp = 0)



# # Turn tsv into bitlocus.
origfa = '../vignettes/simulated_seq_10kb_4SDs.fa'
mutfa = '../vignettes/simulated_seq_10kb_dup.fa'
# #mutfa = '../vignettes/mut.fa'
outpaf = as.character(runif(1,1e10,1e11))#'../vignettes/bitlocus8.paf'

# #samplefasta_link = system.file('extdata', '10ktest.fa', package='nahrtoolkit')
#
# #outpaf_link = '/Library/Frameworks/R.framework/Versions/4.1/Resources/library/nahrtoolkit/extdata/10ktest.fa62232.paf'
# outpaf_link = '/Users/hoeps/phd/projects/nahrcall/nahrchainer/seqbuilder/res/outpaf1d'
# # outpaf_link = '/Users/hoeps/phd/projects/nahrcall/nahrchainer/seqbuilder/res/outpafsds13335'
# # outpaf_link = paste0(samplefasta_link, '62232122.paf')
# # inpaf = outpaf_link

mutfa_end = '../vignettes/simulated_seq_10kb_dup_end.fa'

make_chunked_minimap_alnment(origfa, mutfa, outpaf, 
                              chunklen = 500, minsdlen = 500, saveplot=F,
                              hllink = outpaf, hltype = 'paf', quadrantsize = 1000)
grid = wrapper_paf_to_bitlocus(outpaf, minlen = 100, gp = 10)

# Read paf





enforce_slope_one <- function(df){
  
  df$qlen_aln = abs(df$qend - df$qstart)
  df$tlen_aln = abs(df$tend - df$tstart)
  
  df$squaresize = -(apply(-df[,c('qlen_aln','tlen_aln')], 1, function(x) max(x)))
  
  df$qend = df$qstart + df$squaresize
  df$tend = df$tstart + df$squaresize
  
  df[,c('qlen_aln','tlen_aln','squaresize')] = NULL
  
  return(df)
  
}





#ggplot2::ggplot(as.data.frame(grid)) + ggplot2::geom_tile()

