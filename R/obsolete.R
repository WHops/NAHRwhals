# Collection of obsolete functions probably no longer needed #

#' find_x_intersection

#' @author Wolfram Hoeps
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


#' get_aln_overlap_in_sector

#' @author Wolfram Hoeps
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

