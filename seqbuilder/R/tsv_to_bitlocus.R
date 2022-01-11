library(nahrtoolkit)

# Whoeps, 8th Jan 2021



find_x_intersection <- function(vectors, point){
  
  # Find vectors that overlap with the point along y axis
  overlap_vecs = vectors[(vectors$qstart < point[2]) & (vectors$qend > point[2]),]
  if (dim(overlap_vecs)[1] == 0){
    # Early return if there is no overlap
    return(point[1])
  }
  
  # If there is overlap, calculate overlap of point with every overlapping line
  x_overlap = c(point[1])
  for (i in 1:dim(overlap_vecs)[1]){
    
    vec = overlap_vecs[i,]
    x_i = (point[2] - vec$y_intercept) / vec$slope
    x_overlap = c(x_overlap, x_i)
  }
  
  return(x_overlap)
}

find_y_intersection <- function(vectors, point){


  overlap_vecs = vectors[(vectors$tstart < point[1]) & (vectors$tend > point[1]),]
  if (dim(overlap_vecs)[1] == 0){
    # Early return if there is no overlap
    return(point[2])
  }
  
  # If there is overlap, calculate overlap of point with every overlapper
  y_overlap = c(point[2])
  for (i in 1:dim(overlap_vecs)[1]){
    
    vec = overlap_vecs[i,]
    y_i = (point[1] - vec$x_intercept) / vec$slope_inv
    y_overlap = c(y_overlap, y_i)
  }
  
  return(y_overlap)
}

add_slope_intercept_info <- function(vector_df, 
                                     xstart = 'tstart',
                                     xend = 'tend', 
                                     ystart = 'qstart',
                                     yend = 'qend'){
  
  dx = vector_df[xend] - vector_df[xstart]
  dy = vector_df[yend] - vector_df[ystart]
  
  # Slope
  m = dy/dx
  
  # m_inv explicitly
  m_inv = 1/m
  
  # y = mx + c    -->    c = y - mx
  y_intercept = vector_df[ystart] - (m * vector_df[xstart])
  
  # x intercept
  x_intercept = y_intercept / m
  
  slope_df = data.frame(m, m_inv, y_intercept, x_intercept)
  colnames(slope_df) = c('slope', 'slope_inv', 'y_intercept', 'x_intercept')
  
  return(slope_df)
}

get_gridlines_x <- function(paf, gp = 10){
  # Let's do it slow and bad for now. This doesn't seem to be very time critical anyway
  # because it only has to be run once per locus. 
  gridlines_x = c()
  for (i in 1:dim(paf)){
    gridlines_x = c(gridlines_x, find_x_intersection(paf, as.numeric(paf[i,c('tstart', 'qstart')])))
    gridlines_x = c(gridlines_x, find_x_intersection(paf, as.numeric(paf[i,c('tend', 'qend')])))
  }
  gridlines_x = sort(unique(gridlines_x))
  
  gridlines_x = gridlines_x[diff(gridlines_x) > gp]
  return(as.integer(gridlines_x[gridlines_x > 0]))
}

get_gridlines_y <- function(paf, gp = 10){
  gridlines_y = c()
  for (i in 1:dim(paf)){
    gridlines_y = c(gridlines_y, find_y_intersection(paf, as.numeric(paf[i,c('tstart', 'qstart')])))
    gridlines_y = c(gridlines_y, find_y_intersection(paf, as.numeric(paf[i,c('tend', 'qend')])))
  }
  gridlines_y = sort(unique(gridlines_y))
  gridlines_y = gridlines_y[diff(gridlines_y) > gp]
  
  return(as.integer(gridlines_y[gridlines_y > 0]))
}

wrapper_paf_to_bitlocus <- function(inpaf, realplot = T, bitlocusplot = T, minlen=2000, gp=10){
  
  # Read paf
  paf = read.table(inpaf)
  colnames(paf) = c('qname','qlen','qstart','qend',
                    'strand','tname','tlen','tstart',
                    'tend','nmatch','alen','mapq')
  

  if (dim(paf)[1] > 50){
    print(dim(paf))
    paf = paf[paf$alen > minlen,]
    print(dim(paf))
  }
  # For negative alignments, interchange start and end. So that
  # 
  paf = transform(paf, 
                  tend = ifelse(strand == '-', tstart, tend),
                  tstart = ifelse(strand == '-', tend, tstart))
  
  paf = cbind(paf, add_slope_intercept_info(paf))
  
  gridlines_x = get_gridlines_x(paf, gp = gp)
  gridlines_y = get_gridlines_y(paf, gp = gp)
  
  # Grid fill machinery
  grid = matrix(0, length(gridlines_x), length(gridlines_y))
  # Fill the bitlocus. Naive approach.
  for (x in 1:length(gridlines_x)-1){
    x_mid = (gridlines_x[x] + gridlines_x[x+1]) / 2
    
    for (y in 1:length(gridlines_y)-1){
      y_mid = (gridlines_y[y] + gridlines_y[y+1]) / 2
      
      grid[x, y] = get_aln_overlap_in_sector(paf, 
                                             gridlines_x[x], 
                                             gridlines_y[y], 
                                             x_mid, 
                                             y_mid, 
                                             gridlines_x[x+1], 
                                             gridlines_y[y+1], 
                                             gp = gp)
    }
  }
  
  if (realplot){
    plot = ggplot2::ggplot() + ggplot2::geom_segment(data=paf, 
      ggplot2::aes(x=tstart, xend=tend, y=qstart, yend=qend)) +
      ggplot2::geom_hline(ggplot2::aes(yintercept=gridlines_y)) +
      ggplot2::geom_vline(ggplot2::aes(xintercept=gridlines_x))# +
      #ggplot2::xlim(c(0,10000)) + ggplot2::ylim(c(0,10000))
    print(plot)
  }
  
  if (bitlocusplot){
    image(log2(abs(grid)))
  }
  
  return(grid)
}

get_aln_overlap_in_sector <- function(paf, x_start, y_start, x_mid, y_mid, x_end, y_end, gp = 10){
  
  # We say that a line traverses a gridpart if the mid of the grid
  # is on the line. We only have to make sure (later) that the 
  # line actually extends that long (remember we have vectors, not lines.)
  hits = which( abs(((x_mid * paf$slope) + paf$y_intercept) - y_mid) <= gp )
  
  #x_mid
  #everyones_y = ((x_mid * paf$slope) + paf$y_intercept)
  
  
  # If anything is found, continue
  if (length(hits) > 0){
    hitpaf = paf[hits,]
    if (y_start > 2500){
      browser()
    }
    paf_ystart_at_x = ((x_start * hitpaf$slope) + hitpaf$y_intercept)
    paf_yend_at_x = ((x_end * hitpaf$slope) + hitpaf$y_intercept)
    
    hitpaf = hitpaf[(abs(paf_ystart_at_x - y_start) < gp) &
                    (abs(paf_yend_at_x - y_end) < gp),
          ]
    
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
    
    if (dim(hitpaf)[1] > 0){
      if (dim(hitpaf)[1] != 1){
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





# Turn tsv into bitlocus.
origfa = '../vignettes/simulated_seq_10kb_4SDs.fa'
mutfa = '../vignettes/simulated_seq_10kb_del_trim.fa'
#mutfa = '../vignettes/mut.fa'
outpaf = '../vignettes/bitlocus8.paf'
#samplefasta_link = system.file('extdata', '10ktest.fa', package='nahrtoolkit')

outpaf = '/Library/Frameworks/R.framework/Versions/4.1/Resources/library/nahrtoolkit/extdata/10ktest.fa62232.paf'
outpaf_link = '/Users/hoeps/phd/projects/nahrcall/nahrchainer/seqbuilder/res/outpaf1d'
outpaf_link = '/Users/hoeps/phd/projects/nahrcall/nahrchainer/seqbuilder/res/outpafsds13335'
outpaf_link = paste0(samplefasta_link, '62232122.paf')
inpaf = outpaf_link
grid = wrapper_paf_to_bitlocus(outpaf_link)

# plot = make_chunked_minimap_alnment(origfa, mutfa, outpaf, 
#                              outplot=NULL, chunklen = 1000, minsdlen = 10, saveplot=F, 
#                              hllink = outpaf, hltype = 'paf', quadrantsize = 1000)

# Read paf












#ggplot2::ggplot(as.data.frame(grid)) + ggplot2::geom_tile()




