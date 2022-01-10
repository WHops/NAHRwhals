# Whoeps, 8th Jan 2021
find_x_intersection <- function(vectors, point){
  
  #vectors = paf
  #point = as.numeric(paf[1,c('tend', 'qend')])
  x_overlap = NULL
  #vectors = data.frame(xstart = 2, xend = 4, ystart = 6, yend = 8)
  #point = c(5,7)
  
  overlap_vecs = vectors[(vectors$qstart < point[2]) & (vectors$qend > point[2]),]
  if (dim(overlap_vecs)[1] == 0){
    # Early return if there is no overlap
    return(point[1])
  }
  
  # If there is overlap, calculate overlap of point with every overlapper
  x_overlap = c(point[1])
  for (i in 1:dim(overlap_vecs)[1]){
    vec = overlap_vecs[i,]
    dy = vec$qend - vec$qstart
    dx = vec$tend - vec$tstart
    m = dy/dx
    c = vec$qstart - (m * vec$tstart)
    
    x_i = (point[2] - c) / m
    x_overlap = c(x_overlap, x_i)
  }
  
  return(x_overlap)
}

find_y_intersection <- function(vectors, point){
  
  #vectors = paf
  #point = as.numeric(paf[1,c('tend', 'qend')])
  x_overlap = NULL
  #vectors = data.frame(xstart = 2, xend = 4, ystart = 6, yend = 8)
  #point = c(5,7)
  
  overlap_vecs = vectors[(vectors$tstart < point[1]) & (vectors$tend > point[1]),]
  if (dim(overlap_vecs)[1] == 0){
    # Early return if there is no overlap
    return(point[2])
  }
  
  # If there is overlap, calculate overlap of point with every overlapper
  x_overlap = c(point[2])
  for (i in 1:dim(overlap_vecs)[1]){
    vec = overlap_vecs[i,]
    dy = vec$tend - vec$tstart
    dx = vec$qend - vec$qstart
    m = dy/dx
    c = vec$tstart - (m * vec$qstart)
    
    x_i = (point[1] - c) / m
    x_overlap = c(x_overlap, x_i)
  }
  
  return(x_overlap)
}



# Turn tsv into bitlocus.
origfa = '../vignettes/simulated_seq_10kb_4SDs.fa'
mutfa = '../vignettes/simulated_seq_10kb_del_trim.fa'
#mutfa = '../vignettes/mut.fa'
outpaf = '../vignettes/bitlocus8.paf'
#samplefasta_link = system.file('extdata', '10ktest.fa', package='nahrtoolkit')

plot = make_chunked_minimap_alnment(origfa, mutfa, outpaf, 
                             outplot=NULL, chunklen = 1000, minsdlen = 10, saveplot=F, 
                             hllink = outpaf, hltype = 'paf', quadrantsize = 1000)

# Read paf
paf = read.table(outpaf)
colnames(paf) = c('qname','qlen','qstart','qend',
                 'strand','tname','tlen','tstart',
                 'tend','nmatch','alen','mapq')
paf = transform(paf, 
                tend = ifelse(strand == '-', tstart, tend),
                tstart = ifelse(strand == '-', tend, tstart))

# Let's do it slow and bad for now. This doesn't seem to be very time critical anyway
# because it only has to be run once per locus. 
gridlines_x = c()
for (i in 1:dim(paf)){
  print(i)
  gridlines_x = c(gridlines_x, find_x_intersection(paf, as.numeric(paf[i,c('tstart', 'qstart')])))
  gridlines_x = c(gridlines_x, find_x_intersection(paf, as.numeric(paf[i,c('tend', 'qend')])))
  print(gridlines_x)
}
gridlines_x = sort(unique(gridlines_x))


gridlines_y = c()
for (i in 1:dim(paf)){
  print(i)
  gridlines_y = c(gridlines_y, find_y_intersection(paf, as.numeric(paf[i,c('tstart', 'qstart')])))
  gridlines_y = c(gridlines_y, find_y_intersection(paf, as.numeric(paf[i,c('tend', 'qend')])))
  print(gridlines_y)
}
gridlines_y = sort(unique(gridlines_y))


gridlines = sort(unique(c(gridlines_x, gridlines_y)))

grid = matrix(0, length(gridlines), length(gridlines))
# Fill the bitlocus. Naive approach.
for (x in 1:length(gridlines)){
  for (y in 1:length(gridlines)){
    #if((gridlines[x] %in% paf$tstart) & (gridlines[y] %in% paf$qstart)){
    #print()
    if (any(paf$tstart == gridlines[x] & paf$qstart == gridlines[y])){
      grid[x, y] = 1
    }
    if (any(paf$tstart == gridlines[x+1] & paf$qstart == gridlines[y+1])){
        grid[x, y] = -1  
    }
  }
}






