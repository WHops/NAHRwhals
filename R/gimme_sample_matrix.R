#' gimme_sample_matrix
#'
#' @description Think of this as a function used for development and debugging.
#' @param mode diff/same: should x and y be the same sequence?
#' @return matrix (bitlocus)
#'
#' @author Wolfram Höps
#' @export
#' gimme_sample_matrix
#'
#' @description Think of this as a function used for development and debugging.
#' @param mode diff/same: should x and y be the same sequence?
#' @return matrix (bitlocus)
#'
#' @author Wolfram Höps
#' @export
gimme_sample_matrix <- function(mode = 'diff') {
  #samplefasta_link = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/vignettes/simulated_seq_10kb_4SDs.fa'
  #samplemutfasta_link = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/vignettes/simulated_seq_10kb_del_trim.fa'
  # samplefasta_link = system.file('extdata', 'simulated_seq_10kb_4SDs.fa', package =
  #                                  'nahrtoolkit')
  # samplemutfasta_link = system.file('extdata', 'simulated_seq_10kb_dup.fa', package =
  #                                     'nahrtoolkit')
  # outmutfasta2 = './'
  samplefasta_link = system.file('extdata', 'simulated_seq_twonest.fa', package =
                                   'nahrtoolkit')
  samplemutfasta_link = system.file('extdata', 'simulated_seq_twonest_inv_dup.fa', package =
                                      'nahrtoolkit')
  #samplemutfasta_link = '/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/vignettes/simulated_seq_10kb_del.fa'
  
  if (mode == 'same') {
    samplemutfasta_link = samplefasta_link
  }
  
  samplepaf_link = paste0('blub33.paf')
  make_chunked_minimap_alnment(
    samplefasta_link,
    samplemutfasta_link,
    samplepaf_link,
    chunklen = 500,
    minsdlen = 0,
    saveplot = F,
    hllink = F,
    hltype = F,
    quadrantsize = 10000
  )
  grid = wrapper_paf_to_bitlocus(samplepaf_link)[[3]]
  
  gridmatrix = gridlist_to_gridmatrix(grid)
  
  return(sample)
} 
