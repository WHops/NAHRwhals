rep.row <- function(x, n) {
  matrix(rep(x, each = n), nrow = n)
}
rep.col <- function(x, n) {
  matrix(rep(x, each = n), ncol = n, byrow = TRUE)
}




#' calc_coarse_aln_score
#'
#' A function assigning an alignment score to a condensed dotplot.
#' we are using a base-function from geeksforgeeks.org, which is
#' essentially like needleman-wunsch.
#' Since in our case the costs for traversing a node 'vertically' or
#' 'horizontally' are different, we had to work around this a bit by
#' adapting the algorithm.
#'
#' @param mat matrix. Bitlocus / Gridmatrix.
#' @author  geeksforgeeks.org
#' @export
calc_coarse_grained_aln_score <- function(mat, verbose = F, forcecalc = F) {
  # Save matrix dimensions.
  dim_ = dim(mat)
  row = dim_[1]
  col = dim_[2]
  
  # Decide if it is 'worth' to calculate an aln score, of if it's
  # clear that this result will not be good.
  symmetry = min(dim_) / max(dim_)
  
  if ((symmetry < 0.9) & (forcecalc == F)) {
    return(NA)
  }
  
  # Construct our three cost matrices:
  # cost_u, if you want to go 'up'
  # cost_r, if you want to go 'right'
  # cord_d, if you want to go 'diagonal'
  climb_up_cost = colMax(as.data.frame(t(abs(mat))))
  walk_right_cost = colMax(as.data.frame(abs(mat)))
  cost_u = rep.col(climb_up_cost, dim(mat)[2])
  cost_r = rep.row(walk_right_cost, dim(mat)[1])
  
  # cost_d: Diagonal is only allowed if there is a pos alignment.
  #  0 if we have a positive alignment
  #  Inf otherwise.
  cost_d = -mat
  cost_d[cost_d >= 0] = Inf
  cost_d[cost_d < 0] = 0
  
  # ... and initialize a results matrix.
  cost_res = cost_u - cost_u
  
  
  ### Fill the cost matrix ###
  
  # For 1st column
  for (i in 2:row) {
    cost_res[i, 1] = cost_u[i, 1] + cost_res[i - 1, 1]
  }
  # For 1st row
  for (j in 2:col) {
    cost_res[1, j] = cost_r[1, j] + cost_res[1, j - 1]
  }
  # For rest of the 2d matrix
  for (i in 2:row) {
    for (j in 2:col) {
      cost_res[i, j] =  (min(
        cost_d[i - 1, j - 1] + cost_res[i - 1, j - 1],
        cost_res[i - 1, j] + cost_u[i - 1, j],
        cost_res[i, j - 1] + cost_r[i, j - 1]
      ))
    }
  }

  if (verbose) {
    print(cost_res)
    print(paste0(
      'Unmatching: ca ',
      cost_res[row, col],
      ' out of ',
      sum(climb_up_cost, walk_right_cost),
      ' (',
      round((100 * cost_res[row, col]) / sum(climb_up_cost, walk_right_cost), 2),
      '%)'
    ))
    
  }
  
  # Return value: percentage of basepairs not crossed along alignments. 
  return((cost_res[row, col] * 100) / sum(climb_up_cost, walk_right_cost))
}





#' find_maxdiag
#'
#' Scan a matrix for the length of the longest (positive) diagonal. This is part
#' of the 'eval' method.
#' @param m A matrix
#' @author Unknown (https://stackoverflow.com/questions/47330812/find-the-longest-diagonal-of-an-element-in-a-matrix-python)
#' @export
find_maxdiag <- function(m) {
  # Pre-compute dimension, because we need this often.
  dim_m = dim(m)
  dim_m_1 = dim_m[1]
  dim_m_2 = dim_m[2]
  
  maxdiag = 0
  for (i in 1:(dim_m_1 + 1)) {
    # Iterate x points
    for (j in 1:(dim_m_2 + 1)) {
      # Iterate y points
      k = 0 # k = current diag size
      
      while ((i + k < dim_m_1) &
             # make sure i+k is not out of bounds
             (j + k < dim_m_2)) {
        # make sure j+k is not out of bonds
        direction = sign(m[i, j])
        
        if ((direction != 0) &
            (sign(m[i + k, j + k]) == direction)) {
          # that next point is not zero
          k = k + 1 # Diag has become longer
          if (k > maxdiag) {
            maxdiag = k
          }
        } else {
          k = k + 1e10
        }
        
      }
    }
  }
  return(maxdiag + 1)
}



#' eval_mutated_seq
#'
#' @description Evaluation function: how similar are the input and output?
#' This needs work still...
#'
#' @param bitlocus  an nxm matrix (n=gridpoints_x, m=gridpoints_y)
#' @return numeric [0-1] with evaluation score.
#'
#' @author Wolfram Höps
#' @export
eval_mutated_seq <- function(bitlocus) {
  # Old stuff here:
  #coverage_x =  sum(colMax(as.data.frame(bitlocus)) > 0) / dim(bitlocus)[2]
  #coverage_y =  sum(colMax(as.data.frame(t(bitlocus))) > 0) / dim(bitlocus)[1]
  
  
  # A matrix is perfect if it has a non-zero, non-negative diagonal.
  # We also want the matrix to be symmetrical.
  
  
  symmetry = min(dim(bitlocus)) / max(dim(bitlocus))
  # diag_filled = sum(diag(bitlocus) > 0) / min(dim(bitlocus))
  # minusdiag_filled =
  
  # Implement properly here a function defining if the overall alignment is pos
  # or negative.
  if (symmetry > 0.9) {
    if (sign(bitlocus[1, 1]) == 1) {
      maxdiag = find_maxdiag(bitlocus)
    } else if (sign(bitlocus[1, 1]) == -1) {
      maxdiag = find_maxdiag(apply(bitlocus, 2, rev))
    } else{
      maxdiag = max(find_maxdiag(apply(bitlocus, 2, rev)), find_maxdiag(bitlocus))
    }
  } else {
    maxdiag = 0
  }
  #return(maxdiag)
  return(round(maxdiag * symmetry, 3))
  #  return(round(maxdiag/dim(bitlocus)[2], 3))
}


#' add_eval
#' CRITICAL: NEEDS IMPROVEMENTS
#' @description A helperfunction, wrapping the evaluation function.
#' @param bitlocus matrix, nxm
#' @param depth How many SVs in sequence should be simulated?
#' @return evaluation matrix
#'
#' @author Wolfram Höps
#' @export
add_eval <-
  function(res,
           m,
           layer,
           pair_level1,
           pair_level2,
           pair_level3) {
    
    # Ich kotz.
    
    if (layer == 1) {
      res_add = unlist(c(
        calc_coarse_grained_aln_score(m),
        paste(pair_level1$p1, pair_level1$p2, pair_level1$sv, sep = '_')
      ))
      
      
    } else if (layer == 2) {
      res_add = unlist(c(
        calc_coarse_grained_aln_score(m),
        paste(pair_level1$p1, pair_level1$p2, pair_level1$sv, sep = '_'),
        paste(pair_level2$p1, pair_level2$p2, pair_level2$sv, sep =
                '_')
      ))
      
    } else if (layer == 3) {
      # add
      res_add = unlist(c(
        calc_coarse_grained_aln_score(m),
        paste(pair_level1$p1,
              pair_level1$p2,
              pair_level1$sv,
              sep = '_'),
        paste(pair_level2$p1,
              pair_level2$p2,
              pair_level2$sv,
              sep = '_'),
        paste(pair_level3$p1,
              pair_level3$p2,
              pair_level3$sv,
              sep = '_')
      ))
      
    }
    return(res_add)
  }
