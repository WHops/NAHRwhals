
#' rep.row
#' 
#' Tiny helperfunction to repeat rows in a matrix
#' @author  Wolfram Hoeps
#' @export
rep.row <- function(x, n) {
  matrix(rep(x, each = n), nrow = n)
}

#' rep.col
#' 
#' Tiny helperfunction to repeat columns in a matrix
#' @author  Wolfram Hoeps
#' @export
rep.col <- function(x, n) {
  matrix(rep(x, each = n), ncol = n, byrow = TRUE)
}

#' absmin
#' 
#' Tiny helperfuction to return the absolute minimum
#' @author  Wolfram Hoeps
#' @export
absmin <- function(x) { x[which.min( abs(x) )]}

#' calculate_estimated_aln_score
#' 
#' Quick & Dirty estimation of the score of a condensed dotplot. 
#' Is used to filter out hopeless cases. 
#' @param mat A matrix object representing a condensed dotplot
#' @output a percentage estimate for the score
#' @author  Wolfram Hoeps
#' @export
calculate_estimated_aln_score <- function(mat){
  
  matdiag = diagonalize(mat, 0.25)
  sum_of_rows_with_pos_aln_bp = sum(matrixStats::rowMaxs(matdiag))
  alt_seq_len_bp = sum(as.numeric(rownames(mat)))
  
  sum_of_cols_with_pos_aln_bp = sum(matrixStats::colMaxs(matdiag))
  ref_seq_len_bp = sum(as.numeric(colnames(mat)))
  
  pct_covered = ((sum_of_rows_with_pos_aln_bp + sum_of_cols_with_pos_aln_bp) / 
                   (alt_seq_len_bp + ref_seq_len_bp)) * 100 
  
  return(pct_covered) # This replaces the previous version (next line)
  #return(sum(na.omit(rowMeans(replace(matdiag, matdiag == 0, NA), na.rm = TRUE))))
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
#' @author  geeksforgeeks.org; Wolfram Hoeps
#' @export
calc_coarse_grained_aln_score <-
  function(mat,
           old_way_of_calc = F,
           verbose = F,
           forcecalc = F,
           orig_symm = 1) {
    
    # Remove zero-pads.
    mat = matrix_remove_zero_pads(mat)
    
    n_eval_total <<- n_eval_total + 1
    
    # If matrix has no entries (no alignments), return 0
    # If matrix is only one number, report 100%.
    if (is.null(dim(mat))) {
      return(0)
    } else if ((dim(mat)[1] == 1) & (dim(mat)[2] == 1)) {
      return(100)
    } else if ((dim(mat)[1] == 1) & (dim(mat)[2] > 1)) {
      pos_aln = sum(mat[mat>0])
      return(round((pos_aln / sum(
        as.numeric(colnames(mat))
      )) * 100, 3))
    } else if ((dim(mat)[1] > 1) & (dim(mat)[2] == 1)) {
      pos_aln = sum(mat[mat>0])
      return(round((pos_aln / sum(
        as.numeric(row.names(mat))
      )) * 100, 3))
    }
    # Save matrix dimensions.
    dim_ = dim(mat)
    row = dim_[1]
    col = dim_[2]
    # Decide if it is 'worth' to calculate an aln score, of if it's
    # clear that this result will not be good.
    climb_up_cost = as.numeric(row.names(mat))
    walk_right_cost = as.numeric(colnames(mat))
    symmetry = min(sum(climb_up_cost), sum(walk_right_cost)) / max(sum(climb_up_cost), sum(walk_right_cost))
    
    
    
    
    symm_factor = 1
    # if (is.na(symmetry)){
    #   browser()
    # }
    # Run away if there are at least 5 columns, and we have less than 95% similarity
    if ((symmetry < (orig_symm*symm_factor) & row > 5) & (forcecalc == F)) {
      return(1)
    } else if ((dim(mat)[1] == 1) & (dim(mat)[2] == 1)) {
      return(100)
    } else if ((dim(mat)[1] == 1) & (dim(mat)[2] > 1)) {
      pos_aln = sum(mat)
      return(round((pos_aln / sum(
        as.numeric(row.names(mat))
      )) * 100, 3))
    } else if ((dim(mat)[1] > 1) & (dim(mat)[2] == 1)) {
      pos_aln = sum(mat)
      return(round((pos_aln / sum(
        as.numeric(colnames(mat))
      )) * 100, 3))
    }
    

    # If all this, then one more chance is to throw stuff

    est = calculate_estimated_aln_score(mat)
    if (forcecalc){
      est_highest <<- est
    } 
    
    if (est < (est_highest*0.9)){
      return(1)#est)
      #print('hi')
    } else if (est >= (est_highest*0.9)){
      est_highest <<- est
    }
    
    n_eval_calc <<- n_eval_calc+1
    
    # Construct our three cost matrices:
    # cost_u, if you want to go 'up'
    # cost_r, if you want to go 'right'
    # cord_d, if you want to go 'diagonal'

    climb_up_cost = as.numeric(row.names(mat))
    walk_right_cost = as.numeric(colnames(mat))
    
    cost_u = rep.col(climb_up_cost, dim(mat)[2])
    cost_r = rep.row(walk_right_cost, dim(mat)[1])
    
    cumsum_u = cumsum(climb_up_cost)
    cumsum_r = cumsum(walk_right_cost)
    # cost_d: Diagonal is only allowed if there is a pos alignment.
    #  0 if we have a positive alignment
    #  Inf otherwise.
    cost_d = -mat
    cost_d[cost_d >= 0] = Inf
    cost_d[cost_d < 0] = 0
    
    # ... and initialize a results matrix.
    cost_res = cost_u - cost_u
    
    
    ### Fill the cost matrix ###
    if (dim_[1] < 10) {
      calc_corridor_pct = 1
    } else if (dim_[1] < 30) {
      calc_corridor_pct = 0.3
    } else {
      calc_corridor_pct = 0.2
    }
    
    if (forcecalc) {
      calc_corridor_pct = 1
    }

    # 1st value
    cost_res[1, 1] = min(climb_up_cost[1] + walk_right_cost[1], cost_d[1, 1])
    
    
    # For 1st column
    for (i in 2:row) {
      if (i / row > calc_corridor_pct) {
        cost_res[i, 1] = Inf
      } else {
        cost_res[i, 1] = min(cost_u[i, 1] + cost_res[i - 1, 1],
                             cumsum_u[i - 1] + cost_d[i, 1])
        
      }
    }
    # For 1st row
    for (j in 2:col) {
      if (j / col > calc_corridor_pct) {
        cost_res[1, j] = Inf
      } else {
        cost_res[1, j] = min(cost_r[1, j] + cost_res[1, j - 1],
                             cumsum_r[j - 1] + cost_d[1, j])
      }
    }

    #### Out custom speed-up #####
    cost_res[2:row, 2:col] = Inf


    i_mat = rep.col(  ((1:row) / row), col)
    j_mat = rep.row(  ((1:col) / col), row)

    middleval_mat = 1 - ((i_mat > (j_mat + calc_corridor_pct)) + 
                    (j_mat > (i_mat + calc_corridor_pct)))
    middleval_mat[1,] = 0
    middleval_mat[,1] = 0
    
    idxs = as.vector(t(which(t(middleval_mat) == 1, arr.ind = T)))

    for (idx in seq(1,length(idxs),2)){
      i = idxs[idx+1]
      j = idxs[idx]
      cost_res[i, j] =  (min(                                
        cost_res[i - 1, j - 1] + cost_d[i, j],
        cost_res[i - 1, j] +     cost_u[i, j],
        cost_res[i, j - 1] +     cost_r[i, j]
      ))
    }
    ###### OVER ###########
    
    if (cost_res[row, col] == Inf) {
      return(1)
    }
    # For rest of the 2d matrix
    # for (i in 2:row) {
    #   for (j in 2:col) {
    #     cost_res[i, j] =  (min(
    #       cost_d[i - 1, j - 1] + cost_res[i - 1, j - 1],
    #       cost_res[i - 1, j] + cost_u[i - 1, j],
    #       cost_res[i, j - 1] + cost_r[i, j - 1]
    #     ))
    #   }
    # }
    
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
    
    
    #print(cost_res)
    #print(round((1-(cost_res[row, col]) / sum(climb_up_cost, walk_right_cost)) * 100, 3))
    #print('##########################')
    # Return value: percentage of basepairs not crossed along alignments.
    return(round((
      1 - (cost_res[row, col]) / sum(climb_up_cost, walk_right_cost)
    ) * 100, 3))
  }



