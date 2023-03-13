


#' carry_out_compressed_sv
#'
#' @description Carry out an SV IN GRIDSPACE
#'
#' @param bitl  an nxm matrix (n=gridpoints_x, m=gridpoints_y)
#' @param input_ins a vector with instruction: p1, p2, sv.
#'
#' @author Wolfram Höps
#' @export
carry_out_compressed_sv <- function(bitl, input_ins) {
  
  pair = as.numeric(input_ins[1:2])
  action = input_ins[3]
  
  # Modify the matrix accordingly
  if (action == 'del') {
    
    bitl_mut = bitl[, -c(as.numeric(pair[1]):(as.numeric(pair[2]) - 1))]
    
  } else if (action == 'dup') {
    
    bitl_mut = cbind(bitl[, 1:as.numeric(pair[2])],
                     bitl[, as.numeric(pair[1] + 1):dim(bitl)[2]])
    
    if (pair[2] - pair[1] == 1){
      if (pair[1] == 1){
        # W, 10th March 2022, Changed this a bit to account for cases
        # with only 2 columns. 
        colnames(bitl_mut)[dim(bitl_mut)[2]] = colnames(bitl_mut)[2] 
      } else if (pair[2] == dim(bitl)[2]){
        colnames(bitl_mut)[dim(bitl_mut)[2]] = colnames(bitl_mut)[dim(bitl_mut)[2]-1]
      }
    }
    
  } else if (action == 'inv') {
    
    # W, 8th March 2022. 
    # Noticed a problem with one-column inversions. Here, the expression 
    # -bitl[, (as.numeric(pair[2]) - 1):(as.numeric(pair[1]) + 1)])
    # Returns a one-line dataframe, which loses its column name. 
    # Therefore, we now handle one-column inversions differently now. 
    if ((pair[2] - pair[1]) == 1){
      bitl_mut = bitl
    } else if (pair[2] - pair[1] > 1){
      
        # A lot of bordercase handling...
      
        # Border case 1: The first block is part of inversion
        if ((pair[1] == 1) & (pair[2] != dim(bitl)[2])){
          bitl_mut = cbind( -bitl[, (as.numeric(pair[2])):(as.numeric(pair[1]))],
                             bitl[, as.numeric(pair[2]+1):dim(bitl)[2]])
          # Manually edit the last entry of the header. This is necessary if pair[2] == dim-1.
          # This is in response to an error ('error #2, may11th 2022')
          colnames(bitl_mut)[dim(bitl_mut)[2]] = colnames(bitl)[dim(bitl_mut)[2]]
          
          
        # Border case 2: The last block is part of inversion
        } else if ((pair[1] > 1) & (pair[2] == dim(bitl)[2])){
          bitl_mut = cbind( bitl[, 1:(as.numeric(pair[1])-1), drop=F],
                            -bitl[, (as.numeric(pair[2])):(as.numeric(pair[1])), drop=F])
          # Manually edit the last entry of the header. This is necessary if pair[1] == 2
          # This is in response to an error ('error #2, may11th 2022')
          colnames(bitl_mut)[1] = colnames(bitl)[1]
          
        # Border case 3: The inversion goes from first to last block
        } else if ((pair[1] == 1) & (pair[2] == dim(bitl)[2])){
          bitl_mut = -bitl[, (as.numeric(pair[2])):(as.numeric(pair[1]))]
          
        # The standard case
        } else {
        bitl_mut = cbind(cbind(bitl[, 1:as.numeric(pair[1])],
                               -bitl[, (as.numeric(pair[2]) - 1):(as.numeric(pair[1]) + 1)]),
                         bitl[, as.numeric(pair[2]):dim(bitl)[2]])
        } 
      
        # Border case handling. If inversion is exactly one column, 
        # (i.e. pair 3-5), then we assign the colnames manually. 
        if (pair[2] - pair[1] == 2){
          colnames(bitl_mut) = colnames(bitl)
        } else if (pair[2] == dim(bitl)[2]) {
          colnames(bitl_mut)[pair[2]] = colnames(bitl_mut)[pair[1]]
        }
    }
    } else if (action == 'ref'){
      bitl_mut = bitl
    }
  
  # It can happen that only one column is left. In that case...
  
  if (any((class(bitl_mut) != c('matrix', 'array')))){
    bitl_mut = as.matrix(bitl_mut)
    row.names(bitl_mut) = row.names(bitl)
    colnames(bitl_mut) = colnames(bitl)[1]
  }
  
  # if ("" %in% colnames(bitl_mut)){
  #   browser()
  # }
  # W, 7th March 2022. Excluding this line since bitl now has 
  # meaningful rownames and colnames. 
  #colnames(bitl_mut) = 1:dim(bitl_mut)[2]

  # Turn minus-zero into zero. 
  bitl_mut[bitl_mut == 0] = 0
  return(bitl_mut)
  
}

#' modify_gridmatrix
#' 
#' Introduce a chain of mutations to gridmatrix. 
#' Is used after mutation simulation to verify/visualize
#' the result. 
#' @param gmx  GridMatriX. Row/Colnames don't matter. 
#' @param r1 mutation instruction. Not sure about format right now. TODO. 
#' @export
modify_gridmatrix <- function(gmx, r1){
  
  # If there is no mutation to be done, return matrix as it is. 
  if (r1[1,'mut1'] == 'ref'){
    return(gmx)
  }
  
  # How many mutations do we have to run?
  
  #nmut = length(r1[, colSums(is.na(r1)) == 0]) - 1 # Outdated since res has gotten more columns
  nmut = sum(c('mut1','mut2','mut3') %in% colnames(r1[,colSums(is.na(r1))<nrow(r1)]))
  # Run each mutation. 
  for (i in 1:nmut){
    instr = r1[,paste0('mut', i)]
    start = as.numeric(strsplit(instr, '_')[[1]][1])
    end = as.numeric(strsplit(instr, '_')[[1]][2])
    action = strsplit(instr, '_')[[1]][3]
    
    # Run ONE mutation
    gmx = carry_out_compressed_sv(gmx, c(start, end, action))
  }
  
  # Not ENTIRELY sure, but for some reason we have to transpose here to get the 
  # correct result back apparently? Might be wrong. Keep an eye on it. 
  
  # Update W, 26th Jul 2022: 
  # No why the hell would you transpose here? xD
  # gmx = t(gmx)
  
  # Simple row and colnames for easy plotting. 
  colnames(gmx) = 1:ncol(gmx)
  row.names(gmx) = 1:nrow(gmx)
  
  return(gmx)
}


