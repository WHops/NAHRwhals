#' New function
#' 
#' Has a record of setting both pos and neg values. 
#' @export
clean_sweep_matrix <- function(grid_list_f, debug = F){
  
  if (debug){
    grid_list_f = read.table("~/Desktop/grid_list_f_sample", sep='\t', header=T)
    grid_list_f = read.table("~/Desktop/grid_2.tsv", sep='\t', header=T)
  }
  
  grid_list_f = grid_list_f[order(grid_list_f$x),]
  #grid_list_f = data.frame(x=rep(1:3,5), y=c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5), z=c(-1,0,1,0,1,0,1,0,-1,0,1,0,1,0,0))
  grid_list_f_2 = grid_list_f

  grid_list_f = grid_list_f[grid_list_f$z != 0,]
  row.names(grid_list_f) = 1:dim(grid_list_f)[1]
  
  for (npoint in 1:dim(grid_list_f)[1]){
    y_vals = sort(grid_list_f[grid_list_f$x == grid_list_f[npoint,'x'],'y'])
    x_vals = sort(grid_list_f[grid_list_f$y == grid_list_f[npoint,'y'],'x'])
    
    y_vals = y_vals[!y_vals == grid_list_f[npoint,'y']]
    x_vals = x_vals[!x_vals == grid_list_f[npoint,'x']]
    # a = grid_list_f[grid_list_f$x == grid_list_f[npoint,'x'],]
    # b = grid_list_f[grid_list_f$y == grid_list_f[npoint,'y'],]
    
    if ((length(x_vals) > 0) & (length(y_vals) > 0)){
      for (i in 1:length(x_vals)){
        for (j in 1:length(y_vals)){
          #print(npoint)
          #print(paste0('i: ', i))
          #print(paste0('j: ', j))

          #print(grid_list_f[npoint,'x'])
          #print(grid_list_f[npoint,'y'])
          
          
          p1 = grid_list_f[npoint,c('x','y')]
          p2 = data.frame(x=x_vals[i], y=p1$y)
          p3 = data.frame(x=p1$x, y=y_vals[j])
          p4 = data.frame(x=x_vals[i], y=y_vals[j])
          
          all_p = rbind(p1, rbind(p2, rbind(p3, p4)))
          
          
          in_list = dplyr::left_join(all_p, grid_list_f[c('x','y', 'z')], by=c('x','y'))
          in_list_z = in_list$z
          #in_list[is.na(in_list$z),]$z = abs(prod(in_list[!is.na(in_list$z),]$z)) * max((prod(in_list[!is.na(in_list$z),]$z)))
          #missing = sum(is.na(in_list_z)) #s
          missing = dplyr::setdiff(all_p,grid_list_f[c('x','y')])
          if ( dim(missing)[1] > 0){
            #print('Missing point detected! Adding...')
            missing$z = sign(prod(in_list[!is.na(in_list$z),]$z)) * max((in_list[!is.na(in_list$z),]$z))
            #print(missing)
            grid_list_f_2 = rbind(grid_list_f_2, missing)
          }
        }
      }
    }
    }

  grid_list_f_2 = grid_list_f_2[order(grid_list_f_2$x),]
  # p = ggplot2::ggplot(dplyr::distinct(grid_list_f_2)) + ggplot2::geom_tile(ggplot2::aes(
  #   x = x,
  #   y = y,
  #   fill = sign(z) * log10(abs(z))
  # )) +
  #   ggplot2::scale_fill_gradient2(
  #     low = 'red',
  #     mid = 'black',
  #     high = 'blue',
  #   ) 
  # p
  return(dplyr::distinct(grid_list_f_2))
}

# grid_list_f = read.table("~/Desktop/grid_2.tsv", sep='\t', header=T)
# g2 = clean_sweep_matrix(grid_list_f, debug=F)
# g3 = clean_sweep_matrix(g2, debug=F)
# g4 = clean_sweep_matrix(g3, debug=F)
# 
# p = ggplot2::ggplot(dplyr::distinct(g4)) + ggplot2::geom_tile(ggplot2::aes(
#   x = x,
#   y = y,
#   fill = sign(z) * log10(abs(z))
# )) +
#   ggplot2::scale_fill_gradient2(
#     low = 'red',
#     mid = 'black',
#     high = 'blue',
#   )
# p
# 
# 
# g2[duplicated(g2[, c('x', 'y')]),]
# 
# ggplot(g2) + geom_point(aes(x=x, y=y, col=sign(z)))
# g3 = clean_sweep_matrix(g2, debug=F)
# ggplot(g3) + geom_point(aes(x=x, y=y, col=sign(z)))
# 


