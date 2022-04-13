#' New function
#' 
#' Has a record of setting both pos and neg values. 
#'#' @export
clean_sweep_matrix <- function(grid_list, debug = F){
  
  if (debug){
    grid_list = read.table("~/Desktop/grid_list_sample", sep='\t', header=T)
  }
  
  grid_list = grid_list[order(grid_list$x),]
  #grid_list = data.frame(x=rep(1:3,5), y=c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5), z=c(-1,0,1,0,1,0,1,0,-1,0,1,0,1,0,0))
  grid_list_2 = grid_list

  grid_list = grid_list[grid_list$z != 0,]
  row.names(grid_list) = 1:dim(grid_list)[1]
  
  for (npoint in 1:dim(grid_list)[1]){
    y_vals = sort(grid_list[grid_list$x == grid_list[npoint,'x'],'y'])
    x_vals = sort(grid_list[grid_list$y == grid_list[npoint,'y'],'x'])
    
    y_vals = y_vals[!y_vals == grid_list[npoint,'y']]
    x_vals = x_vals[!x_vals == grid_list[npoint,'x']]
    # a = grid_list[grid_list$x == grid_list[npoint,'x'],]
    # b = grid_list[grid_list$y == grid_list[npoint,'y'],]
    
    if ((length(x_vals) > 0) & (length(y_vals) > 0)){
      for (i in 1:length(x_vals)){
        for (j in 1:length(y_vals)){
          #print(npoint)
          #print(paste0('i: ', i))
          #print(paste0('j: ', j))

          #print(grid_list[npoint,'x'])
          #print(grid_list[npoint,'y'])
          
          
          p1 = grid_list[npoint,c('x','y')]
          p2 = data.frame(x=x_vals[i], y=p1$y)
          p3 = data.frame(x=p1$x, y=y_vals[j])
          p4 = data.frame(x=x_vals[i], y=y_vals[j])
          all_p = rbind(p1, rbind(p2, rbind(p3, p4)))
          
          in_list = left_join(all_p, grid_list[c('x','y', 'z')], by=c('x','y'))
          in_list_z = in_list$z
          #in_list[is.na(in_list$z),]$z = abs(prod(in_list[!is.na(in_list$z),]$z)) * max((prod(in_list[!is.na(in_list$z),]$z)))
          
          missing = setdiff(all_p,grid_list[c('x','y')])
          if ((dim(missing)[1]) > 0){
            #print('Missing point detected! Adding...')
            missing$z = sign(prod(in_list[!is.na(in_list$z),]$z)) * max((in_list[!is.na(in_list$z),]$z))
            #print(missing)
            grid_list_2 = rbind(grid_list_2, missing)
          }
        }
      }
    }
    }


  return(dplyr::distinct(grid_list_2))
}

grid_list = read.table("~/Desktop/grid_list_sample", sep='\t', header=T)
g2 = clean_sweep_matrix(grid_list, debug=F)
g2[duplicated(g2[, c('x', 'y')]),]

# ggplot(g2) + geom_point(aes(x=x, y=y, col=sign(z)))
# g3 = clean_sweep_matrix(g2, debug=F)
# ggplot(g3) + geom_point(aes(x=x, y=y, col=sign(z)))



