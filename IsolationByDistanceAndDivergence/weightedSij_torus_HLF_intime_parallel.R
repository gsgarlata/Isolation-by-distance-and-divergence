#This script parallelizes the"weightedSij_torus_HLF_intime" function, as to compute the weighted Sij over a list of geographical
#distances and different time points since HL&F.

#migr: migration in a torus (i.e. mx = migr/2 and my = migr/2).
#nd: number of demes over one dimension (x- or y-axis) of the torus before HL&F (i.e., ndx = nd and ndy = nd).
#dx: number of demes in the x-axis (it is the same for all habitat fragments).
#dy: number of demes in the y-axis (it is the same for all habitat fragments).
#rel_x: it defines the x-coordinate (relative to the coordinate in the torus before HL&F) of the center of the fragment.
#rel_y: it defines the x-coordinate (relative to the coordinate in the torus before HL&F) of the center of the fragment.
#X1: x-coordinate of one of the two lineages sampled in the habitat fragment.
#Y1: y-coordinate of one of the two lineages sampled in the habitat fragment.
#X2_list: list of x-coordinates, each corresponding to a lineage sampled in a given fragment.
#Y2_list: list of y-coordinates, each corresponding to a lineage sampled in a given fragment.
#time_list: vector of discrete times since HL&F.
#cpu_n: number of cpu for parallelization 
#path_function: path to needed functions


weightedSij_torus_HLF_intime_parallel<-function(migr,nd,dx,dy,rel_x,rel_y,X1,Y1,X2_list,Y2_list,time_list,cpu_n, path_function){
  
  require(parallel)
  require(foreach)
  require(doParallel)
  require(ggplot2)
  
  ####START functions####
  from_xy<-function(x,y,nside){
    return(x+(nside*y))
  } #it converts x and y coordinates in index of a 2D migration matrix
  
  to_xy<-function(k,nside){  
    return(c(k%%nside,k%/%nside))
  } #it converts an index of a 2D migration matrix into x and y coordinates
  
  make_2dmat<-function(mx,my,dx,dy){
    
    n_matrix = dx*dy #it defines the size of the 2D migration matrix considering all x and y pairwise combinations
                     #this is done by multiplying the number of demes in x and y direction.
    
    mat = matrix(0,nrow=n_matrix,ncol=n_matrix) #it creates a zero matrix of size n_matrix * n_matrix
    names_mat<-NULL
    for(i in 1:n_matrix){
      
      coord = to_xy(i-1,dx) #it convert matrix indexes in x and y coordinates. (It is sets to i-1 because matrix indexes in R start from '1')
      x = coord[1]
      y = coord[2]
      id = paste0(x,'_',y) #set the id based on the x-y coordinate of each index.
      pos1 = from_xy((x+1)%%dx,y,dy)+1 ##it converts x and y coordinates in index of a 2D migration matrix.
                                       #(We add '+1' because matrix indexes in R start from '1')
      pos2 = from_xy((x-1)%%dx,y,dy)+1
      pos3 = from_xy(x,(y+1)%%dy,dx)+1
      pos4 = from_xy(x,(y-1)%%dy,dx)+1
      
      mat[i,i] = 1-mx-my #it sets values on the diagonal equivalent to the probability of not moving both in the x- AND y- directions.
      mat[i,pos1] = mx/2 #it gives the probability of moving on the RIGHT (x-direction)
      mat[i,pos2] = mx/2 #it gives the probability of moving on the LEFT (x-direction)
      
      mat[i,pos3] = my/2 #it gives the probability of moving UP (y-direction)
      mat[i,pos4] = my/2 #it gives the probability of moving DOWN (y-direction)
      
      names_mat<-c(names_mat,id)
      
    }
    colnames(mat)<-names_mat
    rownames(mat)<-names_mat
    return(mat)
    
  } #it creates the 2D migration matrix for a torus-like space, given
                                         #'mx': dispersal rate in the 1d x-direction
                                         #'my': dispersal rate in the 1d y-direction
  ####END functions####
  
  dxToroid<-2*(dx-1) #it derives the number of demes in a 2D torus-like space equivalent to a 2D plane-like space with x-length equal to 'dx' 
  dyToroid<-2*(dy-1) #it derives the number of demes in a 2D torus-like space equivalent to a 2D plane-like space with y-length equal to 'dy'
  
  mx<-migr/2
  my<-migr/2
  
  mat_migr<-make_2dmat(mx,my,dxToroid,dyToroid) #it creates a 2D torus-like migration matrix
  
  worker.init <- function(path_function){
    require(Rcpp)
    sourceCpp(paste0(path_function,'/Cpp_functions_Toroidal_barrier.cpp'))
    sourceCpp(paste0(path_function,'/matrix_power.cpp'))
    source(paste0(path_function,'/weightedSij_torus_HLF_intime.R'))
  }
  
  cl <- makePSOCKcluster(cpu_n)
  clusterCall(cl, worker.init, path_function)
  registerDoParallel(cl)
  
  start.time <- Sys.time()
  #it parallelize the computation of weighted Sij over all the t* values tested (specified in time_list) (second term in Eq. 4)
  res<-foreach(i = 1:length(time_list), .combine = 'rbind',.packages=c('Rcpp','matrixcalc','reshape2'),.noexport = c('series_Sij','series_walk')) %dopar% {
    
    temp<-weightedSij_torus_HLF_intime(X1,Y1,X2_list,Y2_list,mat_migr,time_list[i],dx,dy,nd,rel_x,rel_y,migr)
    temp<-data.frame(geo=temp[,1],Sij=temp[,2])
    temp$time<-time_list[i]
    return(temp)
  }
  
  end.time <- Sys.time()
  parallel::stopCluster(cl)
  
  return(res)
  
}

