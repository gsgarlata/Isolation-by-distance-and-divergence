#This script is used to compute the Sij term (Eq. S2) in a continuous torus (before HL&F) 
#over a list of geographical locations (X2_list, Y2_list) in comparison to point X1, Y1.


#X1: x-coordinate of one of the two lineages sampled in the continuous torus.
#Y1: y-coordinate of one of the two lineages sampled in the continuous torus.
#X2_list: list of x-coordinates, each corresponding to a lineage sampled at a different geographical location than X1, Y1.
#Y2_list: list of y-coordinates, each corresponding to a lineage sampled at a different geographical location than X1, Y1.
#migr: migration in a torus (i.e. mx = migr/2 and my = migr/2).
#nd: number of demes over one dimension (x- or y-axis) of the torus before HL&F (i.e., ndx = nd and ndy = nd).
#path_function: path to needed functions


Sij_torus_beforeHLF<-function(X1,Y1,X2_list,Y2_list,migr,nd,path_function){
  require(Rcpp)
  sourceCpp(paste0(path_function,'/Cpp_functions_Toroidal_barrier.cpp'))
  
  
  final<-NULL
  
  for(a in 1:length(Y2_list)){
    Y2<-Y2_list[a]
  for(s in 1:length(X2_list)){
   X2<-X2_list[s]
    
    temp<-sumSij(abs(Y2-Y1),abs(X2-X1), nd, nd,migr)
    geo_dists<-abs(dist(rbind(c(X1,Y1),c(X2,Y2)),method='euclidean'))
    final<-rbind(final,c(geo_dists,temp))
  }
  }

  data_test<-data.frame(geo.dist=final[,1],Sij=final[,2])
  
  return(data_test)
}
