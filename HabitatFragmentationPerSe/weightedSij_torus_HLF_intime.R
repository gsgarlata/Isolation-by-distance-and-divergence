#This script computes the second term in Eq. 4. That is, given a sampled lineage in a habitat fragment, it computes the weighted Sij 
#over all possible locations of the ancestor (within the habitat fragment) at the time of HL&F.
#Note that Fst is calculated between the lineage sampled at X1,Y1 and all other lineages sampled in different fragments (x-coord list:  X2_list; y-coord list:  Y2_list).

#X1: x-coordinate of one of the two lineages sampled in the habitat fragment.
#Y1: y-coordinate of one of the two lineages sampled in the habitat fragment.
#X2_list: list of x-coordinates, each corresponding to a lineage sampled in a given fragment.
#Y2_list: list of y-coordinates, each corresponding to a lineage sampled in a given fragment.
#mat_migr: a 2D torus-like migration matrix.
#time_barrier: time of HL&F.
#dx: number of demes in the x-axis (it is the same for all habitat fragments).
#dy: number of demes in the y-axis (it is the same for all habitat fragments).
#nd: number of demes over one dimension (x- or y-axis) of the torus before HL&F (i.e., ndx = nd and ndy = nd).
#rel_x: it defines the x-coordinate (relative to the coordinate in the torus before HL&F) of the center of the fragment.
#rel_y: it defines the x-coordinate (relative to the coordinate in the torus before HL&F) of the center of the fragment.
#migr: migration in a torus (i.e. mx = migr/2 and my = migr/2).

weightedSij_torus_HLF_intime<-function(X1,Y1,X2_list,Y2_list,mat_migr,time_barrier,dx,dy,nd,rel_x,rel_y,migr){
  
  require(reshape2)
  require(matrixcalc)
  
  ####FUNCTIONS####
  
  combo_mat<-function(x,df1,df2){
    df1[,3]<-df2[x,1]
    df1[,4]<-df2[x,2]
    return(df1)
  }
  
  combo_mat_evol<-function(x,df1,df2,df3){
    df1[,5]<-df2[x,1]
    df1[,6]<-df2[x,2]
    df1[,7]<-df3[x,1]
    df1[,8]<-df3[x,2]
    
    return(df1)
  }
  
  
  manhattan_dist <- function(a, b){
    dist <- abs(a-b)
    dist <- sum(dist)
    return(dist)
  }
  
  euclidean_dist<- function(a, b) sqrt(sum((a - b)^2))
  ####END FUNCTIONS###
  
  dxToroid<-2*(dx-1) #it defines the number of demes in a 1d circular 
                     #stepping stone equivalent to a linear 1d stepping stone (half size). 
                     #Maruyama, 1971 showed that L*L plane is equivalent to a 2L * 2L toroid.
  dyToroid<-2*(dy-1)
  

  time_matrix<-PowerMatrix(mat_migr,time_barrier) #It compute the time-convolution of migration probabilities 
  colnames(time_matrix)<-colnames(mat_migr)                                     #for each of the two axes
  rownames(time_matrix)<-rownames(mat_migr)                                     #for each of the two axes
  
  mat_deme<-time_matrix[which(rownames(time_matrix)==paste0(rel_x,'_',rel_y)),] #Subset only the row corresponding to the sampled
                                                                                #deme, which in this case is at the center of the 2D fragment
                                                                                #with coordinate x=rel_x and y=rel_y.
  #This series of lines help to give geographical coordinates
  #to the torus-like space, considering that it includes back
  #and forward movement 'reflective edges'. 
  toroidal<-matrix(NA,nrow=dxToroid,ncol=dyToroid) #matrix with row X col corresponding to all possible X and Y pairwise combinations
  colnames(toroidal)<-0:(ncol(toroidal)-1) #give column names corresponding to the numbering in a 1d-circular model [0:(d-1)]
  rownames(toroidal)<-0:(nrow(toroidal)-1) #give row names corresponding to the numbering in a 1d-circular model [0:(d-1)]

  mat1x<-toroidal #make a copy of the toroidal matrix for further editing
  first_row_1x<-((dxToroid/2) + 2):nrow(mat1x) #get row indexes for the 'reflective' ('moving backward') demes in a 1d circular model.
  sec_row_1x<-((dxToroid/2)-1):0 #get row indexes for the 'moving forward' demes in a 1d circular model.
  
  #the same as for the x-direction but for the y-direction.
  first_row_1y<-((dyToroid/2) + 2):ncol(mat1x)
  sec_row_1y<-((dyToroid/2)-1):0
  
  #re-number the row and columns of the matrix considering that a 2D torus includes back
  #and forward movement due to 'reflective edges'.  
  rownames(mat1x)[first_row_1x]<-sec_row_1x[1:length(first_row_1x)]
  colnames(mat1x)[first_row_1y]<-sec_row_1y[1:length(first_row_1y)]
  
  #change format of the matrix to dataframe
  coord_tor<-melt(mat1x)
  coord_tor$value<-NULL
  names(coord_tor)<-c('x1','y1')
  
  #change format of the matrix to dataframe
  tor1x<-melt(toroidal)
  tor1x$value<-NULL
  
  #keep information on the 'absolute' inidexing of the matrix [0:(d-1)]
  coord_tor$idx1<-tor1x$Var1
  coord_tor$idy1<-tor1x$Var2
  
  #it replicate for habitat fragment A2, the steps done for the habitat fragment A1 
  coord_tor$x2<-NA
  coord_tor$y2<-NA
  coord_tor$idx2<-NA
  coord_tor$idy2<-NA
  
  coord_tor2<-coord_tor[,c(1,2)]
  tor2x<-coord_tor[,c(3,4)]
  
  mat_cb<-do.call('rbind',lapply(1:nrow(coord_tor2),combo_mat_evol,df1=coord_tor,df2=coord_tor2,df3=tor2x))
  
  #keep information on the absolute indexes of the matrix for identifying the migration probabilities from 'mat_migr'
  mat_cb$xy_A1<-paste0(mat_cb$idx1,'_',mat_cb$idy1)
  mat_cb$xy_A2<-paste0(mat_cb$idx2,'_',mat_cb$idy2)
    
  #iterate over all Y and X coordinates
  final<-NULL
  for(k in 1:length(Y2_list)){
    
    Y2<-Y2_list[k]
    
    for(i in 1:length(X2_list)){
      X2<-X2_list[i]
      
      mat_cb$i<-abs((mat_cb$x1 - rel_x + X1) - (mat_cb$x2 - rel_x + X2)) #it compute the absolute distance between all possible unknown ancestors
                                                                         #for both fragment A1 and A2. (x-direction)
      mat_cb$j<-abs((mat_cb$y1 - rel_y + Y1) - (mat_cb$y2 - rel_y + Y2)) #it compute the absolute distance between all possible unknown ancestors
                                                                         #for both fragment A1 and A2. (y-direction)
      
      mat_cb$i_j<-paste0(mat_cb$i,'_',mat_cb$j)
      ij_unique<-unique(mat_cb$i_j)
      temp_cb<-mat_cb[match(ij_unique, mat_cb$i_j),]
      temp_cb$Sij<-series_Sij(temp_cb$i, temp_cb$j, nd, nd, migr)#it computes the average amount of time for A1 and A2 alleles at distance
                                                                 #'i' and 'j', in the x and y direction respectively,  to be in the same deme.
      mat_cb$Sij<-temp_cb$Sij[match(mat_cb$i_j,temp_cb$i_j)]
      
      mat_cb$pA1<-mat_deme[match(mat_cb$xy_A1,names(mat_deme))] #extract the transition probabilities from the toroidal matrix for fragment A1
      mat_cb$pA2<-mat_deme[match(mat_cb$xy_A2,names(mat_deme))] #extract the transition probabilities from the toroidal matrix for fragment A2
      
      mat_cb$p_tot<-(mat_cb$pA1*mat_cb$pA2) #It multiplies two two probabilities of A1 and A2
      
      mat_cb$Et<-mat_cb$p_tot*mat_cb$Sij #It weights Sij by the probability of the unknown ancestors for A1 and A2 at time t*
      Et<-sum(mat_cb$Et) #It sums all weighted Sij values across all possible unknown ancestors at time t*
      geo<-euclidean_dist(c(X1,Y1),c(X2,Y2)) #it computes euclidean distance between the two sampled A1 and A2 demes
      res<-c(geo,Et)
      final<-rbind(final,res)
      
    }
  }
  
  return(final)
  
}
