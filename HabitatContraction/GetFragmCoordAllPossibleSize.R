

GetFragmCoordAllPossibleSize<-function(nd_tot,numb_fragm){
  
  #"available_demes" calculation is done since teh genetic differentiation expected in a plane of length "nd_tot"
  #can be reproduced by considering the half of the demes in a torus. This is because the torus is a "wrapped"
  #two-dimensional plane.
  available_demes<-nd_tot/2
  
  #Given the size of half torus and the number of modeled habitat fragments,
  #we compute the maximum number of demes that can be included in a fragment.
  #For instance, in the case of available_demes = 39 and numb_fragm = 3 (in x-direction),
  #we can have 3 (x-direction) x 3 (y-direction) (=9) fragments of maximum 13 x 13 demes (39/3 = 13).
  max_demes_fragm<-(available_demes/numb_fragm)
  
  #It considers some adjustments depending on the full size of the torus and the modeled number of fragments.
  #It computes the (0-based) spatial coordinates of the demes in the habitat fragment, relative to the fragment size.
  if((max_demes_fragm %% 2) == 0){
    max_demes_fragm<-((available_demes-3)/numb_fragm)
    coord_max<-max_demes_fragm-1
    series_demes<-seq(0,coord_max)
    halfpos<-median(series_demes)
    
  }else{
    coord_max<-max_demes_fragm-1
    series_demes<-seq(0,coord_max)
    halfpos<-median(series_demes)
  }
  
  #It extracts the spatial coordinate (equal for the x- and y-axis) of the central deme
  #within the habitat fragment used as focal habitat fragment (X1) RELATIVE to the FULL TORUS BEFORE HL&F.
  X1_list<-median(series_demes)
  
  #It loops over each of the remaining habitat fragments (i.e., excluding the focal habitat fragment X1)
  for(u in 2:numb_fragm){
    
    #"minpos" is the spatial coordinate of the deme that is located on the LEFT edge of the habitat fragment.
    #"maxpos" is the spatial coordinate of the deme that is located on the RIGHT edge of the habitat fragment.
    minpos<-series_demes[length(series_demes)]+1
    maxpos<-minpos+coord_max
    #It lists the spatial coordinate of all demes within the fragment from the LEFT to the RIGHT edge.
    series_demes_fragm<-seq(minpos,maxpos)
    
    #It derives the spatial coordinate of the central deme RELATIVE to the HABITAT FRAGMENT.
    if((max_demes_fragm %% 2) == 0){
      halfpos_fragm<-(length(series_demes_fragm))/2
    }else{
      halfpos_fragm<-(length(series_demes_fragm)-1)/2
    }
    X1_fragm<-median(series_demes_fragm)
    
    X1_list<-c(X1_list,X1_fragm)
    series_demes<-series_demes_fragm
    
  }
  
  res_list<-list(X1=X1_list[1],Y1=X1_list[1],X2_list=X1_list,Y2_list=X1_list)
  
  fragm1<-data.frame(rel_x=halfpos,rel_y=halfpos,dx=max_demes_fragm,dy=max_demes_fragm)
  
  res_list$fragm1<-fragm1
  
  count = 2
  
  #It loops over all the possible sizes of habitat fragments, given the full size of the torus before HL&F
  #and the modeled number of fragments. It considers decreasing habitat fragment sizes, by removing demes on both
  #LEFT and RIGHT edges of the habitat fragment.
  for(i in 1:halfpos){
    if(halfpos_fragm>1){
      
      
      if((max_demes_fragm %% 2) == 0){
        max_demes_fragm<-max_demes_fragm-1
        coord_max<-max_demes_fragm
      }else{
        
        max_demes_fragm<-max_demes_fragm-2
        coord_max<-max_demes_fragm-1
      }
      
      series_demes<-seq(0,coord_max)
      
      halfpos_fragm<-median(series_demes)
      
      temp_fragm<-data.frame(rel_x=halfpos_fragm,rel_y=halfpos_fragm,dx=max_demes_fragm,dy=max_demes_fragm)
      
      res_list[[paste0('fragm',count)]]<-temp_fragm
      
      count = count + 1
    }
  }
  
  return(res_list) 
}
