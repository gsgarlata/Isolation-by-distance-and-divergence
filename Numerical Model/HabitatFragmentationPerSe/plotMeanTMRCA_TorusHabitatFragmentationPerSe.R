path_func = args[1]

source(paste0(path_func,'/weightedSij_torus_HLF_intime_parallel.R'))
source(paste0(path_func,'/Sij_torus_beforeHL&F.R'))
source(paste0(path_func,'/GetFragmCoordAllPossibleSize.R'))

#It defines the number of demes in the x-direction (which is equal to the y-direction) in a toroidal
#stepping stone model.
nd<-78

#"GetFragmCoordAllPossibleSize.R" is used to compute the relative spatial coordinates of the central deme within fragments of a specific size.
#For instance, for a habitat fragment of size 13 x 13, the central deme has spatial coordinates equal to x = 6 and y = 6.
res_fragm<-GetFragmCoordAllPossibleSize(nd_tot=nd,numb_fragm=3)

#From "res_fragm", we extract information regarding a HL&F configuration with 3x3 = 9 habitat fragments,
#each with 13 x 13 demes (fragm1). If one is interested in smaller fragment sizes,
#by setting fragm_id<-'fragm2', for instance, it would model fragment size of 11 x 11 demes.
fragm_id<-'fragm1'

#It extract the spatial coordinate of the central deme RELATIVE TO THE FRAGMENT SIZE, that is
#if the fragment size is 13 x 13, the central deme will have coordinates x = 6 and y = 6.
rel_x<-res_fragm[[fragm_id]]$rel_x
rel_y<-res_fragm[[fragm_id]]$rel_y

#It extracts the number of demes in the x- and y-direction of the habitat fragment.
dx<-res_fragm[[fragm_id]]$dx
dy<-res_fragm[[fragm_id]]$dy

#It extract the spatial coordinate of the central deme of the focal habitat fragment (X1) RELATIVE TO THE TORUS BEFORE HL&F. 
X1<-res_fragm$X1
Y1<-res_fragm$Y1
#It defines the deme effective population size
Ne<-50
#This must be a vector migration rate values in a two-dimensional (torus) stepping-stone model. 
migr_list<-c(0.02,0.04,0.06,0.1)

#It extracts the list of x-coordinates ("X2_list") and y-coordinates ("Y2_list"), each corresponding to a lineage sampled in a given fragment.
Y2_list<-res_fragm$Y2_list
X2_list<-res_fragm$X2_list

out_name<-paste0('Torus_',nd,'x',nd,'demes_to_',2*(dx-1),'x',2*(dy-1),'demes_K',Ne)

#The list of time points for which one wishes to compute the mean TMRCA of two alleles sampled in different habitat fragments.
time_list<-c(50,100,250,500,1000,2000,3000,4000)
#The number of CPUs to be used in the parallelization of the calculation.
cpu_n<-3

final_num<-NULL

for(migr in migr_list){
  
  #"weightedSij_torus_HLF_intime_parallel" is used to compute the second term in Eq. 4. 
  #That is, given a sampled lineage in a habitat fragment, it computes the weighted Sij 
  #over all possible locations of the ancestor (within the habitat fragment) at the time of HL&F.
  res<-weightedSij_torus_HLF_intime_parallel(migr,nd,dx,dy,rel_x,rel_y,X1,Y1,X2_list,Y2_list,time_list,cpu_n,path_func)
  
  #"Sij_torus_beforeHLF" is used to compute the Sij term (Eq. S2) in a continuous torus (before HL&F).
  res_torus<-Sij_torus_beforeHLF(X1,Y1,X2_list,Y2_list,migr,nd,path_func)
  #It combines the Sij calculations for before and after HL&F.
  temp_torus<-data.frame(geo=res_torus$geo.dist,time=0,Sij=res_torus$Sij)
  
  res_tot<-rbind(temp_torus,res)
  
  #It add the "time since HL&F" (the last two terms in Eq. 4).
  res_tot$Sij_time<-res_tot$time+res_tot$Sij
  
  num_tot<-subset(res_tot,geo>0)
  
  num_tot$scen<-paste0('m: ',migr)
  
  final_num<-rbind(final_num,num_tot)
}

#It add the mean TMRCA for two lineages in the same deme in the torus before HL&F (full calculation of Eq. 4). 
final_num$Tij_time = final_num$Sij_time + 2*Ne*nd*nd

plot_TijTime<-ggplot(data=final_num,aes(x=geo,y=Tij_time,group=factor(time),colour=factor(time)))+
  geom_point()+geom_line()+facet_wrap(vars(scen),scales = 'free',nrow = 1)+
  scale_color_viridis(name='Time since \n   HL&F',discrete=TRUE,direction=-1)+  
  xlab('Geographical Distance')+ylab('Mean Time to the MRCA')+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"),strip.text.x = element_text(size = 12))


