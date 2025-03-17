path_func = args[1]

source(paste0(path_func,'/weightedSij_torus_HLF_intime_parallel.R'))
source(paste0(path_func,'/Sij_torus_beforeHL&F.R'))
source(paste0(path_func,'/GetFragmCoordAllPossibleSize.R'))
source(paste0(path_func,'/AverageCoalTime_PopSizeChange_Torus.R'))

#It defines the number of demes in the x-direction (which is equal to the y-direction) in a toroidal
#stepping stone model.
nd<-78
numb_fragm<-3

#"GetFragmCoordAllPossibleSize.R" is used to compute the relative spatial coordinates of the central deme within fragments of a specific size.
#For instance, for a habitat fragment of size 13 x 13, the central deme has spatial coordinates equal to x = 6 and y = 6.
res_fragm<-GetFragmCoordAllPossibleSize(nd_tot=nd,numb_fragm)

#It extract the spatial coordinate of the central deme of the focal habitat fragment (X1) RELATIVE TO THE TORUS BEFORE HL&F.
X1<-res_fragm$X1
Y1<-res_fragm$Y1
#It defines the deme effective population size.
Ne<-50
#This must be a vector migration rate values in a two-dimensional (torus) stepping-stone model. 
migr_list<-c(0.02,0.04,0.06,0.1)

#It extracts the list of x-coordinates ("X2_list") and y-coordinates ("Y2_list"), each corresponding to a lineage sampled in a given fragment.
X2_list<-res_fragm$X2_list
Y2_list<-res_fragm$Y2_list

#The list of time points for which one wishes to compute the mean TMRCA of two alleles sampled in different habitat fragments.
time_list<-c(50,100,250,500,1000,2000,3000,4000)
#The number of CPUs to be used in the parallelization of the calculation.
cpu_n<-3

fragm_list<-c(13,3)

fragm_rows<-do.call(rbind,res_fragm)

fragm_infos<-fragm_rows[fragm_rows$dx %in% fragm_list,]

tot_num<-NULL

for(i in 1:nrow(fragm_infos)){

  dx<-fragm_infos[i,'dx']
  dy<-fragm_infos[i,'dy']
  rel_x<-fragm_infos[i,'rel_x']
  rel_y<-fragm_infos[i,'rel_y']
  
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
  
  #It computes the non-equilibrium mean TMRCA in a torus undergoing habitat contraction.
  df_temp<-AverageCoalTime_PopSizeChangeT00Torus(nd_full=nd,dx,dy,migr/2,Ne=Ne,time_end=time_list[length(time_list)],path_func)
  df_temp<-rbind(c(0,2*Ne*nd*nd),df_temp)
  
  selected_vals<-df_temp[df_temp$time %in% c(0,time_list),]
  
  num_tot$t00_noneq<-NA
  
  for(k in 1:nrow(selected_vals)){
    
    coord<-which(num_tot$time==selected_vals$time[k])
    
    num_tot$t00_noneq[coord]<-selected_vals$t00[k]
    
  }
  
  
  num_tot$scen<-paste0('m: ',migr)
  
  final_num<-rbind(final_num,num_tot)
}

#It add the mean TMRCA for two lineages in the same deme in the torus before HL&F (full calculation of Eq. 4).
final_num$t00<-2*Ne*nd*nd

final_num$Fst<-1 - (final_num$t00_noneq/(final_num$Sij_time + final_num$t00))
final_num$ni<-final_num$Fst/(1-final_num$Fst)


final_num$fragm_size<-paste0(2*(fragm_list[i]-1),'x',2*(fragm_list[i]-1))

tot_num<-rbind(tot_num,final_num)

}


max_frag_size<-2*(fragm_infos$dx[1]-1)
min_frag_size<-2*(fragm_infos$dx[2]-1)


large_frag_df<-subset(tot_num,fragm_size==paste0(max_frag_size,'x',max_frag_size))

g1 = ggplot(data=large_frag_df,aes(x=time,y=Fst,group=factor(round(geo)),colour=factor(round(geo))))+
  geom_point()+
  geom_line()+
  scale_color_viridis(name='Geographical\n   Distance',discrete=TRUE,direction=-1)+
  facet_wrap(~scen, scales = 'free',ncol=4)+
  xlab('Time since HL&F')+ylab(bquote(F[st]))+
  theme_bw()+
  ggtitle(subtitle= '',label=paste0("\nHabitat fragment size: ",max_frag_size,"x",max_frag_size))+
  theme(strip.background =element_rect(fill="white"),strip.text.x = element_text(size = 12),
        plot.title = element_text(hjust = 0.5))



small_frag_df<-subset(tot_num,fragm_size==paste0(min_frag_size,'x',min_frag_size))

g2 = ggplot(data=small_frag_df,aes(x=time,y=Fst,group=factor(round(geo)),colour=factor(round(geo))))+
  geom_point()+
  geom_line()+
  scale_color_viridis(name='Geographical\n   Distance',discrete=TRUE,direction=-1)+
  facet_wrap(~scen, scales = 'free',ncol=4)+
  xlab('Time since HL&F')+ylab(bquote(F[st]))+
  theme_bw()+
  ggtitle(subtitle = '',label=paste0("\nHabitat fragment size: ",min_frag_size,"x",min_frag_size))+
  theme(strip.background =element_rect(fill="white"),
        strip.text.x = element_text(size = 12),plot.title = element_text(hjust = 0.5))


final_plot<-grid.arrange(g1, g2, nrow = 2)
