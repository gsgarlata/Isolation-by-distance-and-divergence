require(ggplot2)
require(viridis)
require(magrittr)
require(patchwork)

#Defines the directory path where the required R functions are saved.
path_func<-args[1]

#"AverageCoalTime_PopSizeChange_Torus.R" is used to compute the mean TMRCA of two alleles sampled in the same deme,
#in a torus undergoing Habitat Contraction.
source(paste0(path_func,'/AverageCoalTime_PopSizeChange_Torus.R'))

#"GetFragmCoordAllPossibleSize.R" is used to compute the relative spatial coordinates of the central deme within fragments of a specific size.
#For instance, for a habitat fragment of size 13 x 13, the central deme has spatial coordinates equal to x = 6 and y = 6.
source(paste0(path_func,'/GetFragmCoordAllPossibleSize.R'))

###function PLOTTING###
make_plot <- function(data,xlim_sub) {
  
  min_ylim<-min(data$t00)
  max_ylim<-max(data$t00)
  
  min_ylim_sub<-min(data[which(data$time==xlim_sub),'t00'])
  
  main <- ggplot(data=data,aes(x=time,y=t00,group=scen,colour=scen))+
    geom_line(size=2)+
    facet_wrap(~migr, scales = 'fixed',ncol=4)+ylim(min_ylim,max_ylim)+
    scale_color_viridis(name='Fragment\n    Size',discrete=TRUE)+  
    xlab('Time since HL&F')+
    ylab('Mean Time to the MRCA')+
    theme_bw()+
    theme(strip.background =element_rect(fill="white"),strip.text.x = element_text(size = 12),
          plot.title = element_text(hjust = 0.5))
  
  
  sub <- ggplot(data=data,aes(x=time,y=t00,group=scen,colour=scen))+
    geom_line(size=1,show.legend =FALSE)+
    scale_color_viridis(name='',discrete=TRUE)+  
    xlab('')+
    ylab('')+ylim(min_ylim_sub,max_ylim)+
    theme_bw()+scale_x_continuous(limits=c(0,xlim_sub),breaks=c(0,round(xlim_sub/2),xlim_sub))+
    theme(legend.position = "none",
          plot.margin = margin(0, 0 , 0, 0),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size=8),
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          legend.background = element_rect(fill = "transparent"), # get rid of legend bg
          legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    )
  
  
  main + inset_element(sub, -0.05, 0, .35, .35)
}
###function PLOTTING###

#The deme effective population size. 
Ne<-50
#This must be a vector migration rate values in a one directional circular stepping-stone model. 
migr_list<-c(0.01,0.02,0.03,0.05) 
#It defines the number of demes on the x-axis (equal to the number of demes on the y-axis) before Habitat Contraction.
nd_full<-78
#It defines the time duration of the contraction for which one wishes to compute the mean TMRCA.
time_end<-4000
#It defines the number of habitat fragments to be modeled over the x-axis and the y-axis.
#It thus consider that the initial torus undergoes HL&F generating 3 x 3 = 9 habitat fragments.
numb_fragm<-3 

#Ignores the Warnings.
res_fragm<-do.call(rbind,GetFragmCoordAllPossibleSize(nd_tot=nd_full,numb_fragm))

coord_fragm<-grep('fragm',rownames(res_fragm))

dx_fragm<-res_fragm$dx[coord_fragm]

migr_df<-NULL

for(migr in migr_list){
  
  final_df<-NULL
  
  for(u in 1:length(dx_fragm)){
    
    dx<-dx_fragm[u]
    dy<-dx_fragm[u] #dimensions of a plane fragment and not a torus/circular stepping-stone
    
    df_temp<-AverageCoalTime_PopSizeChangeT00Torus(nd_full,dx,dy,migr,Ne,time_end, path_func)
    df_temp<-rbind(data.frame(time=0,t00=2*Ne*nd_full*nd_full),df_temp)
    df_temp$scen<-paste0((2*(dx-1)),'x',(2*(dx-1)))
    
    final_df<-rbind(final_df,df_temp)
  }
  
  final_df$migr<-paste0('m: ',2*migr)
  
  migr_df<-rbind(migr_df,final_df)
  
}

migr_df$scen<-factor(migr_df$scen,levels=unique(migr_df$scen))
migr_df$migr<-factor(migr_df$migr,levels=unique(migr_df$migr))

p <- migr_df %>% 
  split(.$migr) %>% 
  lapply(make_plot,xlim_sub=50) 

#Ignores the Warnings. They appears because we crop the small inner plot.
plotAvgCoalFragmSize<-p %>% 
  wrap_plots(nrow=1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")


