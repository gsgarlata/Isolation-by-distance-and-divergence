#This script is used to plot pairwise Fst over time for deme pairs at different geographical distance.

require(ggplot2)

files<-c('13x13_K50_m004_hlf_mut4') #name of the simulated scenario
path_to_patgen = args[2] #path where the simulated scenarios are saved, typically in ./SINS/output
t_fragm<-args[3] #time of fragmentation event (in Sgarlata et al., simulations is set at 10000 for the iHL&F scenarios)
time_extract<-seq(10000,14000,500) #time points to extract for visualization
sampling<-'classic' # 'classic' or 'random'
fst_stats<-'Rst'# 'Rst' or 'Fst'
id_sims<-c('m004') # a label for the simulated scenario. It is usefull if one wants to compare different scenarios in the same plot.


res_final<-NULL

for(p in 1:length(files)){
  
  path_to_data<-paste0(path_to_patgen,"/",files[p])
  
  fst_file<-list.files(paste0(path_to_patgen,'/',files[p]),pattern=paste0('_FstStats_codom_[0-9]+_[0-9]+_',sampling))
  fst_df<-read.csv(paste0(path_to_data,"/",fst_file))
  
  stat_df<-subset(fst_df,stat==fst_stats)
  
  stat_df$sampling<-sampling
  stat_df$stat<-fst_stats
  stat_df$scen<-id_sims[p]
  stat_df$time<-stat_df$time-t_fragm
  res_final<-rbind(res_final,stat_df)
}

res_final$ni<-res_final$Fst/(1-res_final$Fst)


res_plot<-aggregate(Fst ~  geo.dist + time + stat + scen, data=res_final,mean)
res_plot$sd<-aggregate(Fst ~ geo.dist + time + stat + scen,data=res_final,sd)[,'Fst']

res_plot$scen<-factor(res_plot$scen,levels=id_sims)

if(is.null(time_extract)){
  
res_tot<-res_plot

}else{
time_extract<-time_extract-t_fragm
res_tot<-res_plot[res_plot$time %in% time_extract,]
}

###This plots pairwise Fst by geographical distance over time####
plot_Fst<-ggplot(data=res_tot,aes(x=time,y=Fst,colour=factor(geo.dist)))+geom_point()+
  geom_ribbon(aes(ymin=Fst - sd, ymax=Fst + sd,fill = factor(geo.dist)),colour = NA,alpha= 0.1,show.legend=FALSE)+
  geom_line(aes(colour=factor(geo.dist)))+
  geom_vline(xintercept = 0)+
  scale_color_discrete(name='')+
  theme_bw()+xlab("Time")+ylab(fst_stats)+facet_wrap(vars(scen),scales='free',nrow=1)



