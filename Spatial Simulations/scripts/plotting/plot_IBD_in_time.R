rm(list=ls())

require(ggplot2)
require(ggpubr)


files<-c('13x13_K50_m004_hlf_mut4')
path_to_patgen = args[2] #path where the simulated scenarios are saved, typically in ./SINS/output
t_fragm<-args[3] #time of fragmentation event (in Sgarlata et al., simulations is set at 10000 for the iHL&F scenarios)
time_extract<-seq(10000,14000,500) #time points to extract for visualization
min_geo<-4 #4 in the classical 13x13 HL&F and 12 in the 39x39 HL&F
name<-c('classic') # sampling strategy: 'classic' or 'random'
id_sims<-c('m004') # a label for the simulated scenario. It is usefull if one wants to compare different scenarios in the same plot.
fst_measure = 'Rst' # 'Rst' or 'Fst'

res_final<-NULL

for(p in 1:length(files)){
  
  path_to_data<-paste0(path_to_patgen,"/",files[p])
  
  fst_file<-list.files(paste0(path_to_patgen,'/',files[p]),pattern=paste0('_FstStats_codom_[0-9]+_[0-9]+_',name))
  fst_df<-read.csv(paste0(path_to_data,"/",fst_file))
    
    stat_df<-fst_df[fst_df$stat %in% fst_stats,]
    
    time_fst<-stat_df[stat_df$time %in% time_extract,]
    
    time_fst$sampling<-name
    time_fst$scen<-id_sims[p]
    time_fst$time<-time_fst$time-t_fragm
    res_final<-rbind(res_final,time_fst)
}

res_final$ni<-res_final$Fst/(1-res_final$Fst)

res_plot<-aggregate(ni ~ geo.dist + time + stat + sampling + scen, data=res_final,mean)
res_plot$sd<-aggregate(ni ~ geo.dist + time + stat + sampling + scen,data=res_final,sd)[,'ni']
res_plot$Fst<-aggregate(Fst ~ geo.dist + time + stat + sampling + scen, data=res_final,mean)[,'Fst']
res_plot$Fst_sd<-aggregate(Fst ~ geo.dist + time + stat + sampling + scen,data=res_final,sd)[,'Fst']

res_plot$scen<-factor(res_plot$scen,levels=id_sims)

res_plot<-subset(res_plot,time>=0)

res_fst_measure<-subset(res_plot,stat==fst_measure)

time_uniq<-unique(res_fst_measure$time)

final_data<-NULL

for(u in 1:length(time_uniq)){
  
  res_fst_time<-subset(res_fst_measure, time==time_uniq[u])
  
  norm_Fst_time<-res_fst_time$Fst[res_fst_time$geo.dist==min_geo]
    
  res_fst_time$Fst<-res_fst_time$Fst-norm_Fst_time
  final_data<-rbind(final_data,res_fst_time)
}


plotIBDinTime_classic<-ggplot(data=final_data,aes(x=log(geo.dist),y=Fst/(1-Fst),colour=time))+geom_point()+
  geom_line(aes(group=factor(time)))+
  scale_color_gradient(name='Time since\n HL&F',low = 'grey39', high = '#E69F00',trans = 'reverse')+
  xlab("log(Geographical distance)")+ylab(paste0(fst_measure,'/(1 - ',fst_measure,')'))+
  theme_bw()


