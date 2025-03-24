rm(list=ls())

require(reshape2)
require(ggplot2)

sim_files<-c('13x13_K50_m004_hlf_mut4')
path_to_patgen = args[2] #path where the simulated scenarios are saved, typically in ./SINS/output
id_sims<-c('0.04') # a label for the simulated scenario. It is usefull if one wants to compare different scenarios in the same plot.
t_fragm<-10000 #time of fragmentation event (in Sgarlata et al., simulations is set at 10000 for the iHL&F scenarios)
sampling<-c('classic') # sampling strategy: 'classic' or 'random'
div_stat<-'Hexp' #genetic diversity statistics: 'Hexp' or 'Hobs' or 'loc.n.all'

final<-NULL

for(i in 1:length(sim_files)){
  
  sim_file<-sim_files[i]
  path_to_data<-paste0(path_to_patgen,"/",sim_file)
  
  div_file<-list.files(paste0(path_to_data),pattern=paste0('_diversity_raw_[0-9]+_[0-9]+_',sampling))
  div.array<-read.csv(paste0(path_to_data,"/",div_file))
  
  res.na<-grep('NA',names(div.array))
  if(length(res.na)==0){}else{
    div.array<-div.array[,-grep('NA',names(div.array))]
  }
  div.df<-subset(div.array,div==div_stat)
  sims<-c(1:length(grep("sim",colnames(div.df))))
  cnames<-NULL
  for(j in 1:length(sims)){
    cnames[[j]]<-paste("sim",j,sep="")
  }
  div.df.new<-melt(div.df,measure.vars=cnames)
  work.df<-aggregate(value ~ gen + pops + div + x + y,data = div.df.new,"mean")
  work.df$sd<-aggregate(value ~ gen + pops + div + x + y,data = div.df.new,"sd")[,"value"]
  names(work.df)[1]<-'time'
  work.df$type<-work.df$pops
  work.df$scen<-id_sims[i]
  work.df$sampling<-sampling
  final<-rbind(final,work.df)
  
}


final$time<-final$time-t_fragm

res_mean<-aggregate(value ~ time + div + scen + sampling,data = final,mean)
res_mean$sd<-aggregate(value ~ time + div + scen + sampling,data = final,sd)[,"value"]

res_mean$scen<-paste0('m: ',res_mean$scen)

if(div_stat=='Hobs'){
  label_y<-expression(H["o"])
}else{
  label_y<-expression(H["e"])
}


plot_true_val<-ggplot(data=res_mean,aes(x=time,y=value,group=factor(sampling),colour=factor(scen)))+
  geom_line(aes(linetype=factor(sampling)),size=1)+
  scale_linetype_discrete(name='Sampling\nScheme')+
  scale_colour_manual(name='Migration\n rate',values=c("#999999", "#E69F00", "#56B4E9","red"))+
  geom_vline(xintercept=0)+
  facet_wrap(vars(scen),scale='fixed',nrow = 1)+theme_bw()+
  labs(y=label_y, x='Time')+
  theme(strip.background =element_rect(fill="white"),strip.text.x = element_text(size = 12))


res_mean$combo<-paste0(res_mean$scen,'_',res_mean$sampling)

fina_norm<-NULL

for(migr in unique(res_mean$combo)){
  temp_m<-subset(res_mean,combo==migr)
  coord<-which(temp_m$time==0)
  temp_m$value<-temp_m$value-temp_m$value[coord]
  fina_norm<-rbind(fina_norm,temp_m)
  
}


plot_norm<-ggplot(data=fina_norm,aes(x=time,y=value,group=factor(combo),colour=factor(scen)))+
  geom_line(aes(linetype=factor(sampling)),size=1)+
  scale_linetype_discrete(name='Sampling\nScheme')+
  scale_colour_manual(name='Migration\n rate',values=c("#999999", "#E69F00", "#56B4E9","red"))+
  geom_vline(xintercept=0)+
  facet_wrap(vars(sampling),scale='fixed')+
  theme_bw()+
  labs(y=label_y, x='Time')+
  theme(strip.background =element_rect(fill="white"),strip.text.x = element_text(size = 12))



