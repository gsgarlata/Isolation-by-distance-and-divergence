rm(list=ls())
require(ggplot2)


files<-c('13x13_K50_m004_hlf_mut4')
path_data = args[2] #path where the simulated scenarios are saved, typically in ./SINS/output
time_init<-rep(9000,length(files)) #first sampling time point
time_end<-rep(14000,length(files)) #last sampling time point
t_fragm<-rep(10000,length(files)) #time of fragmentation event (in Sgarlata et al., simulations is set at 10000 for the iHL&F scenarios)
sampling<-c('classic') # sampling strategy: 'classic' or 'random'
id_sims<-c('0.04') # a label for the simulated scenario. It is usefull if one wants to compare different scenarios in the same plot.


final_plot_avg<-NULL

for(u in 1:length(files)){
  
  input_file<-files[u]
  path_input<-paste0(path_data,'/',input_file)
  
  ibd_file<-list.files(path=path_input,pattern=paste0('_ibd_raw_[0-9]+_[0-9]+_',sampling))
  ibd<-read.csv(paste0(path_input,"/",ibd_file))
  
  ibd_avg<-aggregate(cbind(ibd_slope,mantel_pval) ~ time + stat, data=ibd,mean)
  sd_data<-aggregate(cbind(ibd_slope,mantel_pval) ~ time + stat, data=ibd,sd)
  ibd_avg$sd_slope<-sd_data[,'ibd_slope']
  ibd_avg$sd_pval<-sd_data[,'mantel_pval']
  ibd_avg$scen<-id_sims[u]
  ibd_avg$time<-ibd_avg$time-t_fragm[u]
  ibd_avg$time_fragm<-0
  ibd_avg$sampling<-sampling
  final_plot_avg<-rbind(final_plot_avg,ibd_avg)
  
}

Fst_avg_plot<-subset(final_plot_avg,stat=='Fst')
Rst_avg_plot<-subset(final_plot_avg,stat=='Rst')

Fst_avg_plot$combo<-paste0(Fst_avg_plot$scen,'_',Fst_avg_plot$sampling)

Fst_pval<-ggplot(Fst_avg_plot, aes(x=time,y=mantel_pval,colour=factor(scen)))+
  geom_hline(yintercept = 0.05,colour='black',linetype='dashed')+
  geom_ribbon(aes(ymin=mantel_pval-sd_pval,ymax=mantel_pval+sd_pval,fill=factor(scen)),alpha=.1,colour=NA,show.legend = FALSE)+
  geom_point()+geom_line()+
  geom_vline(data=final_plot_avg,mapping=aes(xintercept = time_fragm))+
  facet_wrap(vars(sampling))+
  scale_colour_manual(name='Migration\n rate',values=c("#999999", "#E69F00", "#56B4E9","red"))+
  scale_fill_manual(name='',values=c("#999999", "#E69F00", "#56B4E9","red"))+
  ylab('IBD p-value')+xlab('Time since HL&F')+theme_bw()+xlim(-1000,3000)+
  theme(strip.background =element_rect(fill="white"),strip.text.x = element_text(size = 12))

Rst_pval<-ggplot(Rst_avg_plot, aes(x=time,y=mantel_pval,colour=factor(scen)))+
  geom_hline(yintercept = 0.05,colour='black',linetype='dashed')+
  geom_ribbon(aes(ymin=mantel_pval-sd_pval,ymax=mantel_pval+sd_pval,fill=factor(scen)),alpha=.1,colour=NA,show.legend = FALSE)+
  geom_point()+geom_line()+
  geom_vline(data=final_plot_avg,mapping=aes(xintercept = time_fragm))+
  facet_wrap(vars(sampling))+
  scale_colour_manual(name='Migration\n rate',values=c("#999999", "#E69F00", "#56B4E9","red"))+
  scale_fill_manual(name='',values=c("#999999", "#E69F00", "#56B4E9","red"))+
  ylab('IBD p-value')+xlab('Time since HL&F')+theme_bw()+xlim(-1000,3000)+
  theme(strip.background =element_rect(fill="white"),strip.text.x = element_text(size = 12))

