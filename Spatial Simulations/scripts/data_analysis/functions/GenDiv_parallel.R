#!/bin/sh
#SBATCH -J getArp

require(parallel)
require(foreach)
require(doParallel)
require(data.table)
require(stringr)
require(adegenet)


files<-args[1]
generation<-args[2]
name<-args[3]
path_to_patgen<-args[4]
cpu_n<-args[5]
path_to_func<-args[6]

path_to_data<-paste0(path_to_patgen,"/",files)

    if(generation[1]=="all"){
      generation.rdata<-list.files(path=path_to_data,pattern=paste0("_gen_",name,"_singlepops.Rdata"))
    }else{
      generation.rdata<-list.files(path=path_to_data,pattern=paste0("_gen_",name,"_singlepops.Rdata"))
      generation.rdata<-generation.rdata[grep(paste(paste0('_',generation,'_'),collapse="|"),generation.rdata)]
    }
  
  gen.old<-gsub(".*_codom_","",generation.rdata)
  gen<-sort(as.numeric(gsub("_gen.*","",gen.old)))


    path_to_func<-'/home/lchikhi/work/Gabriele/Spatial_Sims/SINS_Analysis/Functions/HLF'
    source(paste0(path_to_func,"/Functions_GenDiv_parallel.R"))
    DivPerPop<-function(genind.object){
      n.pop<-seppop(genind.object)
      div.arra<-do.call("c",lapply(n.pop, function(x) lapply(summary(x)[c(3,6,7)],function(x) mean(x))))
      return(div.arra)
    }
 
  
 final_res<-NULL
  
for(i in 1:length(generation.rdata)){

    div_temp<-DivPerPop_parallel(generation.rdata[i],time_point=gen[i],cpu_n,path_to_data)
final_res<-rbind(final_res,div_temp) 
 }

  
  time_results<-final_res[order(final_res$gen),]

if(generation=='all'){

init_plot<-gen[1]
end_plot<-gen[length(gen)]

write.csv(as.matrix(time_results),file=paste0(path_to_data,"/",files,'_diversity_raw_',init_plot,'_',end_plot,"_",name,".csv"),row.names=FALSE,quote=FALSE)

}else{

write.csv(as.matrix(time_results),file=paste0(path_to_data,"/",files,'_diversity_raw_',generation,"_",name,".csv"),row.names=FALSE,quote=FALSE)

}
