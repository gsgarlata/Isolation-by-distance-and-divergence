ibd.intime_parallel<-function(file,
                     path_to_patgen,
                     time_init,
                     time_end,
                     name=NULL,
                     nloci,
                     cpu_n){
  

require(ade4)
require(parallel)
require(foreach)
require(doParallel)


ibd_raw_sims<-function(df,sim_id,generation){
  tryCatch({
    temp<-subset(df,simID==sim_id)
    temp$ni<-temp$Fst/(1-temp$Fst)
    linear<-lm(temp$ni ~ log(temp$geo)) #Compute ibd at a specific timepoint
    ibd.slope<-linear$coefficients[2]
    ibd.inter<-linear$coefficients[1]
    names_id<-paste0('pop',1:nrow(temp))
    dist.gen<-structure(temp$ni, Size = nrow(temp), Labels = names_id,Diag = FALSE,Upper = FALSE,method = "user",class = "dist")
    dist.geo<-with(temp, structure(temp$geo, Size = nrow(temp), Labels = names_id,Diag = FALSE,Upper = FALSE,method = "user",class = "dist"))
    
    mantel_ibd<-mantel.randtest(dist.gen, dist.geo, nrepet = 9999)
    mantel.pval<-mantel_ibd$pvalue #record Mantel r value and p-value
    mantel.r<-mantel_ibd$obs
    res<-data.frame(time=generation,simID=sim_id,ibd_slope=ibd.slope,ibd_inter=ibd.inter,mantel_r=mantel.r,mantel_pval=mantel.pval)
    return(res)
  }, error=function(e){})
}



getIBD_info<-function(df,generation,fst_stats,sims){
  
  temp_gen<-subset(df,time==generation)

  temp_ibd<-NULL
  for(i in 1:length(fst_stats)){
    
    fst_temp<-subset(temp_gen,stat==fst_stats[i])
    
    ibd_stat<-do.call(rbind,lapply(sims,ibd_raw_sims,df=temp_gen,generation))
    
    ibd_stat$stat<-fst_stats[i]
    rownames(ibd_stat)<-NULL
    
    temp_ibd<-rbind(temp_ibd,ibd_stat)
  }
  
  return(temp_ibd)

}


path_to_data<-paste0(path_to_patgen,"/",file)

fst_raw<-read.csv(paste0(path_to_data,"/",file,"_FstStats_codom_",time_init,"_",time_end,'_',name,"_",nloci,"STR.csv"))

fst_stats<-c('Rst','Fst')
time_points<-unique(fst_raw$time)
sims<-unique(fst_raw$simID)


cl <- makePSOCKcluster(cpu_n)
clusterEvalQ(cl, .libPaths(c( .libPaths(), "/work/lchikhi/softwares/R/4.0.3/lib/")))
registerDoParallel(cl)

res<-foreach(i = 1:length(time_points), .combine = 'rbind',.packages=c('reshape2','ade4')) %dopar% {
 
  
 int_IBD<-getIBD_info(df=fst_raw,generation=time_points[i],fst_stats=fst_stats,sims=sims)
 
 return(int_IBD)
 
}

parallel::stopCluster(cl)

write.csv(as.matrix(res), file=paste0(path_to_data,"/",file,"_ibd_raw_",time_init,'_',time_end,'_',name,".csv"), row.names=FALSE)


}


