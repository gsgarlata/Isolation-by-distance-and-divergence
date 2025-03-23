

    DivPerPop<-function(genind.object){
      n.pop<-seppop(genind.object)
      div.arra<-do.call("c",lapply(n.pop, function(x) lapply(summary(x)[c(3,6,7)],function(x) mean(x))))
      return(div.arra)
    }


DivPerPop_parallel<-function(data_div,time_point,cpu_n,path_to_data){
require(parallel)
  require(foreach)
  require(doParallel)
require(tidyr)

DivPerPop<-function(genind.object){
      n.pop<-seppop(genind.object)
      div.arra<-do.call("c",lapply(n.pop, function(x) lapply(summary(x)[c(3,6,7)],function(x) mean(x))))
      return(div.arra)
    }

func_list<-get(load(paste0(path_to_data,"/",data_div)))
  
worker.init <- function() {
  require(parallel)
  require(foreach)
  require(doParallel)
  DivPerPop<-function(genind.object){
    n.pop<-seppop(genind.object)
    div.arra<-do.call("c",lapply(n.pop, function(x) lapply(summary(x)[c(3,6,7)],function(x) mean(x))))
    return(div.arra)
  }
}

cl <- makePSOCKcluster(cpu_n)
clusterCall(cl, worker.init)
registerDoParallel(cl)

results<-foreach(i= 1:(length(func_list)),.combine = 'rbind',.packages=c('adegenet','dplyr','reshape2')) %dopar% {
  
  temp<-do.call(rbind,DivPerPop(func_list[[i]]))
  temp<-data.frame(id=rownames(temp),value=temp[,1])
  temp$sim<-i
  return(temp)
}
rownames(results)<-NULL
parallel::stopCluster(cl)

temp.df<-spread(results, sim, value, fill=0)
pop_infos_1st<-str_split_fixed(temp.df$id, "\\.", 2)
pop_infos_2nd<-str_split_fixed(pop_infos_1st[,2], "\\.", 2)

pop_id<-paste0(pop_infos_1st[,1],'.',pop_infos_2nd[,1])
pre_coord<-gsub("pop_0.","",pop_id)
pre_coord<-apply(do.call(rbind,strsplit(pre_coord, "-")),2,function(x) as.numeric(x))

res<-data.frame(gen=time_point,pops=pop_id,div=pop_infos_2nd[,2],x=pre_coord[,1],y=pre_coord[,2],temp.df[,-1])
names(res)[-(1:5)]<-paste0('sim',1:(length(func_list)))

return(res)

}

