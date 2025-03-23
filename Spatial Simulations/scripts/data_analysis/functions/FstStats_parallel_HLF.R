


FstStats_parallel_HLF<-function(file,path_to_data,time_sample,name,cpu_n,n.loci,path_functions){
  require(parallel)
  require(foreach)
  require(doParallel)
  require(adegenet)
  require(reshape2)

##functions###
givePopNames<-function(df,names){
  rownames(df)<-levels(names)
  colnames(df)<-levels(names)
  return(df)
}

edit_fstmat<-function(x){
  x[upper.tri(x,diag=TRUE)]<-NA
  #res<-cbind(which(!is.na(x),arr.ind = TRUE),na.omit(as.vector(x)))
  res<-reshape2::melt(x, varnames = c('row', 'col'), na.rm = TRUE)
  return(res)
}

match_geo<-function(df, geo){
  temp<-df
  temp$combo<-paste0(temp$col,'-',temp$row)
  temp$dist<-geo[match(temp$combo,geo$combo),'geo.dist']
  end<-temp[is.na(temp$dist),]
  end[is.na(end$dist),'combo']<-paste0(end[is.na(end$dist),'row'],'-',end[is.na(end$dist),'col'])
  temp[is.na(temp$dist),'dist']<-geo[match(end$combo,geo$combo),'geo.dist']
  return(temp)
}

arrFstGeo<-function(i,df,gen,stat_fst){
  x<-df[[i]]
  x$combo<-NULL
  x$time<-gen
  x$stat<-stat_fst
  x$sim<-paste0('sim',i)
  names(x)<-c('popA','popB','Fst','geo.dist','time','stat','simID')
  return(x)
}

genind2genambig<-function(x){

  obj<-genind2df(x, sep='/')

  samples_id<-rownames(obj)
  pop_names<-obj[,'pop']
  loci<-colnames(obj)[-grep('pop',colnames(obj))]
  mygenotypes <- new("genambig", samples = samples_id,
                     loci = loci)

  loci_list<-lapply(loci, function(x) strsplit(obj[,x],'/'))

  for(i in 1:length(loci)){
    Genotypes(mygenotypes, loci = loci[i])<-loci_list[[i]]
  }

  PopInfo(mygenotypes) <- pop_names
  PopNames(mygenotypes)<- levels(pop_names)
  Ploidies(mygenotypes) <- rep(2,length(samples_id))
  Usatnts(mygenotypes) <- c(3,3)
  return(mygenotypes)

}


SubsetLoci<-function(genind.object,n.loci){
    new.genind<-genind.object[loc=c(1:n.loci)]
    return(new.genind)
  }




##functions###

func_list<-get(load(paste0(path_to_data,'/',file,'/',file,'_codom_',time_sample,'_gen_',name,'_singlepops.Rdata')))

if(n.loci==0){}else{

func_list<-ApplyFunctionToArray(func_list,SubsetLoci,n.loci)

}

pop_names<-func_list[[1]]$pop

cl <- makePSOCKcluster(cpu_n)
#clusterCall(cl, worker.init)
registerDoParallel(cl)

genambig_list<-foreach(i= 1:(length(func_list)),.packages=c('adegenet','polysat')) %dopar% {

  temp<-genind2genambig(func_list[[i]])
  return(temp)
}
parallel::stopCluster(cl)


##make frequency of alleles
cl <- makePSOCKcluster(cpu_n)
registerDoParallel(cl)
freq_list<-foreach(i= 1:(length(genambig_list)),.packages=c('polysat')) %dopar% {

  temp_freq<-simpleFreq(genambig_list[[i]])
  return(temp_freq)
}
parallel::stopCluster(cl)


list_fst_stats<-c('Fst','Rst')  #c('Fst','Gst',"Jost's D",'Rst')
plot_df<-NULL
for(u in 1:length(list_fst_stats)){
stat_fst<-list_fst_stats[u]
####START: it computes pairwise differentiation##
worker.init<- function(path_func) {
source(paste0(path_func,'/CalcDiff_Gabi.R'))
}

cl <- makePSOCKcluster(cpu_n)
clusterCall(cl, worker.init,path_functions)
registerDoParallel(cl)
stat_list<-foreach(i= 1:(length(freq_list)),.packages=c('adegenet','polysat')) %dopar% {

myStat <- calcDiffGabi(freq_list[[i]], metric = stat_fst, object = genambig_list[[i]])
 return(myStat)
}
parallel::stopCluster(cl)
####END: it computes pairwise differentiation##

#stat_named<-lapply(stat_list,givePopNames,names=pop_names)

####GEO####
int_geo_df<-gsub("pop_0.","",pop_names)
geo_df<-do.call(rbind,strsplit(int_geo_df, "-"))
geo_df<-apply(geo_df,2,as.numeric)
rownames(geo_df)<-pop_names
######
#genepop.one<-genind2genpop(genind.nested.array[[template]], process.other=TRUE) 
##computing geographic distance##
geo.dist<-dist(geo_df, method = "manhattan", diag = FALSE, upper = FALSE)
geo.dist.matrix<-as.matrix(geo.dist)
geo_data<-melt(geo.dist.matrix)[melt(upper.tri(geo.dist.matrix,diag=TRUE))$value,]
names(geo_data)<-c("PopA","PopB","geo.dist")
geo_data$combo<-paste0(geo_data$PopA,'-',geo_data$PopB)
####GEO####

stat_named<-stat_list
melt_stat<-lapply(stat_named,edit_fstmat)
stat_geo<-lapply(melt_stat,match_geo,geo=geo_data)
final_list<-lapply(1:length(stat_geo),arrFstGeo,df=stat_geo,gen=time_sample,stat_fst)

final<-do.call(rbind,final_list)

final$stat<-stat_fst
plot_df<-rbind(plot_df,final)
}

if(n.loci==0){
write.csv(as.matrix(plot_df), file=paste0(path_to_data,"/",file,'/',file,"_FstStats_codom_",time_sample,"_raw_",name,".csv"), row.names=FALSE)
}else{
write.csv(as.matrix(plot_df), file=paste0(path_to_data,"/",file,'/',file,"_FstStats_codom_",time_sample,"_raw_",name,'_',n.loci,"STR.csv"), row.names=FALSE)
}

}
