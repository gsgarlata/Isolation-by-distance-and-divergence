

#sampling: can be either 'random' (individuals are sampled randomly one per deme within the same fragment) or 'classic' meaning the usual single central deme.
#n_samples: represents the n° of individuals to sample in total for each patch.- Important only if 'random' sampling is specified.

ExtractPops<-function(files,
                      file.path,
                      generation="all",
                      path.patch,
                      patch.file,
                      name="classic",
                      n_samples){

require(adegenet)  

      isNested <- function(l){
                      stopifnot(is.list(l))
        for (i in l) {
                          if (is.list(i)) return(TRUE)
                }
                return(FALSE)
                      } #This function checks whether there are multiple simulation


  path_to_data<-paste0(file.path,"/",files)
  
  if(generation=="all"){
  rdata<-list.files(path=path_to_data,pattern="_gen.Rdata")
  gen.old<-gsub(".*_codom_","",rdata)
  gen<-sort(as.numeric(gsub("_gen.*","",gen.old)))
  generation.rdata<-paste0(files,"_codom_",gen,"_gen.Rdata")
  }else{
    rdata<-list.files(path=path_to_data,pattern=paste0(generation,"_gen.Rdata"))
    generation.rdata<-rdata
  }
  
  
  for(d in 1:length(generation.rdata)){
    
    
genind.nested.array<-get(load(paste0(path_to_data,"/",generation.rdata[d])))

if(name=="random"){
  
patch<-read.csv(paste0(path.patch,"/",patch.file,"_patches.csv"),sep=";")
xmin = patch$xmin
xmax = patch$xmax
ymin = patch$ymin
ymax = patch$ymax
patch_names = patch$patch.name
patch.type=patch$Patch.type
  

# obj = genind.nested.array[[1]][[1]]
# 
# 
# xmin_val = xmin[1]
# xmax_val = xmax[1]
# ymin_val = ymin[1]
# ymax_val = ymax[1]
  
extrPops_random<-function(obj,xmin_val,xmax_val,ymin_val,ymax_val){
    require(data.table)
    x_list<-xmin_val:xmax_val
    y_list<-ymin_val:ymax_val
    pop.name<-paste0("pop_0.",x_list,"-")
    pops<-as.vector(sapply(y_list, function(x) hu<-paste0(pop.name,x)))
    df_data<-genind2df(obj, sep="/", oneColPerAll = FALSE) #Convert STRUCTURE data file in a "df" file with one row per individual
    dfgeninPops<-lapply(1:length(pops), function(x) df_data[which(df_data$pop==pops[x]),])
    dfgenin<-do.call(rbind, unname(dfgeninPops))
    
    init_pops<-pops
    tot.df<-NULL
    for(l in 1:n_samples){
      
      temp.pop<-sample(init_pops,1)
      temp.df<-subset(dfgenin,pop==temp.pop)
      
      end.df<-temp.df[sample(nrow(temp.df),1),]
      tot.df<-rbind(tot.df,end.df)
      
      init_pops<-init_pops[!init_pops %in% temp.pop]
      if(length(init_pops)==0){
        init_pops<-pops
      }
    }
    
    
    locus_names<-gsub("\\.", "_", names(tot.df)[-c(1)])
    mat_data<-as.matrix(tot.df)
    gt<-as.vector(tot.df$pop)
    
    pop_coords<-gsub("pop_0.","",gt)

    new.coord<-c(median(x_list),median(y_list))
    tryCatch({
    tot.df$pop<-paste0("pop_0.",new.coord[1],"-",new.coord[2])
    return(tot.df)  }, error=function(e){})
  }
  

ExtractPatches_random<-function(obj,xmin,xmax,ymin,ymax){
  
    onegen<-lapply(1:length(xmin),function(x) extrPops_random(obj,xmin[x],xmax[x],ymin[x],ymax[x]))
    df.save<-do.call("rbind",onegen)
    locus_names<-gsub("\\.", "_", names(df.save)[-c(1)])
    genin<-df2genind(df.save[,-1],sep="/",ncode=3,ind.names=rownames(df.save),loc.names=locus_names,pop=as.vector(df.save$pop),ploidy = 2,NA.char="NA.NA",type="codom")
    
    gt<-as.vector(genin@pop)
    pop_coords<-gsub("pop_0.","",gt)
    genin@other$xy<-do.call(rbind,strsplit(pop_coords, "-"))
    genin@other$xy<-apply(genin@other$xy,2,as.numeric)
    genin@other$generation<-obj@other$generation
    genin@other$simulation.ID<-obj@other$simulation.ID
    
    return(genin)
  }
  
  
  if(isNested(genind.nested.array)==FALSE){
    
    func_list<-ExtractPatches_random(genind.nested.array[[1]],xmin,xmax,ymin,ymax)
    
  }else{
    
    func_list<-lapply(genind.nested.array,function(x) ExtractPatches_random(x[[1]],xmin,xmax,ymin,ymax))
  }
  
  
  filename<-gsub("_gen.Rdata","",generation.rdata[d])
  
  save(func_list,file=paste0(path_to_data,"/",filename,"_gen_",name,"_singlepops.Rdata"))
  
  
  
}else{

patch<-read.csv(paste0(path.patch,"/",patch.file,".csv"),sep=",")
x = patch$x
y = patch$y
patch_names = patch$patch.name
patch.type=patch$Patch.type

  
extrPops<-function(obj,x,y){
require(data.table)
 
pops<-paste0("pop_0.",x,"-",y)
print(pops)
df_data<-genind2df(obj, sep="/", oneColPerAll = FALSE) #Convert STRUCTURE data file in a "df" file with one row per individual
dfgeninPops<-lapply(1:length(pops), function(x) df_data[which(df_data$pop==pops[x]),])
dfgenin<-do.call(rbind, unname(dfgeninPops))
locus_names<-gsub("\\.", "_", names(dfgenin)[-c(1)])
dfgenin_mat<-as.matrix(dfgenin)
genin<-df2genind(dfgenin_mat[,-1],sep="/",ncode=3,ind.names=rownames(dfgenin_mat),loc.names=locus_names,pop=as.vector(dfgenin$pop),ploidy = 2,NA.char="NA.NA",type="codom")

gt<-as.vector(genin@pop)
pop_coords<-gsub("pop_0.","",gt)
genin@other$xy<-do.call(rbind,strsplit(pop_coords, "-"))
genin@other$xy<-apply(genin@other$xy,2,as.numeric)
genin@other$generation<-obj@other$generation
genin@other$simulation.ID<-obj@other$simulation.ID
return(genin)
}

ApplyExtrPopsToArray<-function(nested_list,x,y){
  isNested <- function(l){
    stopifnot(is.list(l))
    for (i in l) {
      if (is.list(i)) return(TRUE)
    }
    return(FALSE)
  }
  
  if(isNested(nested_list)==FALSE){
      func_list<-extrPops(nested_list[[1]],x,y)
    }else{
      func_list<-lapply(nested_list,function(o) extrPops(o[[1]],x,y))
    }
  
  return(func_list)
}


func_list<-ApplyExtrPopsToArray(genind.nested.array,x,y)

filename<-gsub("_gen.Rdata","",generation.rdata[d])

save(func_list,file=paste0(path_to_data,"/",filename,"_gen_",name,"_singlepops.Rdata"))
}

}

  
} 
  

