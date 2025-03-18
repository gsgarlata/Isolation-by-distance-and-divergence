#This script allows to simulate microsatellite data for a panmictic population at mutation-drift equilibrium.
#This can be useful in order to have some initial level of genetic diversity when performing genetic simulations
#in SINS.


#pop.size: total number of individuals (we consider a panmictic population size equal to the full population size in the stepping-stone model)

#sample.size: how many individuals do you want to sample from the simulated panmictic population.

#num.loci: n° of microsatellite loci to simulate.

#mut.rate: mutation rate for each microsatellite locus.

#exec: name of the fastsimcoal executer (for example "./fsc26_mac64" if working on MacOS)

#path_to_exec: path where the fastsimcoal executer is located.


##START: Function code##
Init_freq<-function(pop.size,sample.size,locus.type="msat",num.loci,mut.rate,exec,path_to_exec){

  require("strataG")
  require("adegenet")
  current_dir<-getwd()
  setwd(path_to_exec)
  
  markers<-seq(1:num.loci)
  
  if(num.loci>10){
    for(p in 1:9){
      markers[[p]]<-paste("00",markers[[p]],sep="")
    }
  }
    if(num.loci>99){
    for(s in 10:99){
      markers[[s]]<-paste("0",markers[[s]],sep="")
      
    }
    }
  
  for(k in 1:num.loci){
    markers[[k]]<-paste("Locus_",markers[[k]],sep="")
  }
  
  #It defines the panmictic population size and sampling size
  pop.info<-fscDeme(deme.size = pop.size, sample.size = sample.size)
  #It defines the genetics of the panmictic population
  loc.params<-fscLocusParams(locus.type = locus.type, num.loci = num.loci,mut.rate = mut.rate)
  #It simulate the panmictic population using a wrapper of fastsimcoal 
  gen.data.gtype<-fastsimcoal(pop.info, loc.params,exec = exec)
  #It convert a "gtypes" object to "genind"
  gen.data.genind<-gtypes2genind(gen.data.gtype)
  #It converts genind object to genepop object
  gen.data.pop<-genind2genpop(gen.data.genind)
  #It computes allele frequency for each locus.
  gen.data.freq<-makefreq(gen.data.pop)
  
  setwd(current_dir)
  #It creates the matrix for the autosomal information required by the SINS program.
  genotype<-matrix(nrow=num.loci*2+1,ncol=2)
  #It writes the number of autosomal loci.
  genotype[1,]<-c("nbAutosomes",num.loci)
  #It creates the matrix for the autosomal mutation rate.
  mut_rate<-matrix(nrow=num.loci,ncol=2)
  #It is usefull for the following loop.
  t<-1
  
  for(i in 1:num.loci){
    
    #It allows to not have abbreviations for decimal numbers.
    options("scipen"=10)
    #It lists the different alleles of a locus.
    alleles<-grep(paste0(markers[i],'\\.'), colnames(gen.data.freq), value=TRUE) 
    #It creates the matrix for the single file about the autosomal information.
    all_freq<-matrix(nrow=1,ncol=2)
    #It writes the number of alleles for a locus.
    all_freq[1,]<-c("nbAlleles",length(alleles))
    #It extracts the allele frequencies for each locus.
    ty<-unlist(lapply(alleles, function(x) gen.data.freq[,which(colnames(gen.data.freq)==x)]))
    #It extracts the names of the alleles.
    size<-gsub(paste0(markers[i],"\\."),"",alleles)
    
    for(j in 1:length(size)){
    size[[j]]<-sub("........","",size[[j]])
    }
    
    #It binds allele frequencies and the corresponding allele type.
    u<-cbind(ty,size) 
    #It links the info on the number of alleles, allele frequencies and allele type.
    STR_marker<-rbind(all_freq,u)
    #It writes the autosomal allele file.
    write.table(STR_marker, paste0("allelesA", i,".txt"), quote=FALSE, sep=" ", row.names=FALSE, col.names=FALSE)
    #It writes in the genotype matrix the length info for each autosomal marker.
    genotype[1+t,]<-c(paste0("A",i,"length"),1)
    #It specifies the type of genetic marker.
    genotype[1+t+1,]<-c(paste0("typeA",i),"microsat")
    #It is useful for the loop.
    t<-t+2
    #It defines the mutation rate for each autosomal locus.
    mut_rate[i,]<-c(paste0("A",i,"mutationRate"), mut.rate)
    
  }
  
  #It takes a genotype file for which we don't change info on X,Y, mtDNA markers.
  template<-as.data.frame(read.table(file="genotype.txt", header=FALSE, sep = " ",as.is=TRUE, check.names=FALSE))
  #It tells in which row we have "Xmutationrate"
  temp_mut_rate_a<-which(template=="XmutationRate")
  #It tells in which row we have "mtDNAmutationrate"
  temp_mut_rate_b<-which(template=="mtDNAmutationRate")
  #It extract the lines of the X,Y, mtDNA marker type and length that we don't change.
  begin_geno<-template[1:6,]
  #It extract the lines of the X,Y, mtDNA mutation rate that we don't change.
  begin_mut<-template[temp_mut_rate_a:temp_mut_rate_b,]
  
  input<-rbind(begin_geno,genotype,begin_mut,mut_rate) #It binds everything related to the "genotype" file.
  write.table(input, "genotype_new.txt", quote=FALSE, sep=" ", row.names=FALSE, col.names=FALSE) #It writes a new genotype file.
  
}
##END: Function code##


