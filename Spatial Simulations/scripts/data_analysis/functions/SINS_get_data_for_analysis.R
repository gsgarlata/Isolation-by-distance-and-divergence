

###########################################################################
###########################################################################
###########################################################################
# Functions in this file are used to get data from SINS into R so that
# we can run analysis on this data
###########################################################################
###########################################################################
###########################################################################






# Use code below only if development version of adegenet is necessary
#if(!require("devtools")){
#  install.packages("devtools", dependencies = TRUE)
#  install_github("thibautjombart/adegenet")
#  library("adegenet")
#  }


if (suppressPackageStartupMessages(!require("adegenet")))
{install.packages("adegenet", dependencies = TRUE)}
if (suppressPackageStartupMessages(!require("ade4")))
{install.packages("ade4", dependencies = TRUE)}
if (suppressPackageStartupMessages(!require("hierfstat")))
{install.packages("hierfstat", dependencies = TRUE)}
if (suppressPackageStartupMessages(!require("pegas")))
{install.packages("pegas", dependencies = TRUE)}
if (suppressPackageStartupMessages(!require("ggplot2")))
{install.packages("ggplot2", dependencies = TRUE)}
if (suppressPackageStartupMessages(!require("reshape2")))
{install.packages("reshape2", dependencies = TRUE)}
if (suppressPackageStartupMessages(!require("gridExtra")))
{install.packages("gridExtra", dependencies = TRUE)}
#if(!require("cowplot")) install.packages("cowplot", dependencies = TRUE)
#if(!require("plotly")) install.packages("plotly", dependencies = TRUE)
if (suppressPackageStartupMessages(!require("Hmisc")))
{install.packages("Hmisc", dependencies = TRUE)}
if (suppressPackageStartupMessages(!require("plyr")))
{install.packages("plyr", dependencies = TRUE)}
#if (suppressPackageStartupMessages(!require("PopGenReport")))
# {install.packages("PopGenReport", dependencies = TRUE)}
if (suppressPackageStartupMessages(!require("data.table")))
{install.packages("data.table", dependencies = TRUE)}
if (suppressPackageStartupMessages(!require("StAMPP")))
{install.packages("StAMPP", dependencies = TRUE)}


###########################################################################
###########################################################################
#++++SECTION 1 - FUNCTIONS TO GET DATA
###########################################################################
###########################################################################

################
#SET NUMBER OF SIMULATIONS
################
#+++++++++++++++
Set_number_of_simulations = function(number_of_simulations) {
  return(number_of_simulations)
}


################
#SET MARKERS
################
#+++++++++++++++
Set_markers_list = function(...) {
  # Make list of the given markers (e.g. "A1","A2","A3", ...)
  the_markers_list = NULL
  the_markers_list = c(...)
  return(the_markers_list)
}

################
#SET GENERATIONS
################
#+++++++++++++++
Set_generations_list = function(initialTimeStep,
                                finalTimeStep,
                                timeStepIncrement) {
  # Make sequence of all the generations that we want to study (e.g. from gen 1 to 100 by increments of 5)
  the_generations_list = NULL
  the_generations_list = seq(initialTimeStep, finalTimeStep, timeStepIncrement)
  return(the_generations_list)
}

################
#SET LAYERS
################
#+++++++++++++++
Set_layers_list = function(...) {
  # Make list of the given layers (e.g. "layer0", "layer1")
  the_layers_list = NULL
  the_layers_list = c(...)
  return(the_layers_list)
}

################
#GET FILE NAMES
################
#+++++++++++++++
Get_file_names = function(generations_list, layers_list) {
  # Build the names of the files that we want to analyse
  the_file_names = NULL
  
  for (j in generations_list) {
    the_file_names = rbind(the_file_names, paste(j, "_", layers_list, "_", sep = ""))
    
  }
  return(the_file_names)
}

################
#READ DATA
################
#+++++++++++++++

Read_data_func = function(name_of_file,
                          markers_list,
                          is_marker_diploid = TRUE,
                          simulation_ID = 0,
                          subind) {
  #name_of_file = "11000_layer0_"
  #markers_list = paste0("A",1:30)
  #is_marker_diploid=T
  #simulation_ID=1
  #setwd(dir = paste(the_path_to_data,"Adegenet_sim_1","/",sep = ""))
  
  if (FALSE) {
    #***fread is much faster than an optimized read.table***
    #read first 5 rows of the table
    #tab5rows <- read.table(file = paste(name_of_file,markers_list[1],sep = ""), header = FALSE, nrows = 100,row.names = 1,na.strings = "NA")
    
    #get the classes of the columns of the table
    #tabClasses <- sapply(tab5rows, class)
  }
  data_f <-
    fread(
      input = paste(name_of_file, markers_list[1], sep = ""),
      encoding = "UTF-8",
      header = FALSE,
      na.strings = "NA",
      data.table = F
    )
  
  if(is.null(subind)==FALSE){
    
    subpops<-unique(data_f$V5)

    data_sub<-NA
    errorSampl<-NA
    
    for(t in 1:length(subpops)){
      
      data_temp<-data_f[which(data_f$V5==subpops[t]),]
      
      if(nrow(data_temp)<subind){
        
        newdata<-data_temp
      }else{
        newdata<-data_temp[sample(1:nrow(data_temp),subind),]
      }
      
      data_sub<-rbind(data_sub,newdata)
      
    }
    
    tot_ind<-data_sub[-1,]$V1
    
  }else{}
  
  
  # get the number of rows that the table will have, use it later to speed up process
  # not sure if it is that much faster, consider erasing it
  
  if(is.null(subind)==FALSE){
    sizeOfTable = length(tot_ind)
  }else{
    sizeOfTable = nrow(data_f)
  }
  
  markerCondition = markers_list == "MT" | markers_list == "Y"
  
  current_generation = NULL
  current_generation = strsplit(x = name_of_file, split = "_")[[1]][1]
  
  data_f = NULL
  marker_data.haplo = NULL
  marker_data.diplo = NULL
  
  a_file_name = NULL
  
  #bind_markers.diplo = NULL
  bind_markers.diplo = matrix(NA, nrow = sizeOfTable, ncol = length(markers_list))
  #bind_markers.haplo = NULL
  bind_markers.haplo = matrix(NA, nrow = sizeOfTable, ncol = length(markers_list))
  
  
  ##KEY##
  for (marker in 1:length(markers_list)) {
    a_file_name = paste(name_of_file, markers_list[marker], sep = "")
    
    #read the rest (or a bigger part) of the table knowing the class of the columns speeds up the reading process
    #data_file <- read.table(a_file_name, header = FALSE, colClasses = tabClasses, nrows = sizeOfTable, row.names = 1,na.strings = "NA")
    #fread is much faster than an optimized read.table
    data_fm <-
      fread(
        input = a_file_name,
        header = FALSE,
        encoding = "UTF-8",
        na.strings = "NA",
        data.table = F
      )
    
    if(is.null(subind)==FALSE){
      
      data_file<-subset(data_fm, data_fm$V1 %in% tot_ind) 
      
    }else{data_file<-data_fm}
    
    #if marker is Y or MTdna then it is haploid data, save that markers col (instead of 2 cols for diploid)
    if (markerCondition[marker]) {
      marker_data.haplo = as.matrix(data_file$V2)
      #bind_markers.haplo = cbind(bind_markers.haplo,marker_data.haplo)
      bind_markers.haplo[, marker] = marker_data.haplo
    } else{
      #marker_data.diplo = as.matrix(paste(data_file$V2,data_file$V3,sep = "/"))
      
      marker_data.diplo = as.vector(data_file$V2)
      
      
      #length(marker_data.diplo)
      #bind_markers.diplo = cbind(bind_markers.diplo,marker_data.diplo)
      bind_markers.diplo[, marker] = marker_data.diplo
    }
  }
  
  
  
  if (is_marker_diploid) {
    # note ncode set to 3 - alleles are considered to be coded by three characters
    markers.diplo_df2genind = df2genind(
      X = bind_markers.diplo,
      ploidy = 2,
      ncode = 3,
      ind.names = data_file$V1,
      sep = "/"
    )
    # define pops
    markers.diplo_df2genind@pop = factor(data_file$V5)
    # define coords
    xy_coords = cbind(data_file$V3, data_file$V4)
    markers.diplo_df2genind@other$xy = xy_coords
    # define gen
    markers.diplo_df2genind@other$generation = current_generation
    # define simulation ID
    markers.diplo_df2genind@other$simulation.ID = simulation_ID
    # return genind object
    return(markers.diplo_df2genind)
  } else{
    markers.haplo_df2genind = df2genind(
      X = bind_markers.haplo,
      ploidy = 1,
      ncode = 3,
      ind.names = data_file$V1,
      sep = "/",
      NA.char = "NA"
    )
    # define pops
    markers.haplo_df2genind@pop = factor(data_file$V5)
    # define coords
    xy_coords = cbind(data_file$V3, data_file$V4)
    markers.haplo_df2genind@other$xy = xy_coords
    # define gen
    markers.haplo_df2genind@other$generation = current_generation
    # define simulation ID
    markers.haplo_df2genind@other$simulation.ID = simulation_ID
    # return genind object
    return(markers.haplo_df2genind)
  }
  

}

################
#GET DATA FROM FILES FROM A SINGLE SIMULATION
################
#+++++++++++++++
Get_data_from_files_single_sim = function(the_file_names,
                                          set_markers,
                                          generations_list,
                                          path_to_data,
                                          simulation_number,
                                          is_marker_diploid = TRUE,
                                          subind) {
  # Gets data from the files with the parameters gotten from the previous functions
  
  setwd(paste(path_to_data, "Adegenet_sim_", simulation_number, "/", sep = ""))
  
  the_data_list = NULL
  
  for (i in the_file_names) {
    the_data_list = c(
      the_data_list,
      Read_data_func(
        i,
        set_markers,
        is_marker_diploid = is_marker_diploid,
        simulation_ID = simulation_number,
        subind)
    )
  }
  #set names of elements in list to the generations that they belong to
  names(the_data_list) = generations_list
  
  
  return(the_data_list)
}

################
#GET DATA FROM FILES FROM MULTIPLE SIMULATIONS
################
#+++++++++++++++
#marker_type: codom; mtDNA; Y_chrom. 
#path_to_genind: If this argument is specified, the function will save the genind object into the path_to_genind folder.

Get_data_from_files_multi_sim = function(the_file_names,
                                         set_markers,
                                         generations_list,
                                         path_to_data,
                                         number_of_simulations,
                                         is_marker_diploid, 
                                         marker_type,   # marker_type should be set to "codom" for autosomes
                                         path_to_genind=NULL,
                                         subind) {
  # Gets data from the files with the parameters gotten from the previous functions
  # Returns list of lists of genind objects
  
  ####test data
  #set_markers = c("A1")
  #generations_list = seq(0,100,10)
  #number_of_simulations=1
  #the_file_names= Get_file_names(generations_list,"layer0")
  #path_to_data = the_path_to_data
  #is_marker_diploid = TRUE
  ####
  
  all_sim = vector(mode = "list", length = number_of_simulations)
  
  
  for (j in 1:number_of_simulations) {
    setwd(paste(path_to_data, "/", "Adegenet_sim_", j, "/", sep = ""))
    
    the_data_list = vector(mode = "list", length = length(the_file_names))
    for (i in 1:length(the_file_names)) {
      #the_data_list = c(the_data_list,Read_data_func(i, set_markers,is_marker_diploid=is_marker_diploid, simulation_ID = j))
      the_data_list[[i]] <-
        Read_data_func(
          the_file_names[i],
          set_markers,
          is_marker_diploid = is_marker_diploid,
          simulation_ID = j,
          subind)
      #the_data_list = lapply(the_file_names,FUN = Read_data_func, name_of_file = the_file_names, markers_list = set_markers, is_marker_diploid = T, simulation_ID = j)
      
    }
    #set names of elements in list to the generations that they belong to
    names(the_data_list) = generations_list
    
    all_sim[[j]] = the_data_list
  }
  
  ###test
  #str(list(the_data_list,the_data_list),max.level = 1)
  #str(all_sim, max.level = 1)
  #dim(all_sim)
  #class(all_sim)
  #all_sim[[2]][[5]]@other$generation
  ###
  if(is.character(path_to_genind)==TRUE){
    if(length(generations_list)>1){  
      save(all_sim,file=paste0(path_to_genind,"/",nameOfSimulation,"_",marker_type,"all_generations.Rdata"))
    }else{
      save(all_sim,file=paste0(path_to_genind,"/",nameOfSimulation,"_",marker_type,"_",generations_list,"_gen.Rdata"))
    }
  } else{}
  return(all_sim)
}

#*****************************
# GET SEQUENCE OR SNP DATA FROM FILES
#*****************************
#*****************************
Get_seqSNP_data = function(number_of_simulations,
                           path_to_data,
                           generations_list,
                           layer_ID,
                           marker_type,
                           path_to_genind=NULL
){
  
  #number_of_simulations = 2
  #path_to_data = asd
  #generations_list = seq(10,50,10)
  #layer_ID="layer0"
  
  raw_snp_data = vector(mode = "list", length = number_of_simulations)
  
  for(simID in 1:number_of_simulations){
    
    setwd(paste(path_to_data, "/Adegenet_sim_", simID, "/", sep = ""))
    
    raw_data_per_sim = vector(mode = "list", length = length(generations_list))
    print(paste("Simulation ",simID," of ",number_of_simulations, sep=""))
    
    for(generationID in 1:length(generations_list)){
      
      currentGeneration = generations_list[generationID]
      
      name_of_file = paste(currentGeneration,layer_ID,"SNP",sep = "_")
      #print(name_of_file)
      
      if(marker_type=="MT"){marker="mtDNA"}
      if(marker_type=="Y"){marker="Y_chrom"}
      if(marker_type=="A"){marker="codom"}
      
      raw_data_per_sim[[generationID]] = read.PLINK(file = paste(name_of_file,"_seqSNP_",marker_type,".raw",sep = ""),
                                                    map.file = paste(name_of_file,"_seqSNP_",marker_type,".map",sep = ""),
                                                    saveNbAlleles = T,
                                                    quiet = T)
      
      # this block processes the xy coordinates from each individuals name
      # since they do not come explicitly in the data file
      xy_coords = NULL
      for(individual in raw_data_per_sim[[generationID]]@ind.names){
        xy_coords= rbind(xy_coords,
                         cbind(unlist(strsplit(individual,split = "_"))[4], #x position
                               unlist(strsplit(individual,split = "_"))[5])) #y position
      }
      xy_coords = apply(xy_coords, 2, as.numeric) #xy coords as numerics instead of characters
      raw_data_per_sim[[generationID]]@other$xy = xy_coords
      
      # adding generation to genlight obj
      raw_data_per_sim[[generationID]]@other$generation = currentGeneration
      
      
    }
    
    raw_snp_data[[simID]] = raw_data_per_sim
    
  }
  
  if(is.character(path_to_genind)==TRUE){
    save(raw_snp_data,file=paste0(path_to_genind,"/",nameOfSimulation,"_",marker,".Rdata"))
  }else{}
  return(raw_snp_data)
  
}


