#!/bin/sh
#SBATCH -J GetData

args = commandArgs(trailingOnly=TRUE)

nameOfSimulation <-args[1]
gen_init<-as.numeric(args[2])
gen_end<-as.numeric(args[3])
gen_int<-as.numeric(args[4])
nsimID = as.numeric(args[5])
n_loci = as.numeric(args[6])
nSample = as.numeric(args[7])
path_script = args[8]
path_SamplerOut = args[9]
path_SinsOut = args[10]


source(paste0(path_script,"/SINS_get_data_for_analysis.R"))

currentDirectory<-getwd()

setwd(currentDirectory)


the_number_of_simulations <- Set_number_of_simulations(nsimID)
the_generation_list <-seq(gen_init,gen_end,gen_int)
the_layers_list <- Set_layers_list("layer0")
the_markers_list <- Set_markers_list(paste0("A",1:n_loci))
ind_to_sample<-as.numeric(nSample)

marker_type<-"codom"


for(t in 1:length(the_generation_list)){

theFileNames <- Get_file_names(the_generation_list[t], the_layers_list)
thePathToData <- paste(path_SamplerOut,"/",nameOfSimulation,"/",nameOfSimulation,"/",sep="")
path_to_genind<-paste0(path_SinsOut,"/",nameOfSimulation,"/")
genindNestedArray <- Get_data_from_files_multi_sim(the_file_names = theFileNames,
                                                   set_markers = the_markers_list,
                                                   generations_list = the_generation_list[t],
                                                   path_to_data = thePathToData,
                                                   number_of_simulations = the_number_of_simulations,
                                                   is_marker_diploid = TRUE,
                                                   marker_type,
                                                   path_to_genind,
                                                   subind=ind_to_sample)

}

