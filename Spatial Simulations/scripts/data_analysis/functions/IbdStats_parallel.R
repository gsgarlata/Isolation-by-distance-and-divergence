#!/bin/sh
#SBATCH -J getIbd


file<-args[1]
time_init<-as.numeric(args[2])
time_end<-as.numeric(args[3])
name<-args[4]
nloci = as.numeric(args[5])
cpu_n = as.numeric(args[6])
path_to_data = args[7]
path_functions = args[8]

source(paste0(path_functions,'/ibd_parallel.R'))

ibd.intime_parallel(file,path_to_patgen=path_to_data,time_init,time_end,name=name,nloci,cpu_n)

