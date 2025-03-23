#!/bin/sh
#SBATCH -J getFst

args <- commandArgs(TRUE)

file<-args[1]
generation<-as.numeric(args[2])
name<-args[3]
n_loci<-as.numeric(args[4])
cpu_n<-as.numeric(args[5])
path_functions = args[6]
path_to_data = args[7]

source(paste0(path_functions,'/FstStats_parallel_HLF.R'))
source(paste0(path_functions,'/ApplyFunctionToArray.R'))



FstStats_parallel_HLF(file,path_to_data,time_sample=generation,name=name,cpu_n,n.loci=n_loci,path_functions)


