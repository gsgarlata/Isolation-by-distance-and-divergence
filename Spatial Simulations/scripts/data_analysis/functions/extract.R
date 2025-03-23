#!/bin/sh
#SBATCH -J ExtractPops
args = commandArgs(trailingOnly=TRUE)

file = args[1]
generation = as.numeric(args[2])
name = args[3]
patch.file =  args[4]
nInd = as.numeric(args[5])
pathSinsOut = args[6]
path_scripts = args[7]
path_info = args[8]


source(paste0(path_scripts,"/Extract_pops.R"))



name = 'random' # 'classic' or 'random'


ExtractPops(files,
            file.path=pathSinsOut,
            generation=generation,
            path_info,
            patch.file,
            name=name,
            n_samples=nInd)
