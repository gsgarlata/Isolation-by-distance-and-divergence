#!/bin/sh
#SBATCH -J div_parallel
#SBATCH -o div_parallel.out

file=$1 #name of the simulated scenario
generation=$2 #sampling time point
name=$3 #sampling strategy: sampling across different demes within fragment ("random") or sampling in the central deme within fragment ("classic")
pathSinsOut=$4 #path to SINS output
cpu_n=$5 #number of cpu to be used for parallelization
path_script=$6 #the path where the "process_FstStats.sh" is contained.

Rscript --vanilla $path_script/GenDiv_parallel.R $file $generation $name $pathSinsOut $cpu_n $path_script
