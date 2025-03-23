#!/bin/sh
#SBATCH -p workq

file=$1 #name of the simulated scenario
time_init=$2 #initial sampling time (forward-in-time)
time_end=$3 #last sampling time (forward-in-time)
name=$4 #sampling strategy: sampling across different demes within fragment ("random") or sampling in the central deme within fragment ("classic")
n_loci=$5 #number of microsatellite loci
cpu_n=$6 #number of cpu to be used for parallelization
pathSinsOut=$7 #path to SINS output
path_script=$8 #the path where the "process_FstStats.sh" is contained.

sbatch $path_script/IbdStats_parallel.sh $file $time_init $time_end $name $nloci $cpu_n $pathSinsOut $path_script
