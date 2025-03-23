#!/bin/sh
#SBATCH -p workq

file=$1 #name of the simulated scenario
time_init=$2 #initial sampling time (forward-in-time)
time_end=$3 #last sampling time (forward-in-time)
int_size=$4 #sampling time interval (e.g., how often to sample)
name=$5 #sampling strategy: sampling across different demes within fragment ("random") or sampling in the central deme within fragment ("classic")
n_loci=$6 #number of microsatellite loci
pathSinsOut=$9 #path to SINS output
path_script=$8 #the path where the "process_FstStats.sh" is contained.



sbatch $path_script/process_FstStats.sh $file $time_init $time_end $int_size $name $n_loci $pathSinsOut
