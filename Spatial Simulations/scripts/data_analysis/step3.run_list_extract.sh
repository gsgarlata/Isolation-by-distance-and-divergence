#!/bin/sh
#SBATCH -p workq

file=$1 #name of the simulated scenario
time_start=$2 #initial sampling time (forward-in-time)
time_end=$3 #last sampling time (forward-in-time)
int_size=$4 #sampling time interval (e.g., how often to sample)
name=$5 #sampling strategy: sampling across different demes within fragment ("random") or sampling in the central deme within fragment ("classic")
scen_patch=$6 #file to csv file with information on central demes (for "classic" sampling) or on patches (for "random" sampling).
path_script=$7 #the path where the "extract_generations.sh" is contained.
nInd=$8 #number of sampled individuals
pathSinsOut=$9 #path to SINS output
path_info=${10} #path to info output


sbatch $path_script/extract_generations.sh $file $time_start $time_end $int_size $name $scen_patch $path_script $nInd $pathSinsOut $path_info
