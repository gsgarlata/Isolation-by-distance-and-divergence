#!/bin/sh
#SBATCH -p workq

file=$1 #name of the simulated scenario
time_start=$2 #initial sampling time (forward-in-time)
time_end=$3 #last sampling time (forward-in-time)
nSTR=$4 #number of microsatellite loci
nsims=$5 #number of simulation replicates
int_size=$6 #sampling time interval (e.g., how often to sample)
nind=$7 #number of sampled individuals
path_script=$8 #the path where the "Get_Data.sh" is contained.


time_tot=$(seq $time_start $int_size $time_end)

for time_point in $time_tot; do

$path_script/Get_Data.sh $file $nSTR $nsims $nind $time_point $time_point $int_size $path_script

done
