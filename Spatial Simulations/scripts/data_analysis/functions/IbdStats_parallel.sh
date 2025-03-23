#!/bin/sh
#SBATCH -J run_ibd_parallel
#SBATCH -o ibd_parallel.out
#SBATCH -c 10

file=$1
time_init=$2
time_end=$3
name=$4
nloci=$5
cpu_n=$6
path_to_data=$7
path_script=$8

Rscript $scripts_path/IbdStats_parallel.R $file $time_init $time_end $name $nloci $cpu_n $path_to_data $path_script
