#!/bin/sh
#SBATCH -p workq


file=$1
time_init=$2
time_end=$3
int_size=$4
name=$5
n_loci=$6
path_script=$7
path_functions=$8
cpu_n=$9
path_to_data=${10}


sbatch $path_script/run_FstStats_parallel.sh $file $time_init $time_end $int_size $name $n_loci $path_functions $cpu_n $path_to_data
