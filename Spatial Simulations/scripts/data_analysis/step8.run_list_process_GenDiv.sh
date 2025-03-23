#!/bin/sh
#SBATCH -p workq

file=$1
time_init=$2
time_end=$3
int_size=$4
name=$5
path_input=$6
path_script=$7

sbatch $path_script/process_GenDiv_parallel.sh $file $time_init $time_end $int_size $name $path_input
