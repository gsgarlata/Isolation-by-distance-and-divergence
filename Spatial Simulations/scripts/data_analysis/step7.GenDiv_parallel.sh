#!/bin/sh
#SBATCH -J run_Arp2snp
#SBATCH -o div_parallel.out


file=$1
time_int=$2
time_end=$3
int_size=$4
name=$5
pathSinsOut=$6
cpu_n=$7
path_script=$8

sbatch $path_script/run_GenDiv_parallel.sh $file $time_int $time_end $int_size $name $pathSinsOut $cpu_n $path_script

done
