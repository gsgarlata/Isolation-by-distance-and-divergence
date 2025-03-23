#!/bin/sh
#SBATCH -J run_Arp2snp
#SBATCH -o fst_parallel.out

file=$1
time_int=$2
time_end=$3
int_size=$4
name=$5
n_loci=$6
path_functions=$7
cpu_n=$8
path_to_data=$9

for generation in $(seq $time_int $int_size $time_end); do

sbatch $path_functions/FstStats_parallel.sh $file $generation $name $n_loci $cpu_n $path_functions $path_to_data

done
