#!/bin/sh
#SBATCH -J run_fst_parallel
#SBATCH -o fst_parallel.out
#SBATCH -c 10


file=$1
generation=$2
name=$3
n_loci=$4
cpu_n=$5
path_functions=$6
path_to_data=$7

Rscript $path_functions/FstStats.R $file $generation $name $n_loci $cpu_n $path_functions $path_to_data
