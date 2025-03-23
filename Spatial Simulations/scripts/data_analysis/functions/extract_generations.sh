#!/bin/sh
#SBATCH -o extract.log
#SBATCH -p workq

file=$1
time_start=$2
time_end=$3
int_size=$4
name=$5
scen_patch=$6
path_script=$7
nInd=$8
pathSinsOut=$9
path_info=${10}

for generation in $(seq $time_start $int_size $time_end); do

sbatch $path_script/extract.sh $file $generation $name $nInd $scen_patch $path_script $pathSinsOut $path_info

done
