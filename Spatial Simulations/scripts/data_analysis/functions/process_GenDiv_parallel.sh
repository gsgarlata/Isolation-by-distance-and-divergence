#!/bin/sh
#SBATCH -p workq

file=$1
time_init=$2
time_end=$3
int_size=$4
name=$5
path_input=$6

head -n 1 $path_input/$file/${file}_diversity_raw_${time_init}_${name}.csv > $path_input/$file/${file}_diversity_raw_${time_init}_${time_end}_${name}.csv

for generation in $(seq $time_init $int_size $time_end); do

tail -n +2 $path_input/$file/${file}_diversity_raw_${generation}_${name}.csv >> $path_input/$file/${file}_diversity_raw_${time_init}_${time_end}_${name}.csv

rm $path_input/$file/${file}_diversity_raw_${generation}_${name}.csv

done
