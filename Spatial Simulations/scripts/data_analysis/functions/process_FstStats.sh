#!/bin/sh
#SBATCH -p workq

file=$1
time_init=$2
time_end=$3
int_size=$4
name=$5
n_loci=$6
path_input=$7

if [[ "$n_loci" = 0 ]]; then

name=${name}

else

name=${name}_${n_loci}STR

fi

head -n 1 $path_input/$file/${file}_FstStats_codom_${time_init}_raw_${name}.csv > $path_input/$file/${file}_FstStats_codom_${time_init}_${time_end}_${name}.csv

for generation in $(seq $time_init $int_size $time_end); do

tail -n +2 $path_input/$file/${file}_FstStats_codom_${generation}_raw_${name}.csv >> $path_input/$file/${file}_FstStats_codom_${time_init}_${time_end}_${name}.csv

rm $path_input/$file/${file}_FstStats_codom_${generation}_raw_${name}.csv

done
