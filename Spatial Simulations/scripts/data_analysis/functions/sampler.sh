#!/bin/sh
#SBATCH -p workq

simulation_name=$1
nMarkers=$2
simID=$3
pathToSamplerOut=$4
pathToSinsOut=$5

directory_path=$pathToSinsOut/$simulation_name

simRepPath=$pathToSamplerOut/${simulation_name}/${simulation_name}/Adegenet_sim_${simID}

mkdir -p $simRepPath

cd $simRepPath

for j in $(seq 1 $nMarkers); do

mypath=$directory_path/simulation_${simID}/layer0_A${j}.txt

gen_sim=_layer0_A${j}

awk -v myvar=_layer0_A${j} '{print >> substr($6,1,length($6)-1)myvar; close(substr($6,1,length($6)-1)myvar)}' $mypath

done
