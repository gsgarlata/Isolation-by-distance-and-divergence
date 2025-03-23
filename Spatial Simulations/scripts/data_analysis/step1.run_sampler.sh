#!/bin/sh
#SBATCH -p workq

file=$1 #name of the simulated scenario
nstr=$2 #number of microsatellite loci
nsims=$3 #number of simulation replicates
path_script=$4	#the path where the "sampler.sh" is contained.
pathToSamplerOut=$5 #the path to Sampler output
pathToSinsOut=$6 #the path to SINS output

for simID in $(seq 1 $nsims); do

$path_script/sampler.sh $file $nstr $simID $pathToSamplerOut $pathToSinsOut

done
