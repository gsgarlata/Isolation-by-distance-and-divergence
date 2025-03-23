#!/bin/sh
#SBATCH -p workq

file=$1 #name of the simulated scenario
nSTR=$2 #number of microsatellite loci
nsims=$3 #number of simulation replicates
nind=$4 #number of sampled individuals
gen_init=$5 #sampling time point
gen_end=$6 #sampling time point
gen_int=$7 ##sampling time interval
path_script=$8 #the path where the "Get_Data.R" is contained.
path_SamplerOut=$9 #the path to Sampler output
path_SinsOut=${10} #the path to SINS output

Rscript --vanilla $path_script/Get_Data.R $file $gen_init $gen_end $gen_int $nsims $nSTR $nind $path_script $path_SamplerOut $path_SinsOut
