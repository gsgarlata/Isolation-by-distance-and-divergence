#!/bin/sh
#SBATCH -J extract_pops
#SBATCH -o extract_pops

file=$1
generation=$2
name=$3
nInd=$4
scen_patch=$5
path_script=$6
pathSinsOut=$7
path_info=$8

Rscript --vanilla $path_script/extract.R $file $generation $name $scen_patch $nInd $pathSinsOut $path_script $path_info
