# Run Spatial Simulations in SINS

1) It can be useful to give some initial level of genetic diversity to the population one intends to simulate. This usually reduces the simulation time, since the simulated population reaches mutation-drift equilibrium earlier than if one has to wait for the accumulation of new mutation in the forward-in-time simulator. The script for the `Init_freq` function shown below can be found at the [initial_diversity](initial_diversity) folder. Note that the `Init_freq` function simulates a panmictic population and convert the simulated output into files required for simulations in SINS. Also, note that the working directory must be the one where the genetics files are contained. Only one file is required to be already present in the working directory, that is the "genotype.txt".

`Init_freq(pop.size, sample.size, locus.type="msat", num.loci, mut.rate, exec, path_to_exec)`

2) Once the input files are set up, the simulations can be carried out by typing this command:

`java -jar SINS2.jar -projectName 13x13_K50_m004_hlf_mut4 -outDir output -numberOfSimulations 1 -compress noComp -parallel true -parallelCores 4 -verbose false -outputFormat sins -makeDemographicImages false`

Further details on the meaning of each argument and the possible settings can be found at [SINS](https://github.com/PopConGen/SINS).

3) After the simulation is concluded, one way to analyse SINS output is through the `step1.run_sampler.sh` script in [data_analysis](data_analysis).

`file: name of the simulated scenario`  
`nstr: number of microsatellite loci`  
`nsims: number of simulation replicates`  
`path_script: the path where the "sampler.sh" is contained`  
`pathToSamplerOut: the path to Sampler output`  
`pathToSinsOut: the path to SINS output`  

`./step1.run_sampler.sh $file $nstr $nsims $path_script $pathToSamplerOut $pathToSinsOut`

4) `step2.run_GetData.sh` script in [data_analysis](data_analysis).

`file: name of the simulated scenario`  
`time_start: initial sampling time (forward-in-time)`  
`time_end: last sampling time (forward-in-time)`  
`nSTR: number of microsatellite loci`  
`nsims: number of simulation replicates`  
`int_size: sampling time interval (e.g., how often to sample)`  
`nind: number of sampled individuals`  
`path_script: the path where the "Get_Data.sh" is contained.`  

`./step2.run_GetData.sh $file $time_start $time_end $nSTR $nsims $int_size $nind $path_script`

5) `step3.run_list_extract.sh` script in [data_analysis](data_analysis), `scen_patch` files are stored in the [info](../info) folder.

`file: name of the simulated scenario`  
`time_start: initial sampling time (forward-in-time)`  
`time_end: last sampling time (forward-in-time)`  
`int_size: sampling time interval (e.g., how often to sample)`  
`name: sampling strategy: sampling across different demes within fragment ("random") or sampling in the central deme within fragment ("classic")`  
`scen_patch: file to csv file with information on central demes (for "classic" sampling) or on patches (for "random" sampling).`  
`path_script: the path where the "Get_Data.sh" is contained.`  
`nind: number of sampled individuals`  
`pathSinsOut: path to SINS output`   
`path_info: path to info output`  

`./step3.run_list_extract.sh $file $time_start $time_end $int_size $name $scen_patch $path_script $nind $pathSinsOut $path_info`

6) step4.run_list_FstStats.sh` script in [data_analysis](data_analysis).

file=$1
time_start=$2
time_end=$3
int_size=$4
name=$5
nSTR=$6
path_script=$7
path_functions=$8
cpu_n=$9
path_to_data=${10}

./step4.run_list_FstStats.sh $file $time_start $time_end $int_size $name $nSTR 


7) 
