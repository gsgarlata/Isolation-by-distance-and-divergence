# Run Spatial Simulations in SINS

## 1. Initial Genetic Diversity 
It can be useful to give some initial level of genetic diversity to the population one intends to simulate. This usually reduces the simulation time, since the simulated population reaches mutation-drift equilibrium earlier than if one has to wait for the accumulation of new mutation in the forward-in-time simulator. The script for the `Init_freq` function shown below can be found at the [initial_diversity](initial_diversity) folder. Note that the `Init_freq` function simulates a panmictic population and convert the simulated output into files required for simulations in SINS. Also, note that the working directory must be the one where the genetics files are contained. Only one file is required to be already present in the working directory, that is the "genotype.txt".

  `Init_freq(pop.size, sample.size, locus.type="msat", num.loci, mut.rate, exec, path_to_exec)`

## 2. Simulation in SINS
Once the input files are set up, the simulations can be carried out by typing this command:

`java -jar SINS2.jar -projectName 13x13_K50_m004_hlf_mut4 -outDir output -numberOfSimulations 1 -compress noComp -parallel true -parallelCores 4 -verbose false -outputFormat sins -makeDemographicImages false`

Further details on the meaning of each argument and the possible settings can be found at [SINS](https://github.com/PopConGen/SINS).

## 3. Processing and Analyzing SINS output
After the simulation is concluded, the steps required to reformat SINS output for genetic analyses in `R` are described below. The list of arguments shown here are used across all scripts.

`file: name of the simulated scenario`  
`time_start: initial sampling time (forward-in-time)`  
`time_end: last sampling time (forward-in-time)`  
`int_size: sampling time interval (e.g., how often to sample)`  
`nSTR: number of microsatellite loci`  
`nind: number of sampled individuals`  
`nsims: number of simulation replicates`  
`name: sampling strategy: sampling across different demes within fragment ("random") or sampling in the central deme within fragment ("classic")`  
`scen_patch: file to csv file with information on central demes (for "classic" sampling) or on patches (for "random" sampling).`  
`path_script: the path where the "sampler.sh" is contained`  
`pathToSamplerOut: the path to Sampler output`  
`pathToSinsOut: the path to SINS output` 
`path_info: path to info output`  
`cpu_n: number of CPUs to use for parallelization`  

   - `step1.run_sampler.sh` script in [data_analysis](data_analysis) is used to reformat SINS output in an intermediate text format for analyses.
     
    `./step1.run_sampler.sh $file $nSTR $nsims $path_script $pathToSamplerOut $pathToSinsOut`

   - `step2.run_GetData.sh` script in [data_analysis](data_analysis) is used to extract the genetic data into `adegenet` format in `R` and saved it in an R objejct for each sampled time point.

    `./step2.run_GetData.sh $file $time_start $time_end $nSTR $nsims $int_size $nind $path_script`

   - `step3.run_list_extract.sh` script in [data_analysis](data_analysis) is used to subsample the genetic data extracted in the previous step according to a sampling design specified in the `scen_patch` file (those used in [Sgarlata et al., 2022](https://www.biorxiv.org/content/10.1101/2022.10.26.513874v1) are stored in the [info](../info) folder).

    `./step3.run_list_extract.sh $file $time_start $time_end $int_size $name $scen_patch $path_script $nind $pathToSinsOut $path_info`

   - step4.run_list_FstStats.sh` script in [data_analysis](data_analysis) is used to compute pairwise `Fst` and `Rst` across all pairs of demes. It saves a text file for each sampled time point.

    `./step4.run_list_FstStats.sh $file $time_start $time_end $int_size $name $nSTR $path_script $cpu_n $pathToSinsOut`

   - `step5.run_list_process_Fst.sh` script in [data_analysis](data_analysis) is used to combine in one file the pairwise `Fst` and `Rst` valeus computed for each generation. The separation of this analysis in two steps allows to parallelize the calculation of pairwise Fst for each simulation replicate and generation efficiently.

    `./step5.run_list_process_Fst.sh $file $time_start $time_end $int_size $name $nSTR $pathToSinsOut $path_script`

   - `step6.run_list_IbdStats_parallel.sh` script in [data_analysis](data_analysis) is used to compute isolation-by-distance (IBD) at each generation and simulation replicate by regressing the logarithm of geographical distance with `Fst/(1-Fst)` or `Rst/(1-Rst)`. The results are saved in one text file.

    `./step6.run_list_IbdStats_parallel.sh $file $time_start $time_end $name $nSTR $cpu_n $pathToSinsOut $path_script`

   - `step7.GenDiv_parallel.sh` script in [data_analysis](data_analysis) is used to compute three measures of genetic diversity for each forest fragment, that is expected heterozygosity (`Hexp`), observed heterozygosity (`Hobs`) and number of alleles (`loc.n.all`). This step generates one file for each sampled time point, containig the estimates for all simulation replicates.

    `./step7.GenDiv_parallel.sh $file $time_start $time_end $int_size $name $pathToSinsOut $cpu_n $path_script`

   - `step8.run_list_process_GenDiv.sh` script in [data_analysis](data_analysis) is used to combine in one file the genetic diversity estimated in the previous step.
 
    `./step8.run_list_process_GenDiv.sh $file $time_start $time_end $int_size $name $pathToSinsOut $path_script`

## Plotting Summary Statistics from SINS output

The previous section explained how to analyse SINS output and retrieve summary statistics from the simulated scenario. This section refers to scripts used to plot the summary statistics and how do they change over time, before and after HL&F. All scripts are stored in the [plotting](plotting) folder.

- `plot_IBD_in_time.R` script is used to plot the linear regression of the logarithm of pairwise geographical distance with `Fst/(1 - Fst)` (IBD) over time, as to show the effect of HL&F on IBD.
- `plot_IBDpvalue_in_time.R` script is used to plot the statistical significance of the IBD (measured with the p-value of the Mantel test performed between geographical and genetic distance) over time.
- `plot_pairwiseFst_in_time.R` script is used to plot pairwise Fst over time accounting for geogrpahical distance, that is pairwise Fsts between deme pairs with the sampe geogrpahical distance are averaged together. 
- `plot_withinFragmentGeneticDiversity_in_time.R` script is used to plot within-fragment genetic diversity over time, averaged across all habitat fragments.
  
