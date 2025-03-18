# Run Spatial Simulations in SINS


#The working directory must be the one where the genetics files are contained. Only one file is required to be already present
#in the working directory, that is the "genotype.txt".
#This function will generate files for each autosomal microsatellite locus (inlcudoing all the necessary information)

1) It can be useful to give some initial level of genetic diversity to the population one intend to simulate. This usually reduces the simulation time, since the simulated population reaches mutatiojn-drift equilibrium earlier than if one is supposed to generate initial diversity from the forward-in-time simulator.

`Init_freq(pop.size, sample.size, locus.type="msat", num.loci, mut.rate, exec, path_to_exec)`

2) Once the input files are set up, the simulations can be carried out by typing this command:

`java -jar SINS2.jar -projectName 13x13_K50_m004_hlf_mut4 -outDir output -numberOfSimulations 1 -compress noComp -parallel true -parallelCores 4 -verbose false -outputFormat sins -makeDemographicImages false`

3) 
