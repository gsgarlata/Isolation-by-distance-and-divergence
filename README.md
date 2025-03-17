# Isolation-by-distance-and-divergence

This repository provides the scripts to model habitat loss and fragmentation in a two-dimensional torus following the work performed in [Sgarlata et al., 2022](https://www.biorxiv.org/content/10.1101/2022.10.26.513874v1) study entitled: "The effect of habitat loss and fragmentation on isolation-by-distance-and-divergence".

* [mean TMRCA in Torus Underoing Habitat Contraction](HabitatContraction): It models non-equilibrium mean TMRCA in a torus undergoing habitat contraction. The folder includes:

   * an R function to compute the non-equilibrium mean TMRCA in a torus undergoing habitat contraction ([AverageCoalTime_PopSizeChange_Torus.R](HabitatContraction/AverageCoalTime_PopSizeChange_Torus.R));
   * C++ functions (implemented in the R function using the Rcpp package) to estimate the average amount of time required for two lineages sampldee i and j steps apart, in the x- and y-direction respectively, to be in the same deme ([Cpp_functions_Toroidal_barrier.cpp](HabitatContraction/Cpp_functions_Toroidal_barrier.cpp));
   * an R function used to compute the relative spatial coordinates of the central deme within a fragment of a specific size. For instance, for a habitat fragment of size 13 x 13, the central deme has spatial coordinates equal to x = 6 and y = 6. ([GetFragmCoordAllPossibleSize.R](HabitatContraction/GetFragmCoordAllPossibleSize.R));
   * an R script for plotting the mean TMRCA of two alleles sampled in the same deme in a torus undergoing contraction, by considering different dispersal rate and habitat fragment size ([plotMeanTMRCA_TorusHabitatContraction.R](HabitatContraction/plotMeanTMRCA_TorusHabitatContraction.R));

* [mean TMRCA in Torus Underoing Habitat Fragmentation Per Se](HabitatFragmentationPerSe): It models non-equilibrium mean TMRCA for two alleles sampled in different habitat fragments. The folder includes:

   * an R function to compute the weighted Sij for lineages sampled in two habitat fragments (second term in Eq. 4). In particular, it considers all possible locations of the two ancestral lineages (within each habitat fragment) at the time of HL&F and computes the corresponding Sij (Eq. S2). Each Sij value is then weighted by the probability that the ancestral lineages of the sampled alleles were at a given location (within each habitat fragment) at the time of HL&F ([weightedSij_torus_HLF_intime.R](HabitatFragmentationPerSe/weightedSij_torus_HLF_intime.R));
   * an R function that parallelizes the computation of the weighted Sij ([weightedSij_torus_HLF_intime_parallel.R](HabitatFragmentationPerSe/weightedSij_torus_HLF_intime_parallel.R));
   * C++ functions (implemented in the R function using the Rcpp package) to calculate the power "t" of a matrix (temporal convolution), which is used to compute the probability of dispersal at a given location in time "t" ([matrix_power.cpp](HabitatFragmentationPerSe/matrix_power.cpp));
   * an R function used to compute Sij (Eq. S2), that is the average amount of time required for the ancestral lineages of two alleles sampled i and j steps apart to be in the same deme ([Sij_torus_beforeHL&F.R](HabitatContraction/Sij_torus_beforeHL&F.R));
   * an R script for plotting the mean TMRCA of two alleles sampled in two habitat fragments, by considering different dispersal rate and habitat fragment size ([plotMeanTMRCA_TorusHabitatFragmentationPerSe.R](HabitatFragmentationPerSe/plotMeanTMRCA_TorusHabitatFragmentationPerSe.R));

* [pairwise Fst in a Torus Underoing Habitat Loss and Fragmentation](IsolationByDistanceAndDivergence): It models pairwise Fst in a toroidal stepping-stone model undergoing Habitat Loss and Fragmentation. The folder includes:

  * an R function to compute the weighted Sij for lineages sa
