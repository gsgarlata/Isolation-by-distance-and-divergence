# Isolation-by-distance-and-divergence

This repository provides the scripts to model habitat loss and fragmentation in a two-dimensional torus following the work performed in [Sgarlata et al., 2022](https://www.biorxiv.org/content/10.1101/2022.10.26.513874v1) study entitled: "The effect of habitat loss and fragmentation on isolation-by-distance-and-divergence".

* [mean TMRCA in Torus Underoing Habitat Contraction](HabitatContraction): It models non-equilibrium mean TMRCA in a torus undergoing habitat contraction. The folder includes:
   * an R function to compute the non-equilibrium mean TMRCA in a torus undergoing habitat contraction ([AverageCoalTime_PopSizeChange_Torus.R](HabitatContraction/AverageCoalTime_PopSizeChange_Torus.R));
   * C++ functions (implemented in the R function using the Rcpp package) to estimate the average amount of time required for two lineages sampldee i and j steps apart, in the x- and y-direction respectively, to be in the same deme.
