# Isolation-by-distance-and-divergence

This repository provides the scripts to model habitat loss and fragmentation in a two-dimensional torus following the work performed in [Sgarlata et al., 2022](https://www.biorxiv.org/content/10.1101/2022.10.26.513874v1) study entitled: "The effect of habitat loss and fragmentation on isolation-by-distance-and-divergence". In particular, it includes: the scripts for the numerical model of "isolation-by-distance-and-divergence" 

* [Isolation-by-distance-and-divergence model](Numerical Model): It contains scripts for modelling pairwise Fst in a toroidal stepping-stone model undergoing habitat loss and fragmentation. The model is an extenstion of the original model of isolation-by-distance in a toroidal stepping-stone model by [Slatkin, 1991](https://www.cambridge.org/core/journals/genetics-research/article/inbreeding-coefficients-and-coalescence-times/FCC418CBC6F021B741C83FDE6A0E7558) and [Slatkin, 1993] (https://www.jstor.org/stable/2410134?origin=crossref&seq=1), and it uses an approach like the one in [Duforet-Frebourg and Slatkin, 2016](https://www.sciencedirect.com/science/article/abs/pii/S0040580915001124?via%3Dihub) for modelling the mean time to the most recent common ancestor (TMRCA) and pairwise Fst for two alleles sampled at different time points (e.g., modern vs ancient DNA).

  
* [Spatial Genetic Simulations of HL&F in a two-dimensional plane stepping-stone model](Spatial Simulation): It models pairwise Fst in a toroidal stepping-stone model undergoing Habitat Loss and Fragmentation. The folder includes:
