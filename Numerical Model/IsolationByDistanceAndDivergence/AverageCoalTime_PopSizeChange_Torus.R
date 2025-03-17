#This script is used to compute the non-equilibrium mean TMRCA of two alleles sampled in the same deme in a torus
#undergoing Habitat Contraction.

#nd_full: number of demes on the x-axis (equal to the number of demes on the y-axis) before Habitat Contraction
#dx: number of demes on the x-axis after Habitat Contraction
#dy: number of demes on the y-axis after Habitat Contraction
#migr: deme dispersal rate in a 1D circular stepping stone model (i.e., migr = migr2D / 2, where "migr2D" is teh dispersal rate in a 2D toroidal stepping stone model) 
#Ne: diploid deme effective population size
#time_end: time duration of the contraction for which one wishes to compute the mean TMRCA.
#path_func: path to R functions necessary for the calculations


AverageCoalTime_PopSizeChangeT00Torus<-function(nd_full,dx,dy,migr,Ne,time_end, path_func){

require(reshape2)
require(Rcpp)
require(DescTools)
require(ggplot2)

sourceCpp(paste0(path_func,'/Cpp_functions_Toroidal_barrier.cpp'))
source(paste0(path_func,'/CircularIndexingMatrix2Dataframe.R'))


#START: build a circular matrix#
Ndoogan.Circulant <-function(x) {
  n <- length(x)
  suppressWarnings(
    matrix(x[matrix(1:n,n+1,n+1,byrow=T)[c(1,n:2),1:n]],n,n))
}
#END: build a circular matrix#

#"dxToroid" defines the number of demes in a 1d circular stepping stone, approximately equivalent to a linear 1d stepping stone with "dx" demes.
#Maruyama, (1971) showed that L*L plane is equivalent to a 2L * 2L toroid.
dxToroid<-2*(dx-1) 
dyToroid<-2*(dy-1)

#"a00" describes the probability that two lineages sampled in the same deme (i=0) 
#were also in the same deme (j=0) in the previous generation (Eq. S8a).
a00<-((1-migr)^2)+((migr^2)/2)
#"a01" describes the probability that two lineages sampled in the same deme (i=0) 
#were were at distance j=1 in the previous generation (Eq. S8b).
a01<-migr*(1-migr)
#"a02" describes the probability that two lineages sampled in the same deme (i=0) 
#were were at distance j=2 in the previous generation (Eq. S8c).
a02<-(migr^2)/4

#It defines the number of entries in a migration matrix that have zero values (since migration occurs only with nearby demes).
#It "dxToroid-5" because 5 corresponds to the demes for which there non-zero dispersal (focal deme + 4 neighboring demes).
zero_vals<-dxToroid-5

#It fills the migration matrix entries with the terms defined above.
if(zero_vals<0){
  M_vec<-c(a00,a01,2*a02,a01)
}else{
M_vec<-c(a00,a01,a02,rep(0,zero_vals),a02,a01)
}

#It generates the migration matrix in a 1D circular stepping-stone model with size "dxToroid". 
#In this case we assume dxToroid = dyToroid, so it is sufficient to compute the 1D migration matrix only once.
Mmat<-Ndoogan.Circulant(M_vec)

#"T00" defines the average coalescence time of two alleles sampled in the same deme in a torus with size "nd_full" x "nd_full" (before Habitat Contraction).
T00<-2*Ne*(nd_full*nd_full)

#It generates a 2D dataframe with the x-y coordinates of a torus of size "dxToroid" x "dyToroid" (after Habitat Contraction).
df_Tmat<-CircularIndexingMatrix2Dataframe(dxToroid,dyToroid)

#"series_Sij" computes the average amount of time (Sij) for two alleles sampled at distance i (in the x-axis) and j (in the y-axis) (Sij is computed using Eq. S2).
#This two lines of code are used to compute the matrix T in Eq. 5. The column "value" correspond to the mean TMRCA of two lineages initially separated 
#by a distance i and j on the x- and y-direction of the torus before Habitat Contraction.
df_Tmat$Sij<-series_Sij(df_Tmat$i, df_Tmat$j, nd_full, nd_full, migr)
df_Tmat$value<- T00 + df_Tmat$Sij

#It converts the df_Tmat dataframe into a matrix.
TfinalMat<-as.matrix.xtabs(xtabs(value~x+y, data=df_Tmat))
attr(TfinalMat, "dimnames")<-NULL

#It creates the "T0" matrix in Eq. 5.
Tmat00<-matrix(0,nrow = nrow(TfinalMat),ncol = ncol(TfinalMat))
Tmat00[1,1]<-TfinalMat[1,1]/(2*Ne)

#It creates the "U" matrix in Eq. 5.
Umat<-matrix(1,nrow=nrow(TfinalMat),ncol=ncol(TfinalMat))

#It creates the dataframe where we will save the mean TMRCA of two alleles sampled from the same deme at different time points after Habitat Contraction. 
df_t00<-data.frame(time=1:time_end,t00=rep(NA,time_end))

#Loop over each time point after Habitat Contraction
for(i in 1:time_end){

#it computes the mean coalescence time at the next generation, corresponding to Eq. 9. 
Tp = (Mmat %*% (TfinalMat - Tmat00) %*% Mmat) + Umat

TfinalMat = Tp 
Tmat00[1,1] = TfinalMat[1,1]/(2*Ne)

df_t00$t00[i]<-TfinalMat[1,1]

}

return(df_t00)

}
