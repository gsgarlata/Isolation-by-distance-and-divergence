CircularIndexingMatrix2Dataframe<-function(dxToroid,dyToroid){

toroidal<-matrix(0,nrow=dxToroid,ncol=dyToroid) #matrix with row X col corresponding to all possible X and Y pairwise combinations
colnames(toroidal)<-0:(ncol(toroidal)-1) #give column names corresponding to the numbering in a 1d-circular model [0:(d-1)]
rownames(toroidal)<-0:(nrow(toroidal)-1) #give row names corresponding to the numbering in a 1d-circular model [0:(d-1)]

Tmat<-toroidal #make a copy of the toroidal matrix for further editing
first_row_1x<-((dxToroid/2) + 2):nrow(Tmat) #get row indexes for the 'reflective' ('moving backward') demes in a 1d circular model.
sec_row_1x<-((dxToroid/2)-1):0 #get row indexes for the 'moving forward' demes in a 1d circular model.

#the same as for the x-direction but for the y-direction.
first_row_1y<-((dyToroid/2) + 2):ncol(Tmat)
sec_row_1y<-((dyToroid/2)-1):0

#re-number the row and columns of the matrix considering that a 2D torus includes back
#and forward movement due to 'reflective edges'.  
rownames(Tmat)[first_row_1x]<-sec_row_1x[1:length(first_row_1x)]
colnames(Tmat)[first_row_1y]<-sec_row_1y[1:length(first_row_1y)]

df_xycoord<-melt(toroidal)
names(df_xycoord)<-c('x','y','z')

df_Tmat<-melt(Tmat)
names(df_Tmat)[c(1,2)]<-c('i','j')

df_Tmat$x<-df_xycoord$x
df_Tmat$y<-df_xycoord$y

return(df_Tmat)
}
