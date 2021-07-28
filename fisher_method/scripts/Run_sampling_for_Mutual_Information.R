##Main function to sample dataset for MI and to simulate MI. If no simulations are desired, set replicates to 0. 

#Genes_by_clone_dataframe = raw data in format below
#Ex. 
#   clone1name  clone2name  clone3name ...
#   genex       geney       genez
#   geney       genex       genex
# stringsAsFactors=FALSE!! 

#When importing datasets stringsAsFactors=FALSE!! 

#Requires sourcing Construct_Pairs.R
#Requires sourcing Generate_Pairs_Dataframe.R
#Requires sourcing Probability_Sample.R
#Requires sourcing Build_Probability_Matrix.R

#Required packages: 
  #plyr

Run_Sampling_for_MI<-function(Genes_by_clone_dataframe, replicates){ 
    
  library(plyr)
  
  Construct_Pairs(Genes_by_clone_dataframe) -> All.pairs
  if (replicates>0){
  for (i in 1:replicates){
    permutedpairscount<-Probability_Sample(Genes_by_clone_dataframe)
    permutedpairscount<- permutedpairscount[permutedpairscount$sigma1.sigma2!=0, ]
    P1.1<- permutedpairscount[, 5]
    P0.1<- permutedpairscount[, 6]
    P1.0<- permutedpairscount[, 7]
    P0.0<- permutedpairscount[, 8]
    Psig1<- permutedpairscount[, 9]
    Psig2<- permutedpairscount[, 10]
    I1.1<- P1.1 * log2(P1.1/(Psig1*Psig2))
    I0.1<- P0.1 * log2(P0.1/((1-Psig1)*Psig2))
    I1.0<- P1.0 * log2(P1.0/(Psig1 *(1-Psig2)))
    I0.0<- P0.0 * log2(P0.0/((1-Psig1)*(1-Psig2)))
    I.df<- data.frame(permutedpairscount[1], I1.1, I0.1, I1.0, I0.0)
    I.df$I<- (rowSums(I.df[2:5]))
    I.df<-I.df[c(1, 6)]
    colnames(I.df)[2]<-(sprintf("MI sample %d", i))
    All.pairs<-merge(All.pairs, I.df, by="sigma1.sigma2", all=TRUE)
  }}
  colnames(All.pairs)[1]<-"Gene Pair"
  return(All.pairs)
}