  
#function to generate a matrix of marginal probabilities and Bernoulli r.v. success 
# probabilities from a dataframe wherein colummns are genotypes. 

Build_Probability_Matrix<-function(Genes_by_clone_dataframe){
  library("pracma")
  all.genes<-list() #loop will fill all.genes with all genes in dataset
  for (h in names(Genes_by_clone_dataframe)){
    list<-na.omit(Genes_by_clone_dataframe[, h], stringsAsFactors=FALSE) 
    all.genes<-append(all.genes, list)
  } 
  df<-data.frame(matrix(unlist(all.genes),nrow=length(all.genes), byrow=T), stringsAsFactors=FALSE)
  colnames(df)[1]<- "sigma"
  Probability_Matrix<-count(df, "sigma")
  N<- ncol(Genes_by_clone_dataframe)
  Probability_Matrix$p.sigma<- Probability_Matrix$freq/N
  #defining all the variables needed to find success probability
  #finding the success probability of Bernoulli random variable per gene
  return(Probability_Matrix)
}