##it is critical that there are no spaces in your gene names for this analysis. This code relies on 
  #comparison of gene names as strings, therefore if there are extra spaces or other characters 
  #the program will evaluate the same gene as different strings.
#it is also critical that datasets be imported as tsv or csv tables in this format: 
#   clone1name  clone2name  clone3name ...
#   genex       geney       genez
#   geney       genex       genez
# stringsAsFactors=FALSE!! 

Generate_Pairs_Dataframe<-function(Genes_by_clone_dataframe){ 
  Genes_by_clone_dataframe[Genes_by_clone_dataframe==""]<-NA #pad all empty cells with NA
  all.genes<-list() #loop will fill all.genes with all genes in dataset
    for (i in names(Genes_by_clone_dataframe)){
      list<-na.omit(Genes_by_clone_dataframe[, i], stringsAsFactors=FALSE) 
      all.genes<-append(all.genes, list)
      all.genes<-all.genes[!is.na(all.genes)] #remove NA
      all.genes<-all.genes[!duplicated(all.genes)] #remove duplicates 
      }
  all_possible_pairs<-combn(all.genes, 2, simplify=FALSE, stringsAsFactors=FALSE)
      k=length(all_possible_pairs)
      All.Pairs.df <-data.frame(matrix(unlist(all_possible_pairs),nrow=k, byrow=T), stringsAsFactors=FALSE)
      sigma1.sigma2<- paste(pmin(All.Pairs.df$X1, All.Pairs.df$X2), pmax(All.Pairs.df$X1, All.Pairs.df$X2))
      All.Pairs<-cbind(sigma1.sigma2, All.Pairs.df, stringsAsFactors=FALSE)
      colnames(All.Pairs)[c(2, 3)]<-c("sigma1", "sigma2")
    return(All.Pairs)
      }