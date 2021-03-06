library(plyr)
library(pracma)
#Read In data 


library(plyr)
library(pracma) 


example <- read.delim("Example_genotype_data.tsv",  
                      header=TRUE, stringsAsFactors=FALSE)
  
test <- read.delim("dauer_gene3.tsv", 
                    header=TRUE, stringsAsFactors=FALSE)


####Parallel GI Itot####
######################

##Functions## 

#Run Sampling for Mutual Information 
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


  Construct_Pairs = function(Genes_by_clone_dataframe) {
  #library(plyr)
  
  #first set up an empty dataframe
  
  all.data = data.frame(matrix(vector(), 0, 2,
                               dimnames=list(c(), c("X1", "X2"))),
                        stringsAsFactors=F)
  #count all how many times each gene appears
  all.genes<-list() #loop will fill all.genes with all genes in dataset
  for (h in names(Genes_by_clone_dataframe)){
    list<-na.omit(Genes_by_clone_dataframe[, h], stringsAsFactors=FALSE) 
    all.genes<-append(all.genes, list)
  } 
  df<-data.frame(matrix(unlist(all.genes),nrow=length(all.genes), byrow=T), stringsAsFactors=FALSE)
  colnames(df)[1]<- "sigma2"
  gene_counts.2<-count(df, "sigma2")
  M<- sum(gene_counts.2[, 2])
  epsilon<-1/M
  N<- ncol(Genes_by_clone_dataframe)
  rm(df, all.genes)
  Genes_by_clone_dataframe-> x
  
  for(i in names(Genes_by_clone_dataframe)){
    #combn() generates all pairwise combinations
    list<-na.omit(Genes_by_clone_dataframe[, i])
    generatepairs<-combn(list, 2, simplify=FALSE, stringsAsFactors=FALSE)
    # convert the string generated by combn() to a dataframe 
    #with the 2 columns being the two genes in the pair
    pairs.to.df <-data.frame(matrix(unlist(generatepairs),nrow=length(generatepairs), byrow=T), stringsAsFactors=FALSE)
    all.data<-rbind(all.data, pairs.to.df, stringsAsFactors=FALSE)
  }
  
  All.Pairs<-Generate_Pairs_Dataframe(Genes_by_clone_dataframe)
  ##sig1countdf will store sigma1 frequencies for later  
  sig1countdf<- All.Pairs[1:2]
  gene_counts.1<-gene_counts.2
  colnames(gene_counts.1)[1]<- "sigma1"
  sig1countdf<- merge(sig1countdf, gene_counts.1, by="sigma1")
  colnames(sig1countdf)[3]<-"sigma1.count"
  sig1countdf<-sig1countdf[2:3]
  
  ## take the above dataframe and add a column that 
  #contains both gene names using paste()
  # then use cbind() to add the new column to the dataframe
  sigma1.sigma2<- paste(pmin(all.data$X1, all.data$X2), pmax(all.data$X1, all.data$X2))
  ## old way without alphebetizing
  all.data<-cbind(sigma1.sigma2, all.data)
  pairscount<-count(all.data, "sigma1.sigma2") 
  # merge by pairs
  All.Pairs<-merge(All.Pairs, pairscount, by="sigma1.sigma2", all=TRUE) 
  # ## rename the column
  colnames(All.Pairs)[4]<-"Observed_sigma1.sigma2_Frequncy" 
  All.Pairs<-merge(All.Pairs, gene_counts.2, by="sigma2", all=TRUE) #col 5 is sigma 2 freq
  colnames(All.Pairs)[5]<-"sigma2_Observed_Frequency"
  
  All.Pairs<-All.Pairs[, c(1, 5, 2, 4)]
  All.Pairs[is.na(All.Pairs)]<-0
  #replace na values for non-occuring pairs with 0
  All.Pairs<- All.Pairs[All.Pairs$sigma1.sigma2!= "FALSE FALSE", ]
  All.Pairs<-merge(All.Pairs, sig1countdf, by="sigma1.sigma2")
  
  All.Pairs<-All.Pairs[, c(2,3,1,4,5)]  
  C0<- 1/(N*(1+epsilon))
  C <- 1/(N*(1+epsilon)^2)
  Msig1<- (All.Pairs[,5] + (N*epsilon))
  Msig2<- (All.Pairs[,2] + (N*epsilon))
  Po1<- Msig1 * C0
  Po2<- Msig2 * C0 
  Clones0.0<- N - ((All.Pairs[, 2]+ All.Pairs[, 5])-All.Pairs[, 4])
  
  All.Pairs$"Po1,o2(1,1)"<- C* ((((1+epsilon)^2)*All.Pairs[,4]) + 
                                  (((1+epsilon)*epsilon)* ((All.Pairs[,2]-All.Pairs[,4]) + (All.Pairs[,5]-All.Pairs[,4]))) +
                                  ((epsilon^2)*Clones0.0))
  
  All.Pairs$"Po1,o2(0,1)"<- C* ((All.Pairs[,2]-All.Pairs[,4])+((All.Pairs[,2]-All.Pairs[,4])*epsilon)+(Clones0.0*epsilon))
  All.Pairs$"Po1,o2(1,0)"<- C* ((All.Pairs[,5]-All.Pairs[,4])+((All.Pairs[,5]-All.Pairs[,4])*epsilon)+(Clones0.0*epsilon))
  All.Pairs$"Po1,o2(0,0)"<- C * Clones0.0
  All.Pairs$"Po1"<- Po1
  All.Pairs$"Po2"<- Po2
  All.Pairs<- All.Pairs[c(3,4,1,2,6,7,8,9,10,11)]
  P1.1<- All.Pairs[, 5]
  P0.1<- All.Pairs[, 6]
  P1.0<- All.Pairs[, 7]
  P0.0<- All.Pairs[, 8]
  Psig1<- All.Pairs[, 9]
  Psig2<- All.Pairs[, 10]
  I1.1<- P1.1 * log2(P1.1/(Psig1*Psig2))
  I0.1<- P0.1 * log2(P0.1/((1-Psig1)*Psig2))
  I1.0<- P1.0 * log2(P1.0/(Psig1 *(1-Psig2)))
  I0.0<- P0.0 * log2(P0.0/((1-Psig1)*(1-Psig2)))
  I.df<- data.frame(All.Pairs[1], I1.1, I0.1, I1.0, I0.0)
  I.df$I<- (rowSums(I.df[2:5]))
  I.df<-I.df[c(1, 6)]
  colnames(I.df)[2]<- "Observed  MI"
  return(I.df)
} 


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

    #Probability Sample 
    Probability_Sample<-function(Genes_by_clone_dataframe){
      #library(plyr)
      
      #first step is to generate a matrix of marginal probabilities using the build_probability_matrix function
      Probability_Matrix<-Build_Probability_Matrix(Genes_by_clone_dataframe)
      
      #sample vectors of appropriate sizes (based on original data set) and missing data can be filled with NA's. 
      #empty dataframe and replace with NA
      all.genes<-list() #loop will fill all.genes with all genes in dataset
      for (h in names(Genes_by_clone_dataframe)){
        list<-na.omit(Genes_by_clone_dataframe[, h], stringsAsFactors=FALSE) 
        all.genes<-append(all.genes, list)
      } 
      df_to_sample<-data.frame(matrix(unlist(all.genes),nrow=length(all.genes), byrow=T), stringsAsFactors=FALSE)
      df_to_sample<-df_to_sample[!duplicated(df_to_sample), ]
      col_names<- names(Genes_by_clone_dataframe)
      Sampled_dataframe<- data.frame(matrix(NA, nrow=length(df_to_sample), ncol=length(col_names)+1))
      Sampled_dataframe[,1]<- df_to_sample
      colnames(Sampled_dataframe)[2:ncol(Sampled_dataframe)]<- col_names
      colnames(Sampled_dataframe)[1]<-"sigma"
      Sampled_dataframe<-merge(Probability_Matrix, Sampled_dataframe, by="sigma", all=TRUE)
      # up to this point I set up an empty matrix for sampling and a probability vector
      
      #repeat loop to ensure all genotypes get at least 2 genes
      repeat {
        for (i in 1:nrow(Sampled_dataframe)){
          for(j in 4:ncol(Sampled_dataframe)){
            Sampled_dataframe[i,j]<- rbinom(1, 1, Sampled_dataframe[i,3]) 
          }} # binomial sampling to fill dataframe 
        LessThanTwo<-FALSE
        for (p in 4:ncol(Sampled_dataframe)){
          if (sum(Sampled_dataframe[,p])<2){
            LessThanTwo<-TRUE
          }
        }
        if(!LessThanTwo) break
      }
      for (k in 4:ncol(Sampled_dataframe)){
        for (l in 1:nrow(Sampled_dataframe)){
          if (Sampled_dataframe[l,k]==1){
            Sampled_dataframe[l,k]<-Sampled_dataframe[l,1]
          }
        }
      } 
      Sampled_dataframe[Sampled_dataframe==0]<-NA
      Sampled_dataframe<-Sampled_dataframe[4:ncol(Sampled_dataframe)] 
      # Sampled_dataframe complete
      
      # recount gene frequency     
      all.genes<-list() #loop will fill all.genes with all genes in dataset
      for (h in 1:ncol(Sampled_dataframe)){
        list<-na.omit(Sampled_dataframe[, h], stringsAsFactors=FALSE) 
        all.genes<-append(all.genes, list)
      } 
      df<-data.frame(matrix(unlist(all.genes),nrow=length(all.genes), byrow=T), stringsAsFactors=FALSE)
      colnames(df)[1]<- "sigma2"
      gene_counts.2<-count(df, "sigma2")
      M<- sum(gene_counts.2[, 2])
      epsilon<-1/M
      N<- ncol(Genes_by_clone_dataframe)
      rm(df, all.genes)
      #conduct the normal Pairs analysis on the Sampled_dataframe 
      # code below taken directly from Construct_pairs function 
      all.data = data.frame(matrix(vector(), 0, 2,
                                   dimnames=list(c(), c("X1", "X2"))),
                            stringsAsFactors=F)
      
      for(i in 1:ncol(Sampled_dataframe)){
        #combn() generates all pairwise combinations
        list<-na.omit(Sampled_dataframe[, i])
        generatepairs<-combn(list, 2, simplify=FALSE, stringsAsFactors=FALSE)
        k=length(generatepairs)
        # convert the string generated by combn() to a dataframe 
        #with the 2 columns being the two genes in the pair
        pairs.to.df <-data.frame(matrix(unlist(generatepairs),nrow=k, byrow=T), stringsAsFactors=FALSE)
        all.data<-rbind(all.data, pairs.to.df, stringsAsFactors=FALSE)
      }
      ## take the above dataframe and add a column that 
      #contains both gene names using paste()
      # then use cbind() to add the new column to the dataframe
      
      All.Pairs<-Generate_Pairs_Dataframe(Genes_by_clone_dataframe)
      sigma1.sigma2<- paste(pmin(all.data$X1, all.data$X2), pmax(all.data$X1, all.data$X2))
      ## old way without alphebetizing# pairs<-paste(pairs.to.df[,1], pairs.to.df[,2])
      all.data<-cbind(sigma1.sigma2, all.data)
      pairscount<-count(all.data, "sigma1.sigma2")
      All.Pairs<-merge(All.Pairs, pairscount, by="sigma1.sigma2", all=TRUE)
      colnames(All.Pairs)[4]<-"sigma1.sigma2_Frequncy" 
      #columns 2 and 3 are sigma1 and sigma 2
      All.Pairs<-merge(All.Pairs, gene_counts.2, by="sigma2", all=TRUE) 
      sig1countdf<-All.Pairs[2:3]
      gene_counts.2->gene_counts.1
      colnames(gene_counts.1)<-c("sigma1", "sigma1.count")
      sig1countdf<-merge(sig1countdf, gene_counts.1, by="sigma1")
      sig1countdf<-sig1countdf[2:3]
      colnames(All.Pairs)[5]<-"sigma2_Frequency"
      #replace na values for non-occuring pairs with 0
      All.Pairs<-All.Pairs[, c(1, 5, 2, 4)]
      All.Pairs<- All.Pairs[All.Pairs$sigma1.sigma2!= "FALSE FALSE", ]
      All.Pairs<-merge(All.Pairs, sig1countdf, by="sigma1.sigma2", all=TRUE)
      All.Pairs<-All.Pairs[, c(2,3,1,4,5)]
      All.Pairs[is.na(All.Pairs)]<-0
      ######
      C0<- 1/(N*(1+epsilon))
      C <- 1/(N*(1+epsilon)^2)
      Msig1<- (All.Pairs[,5] + (N*epsilon))
      Msig2<- (All.Pairs[,2] + (N*epsilon))
      Po1<- Msig1 * C0
      Po2<- Msig2 * C0 
      Clones0.0<- N - ((All.Pairs[, 2]+ All.Pairs[, 5])-All.Pairs[, 4])
      
      All.Pairs$"Po1,o2(1,1)"<- C* ((((1+epsilon)^2)*All.Pairs[,4]) + 
                                      (((1+epsilon)*epsilon)* ((All.Pairs[,2]-All.Pairs[,4]) + (All.Pairs[,5]-All.Pairs[,4]))) +
                                      ((epsilon^2)*Clones0.0))
      
      All.Pairs$"Po1,o2(0,1)"<- C* ((All.Pairs[,2]-All.Pairs[,4])+((All.Pairs[,2]-All.Pairs[,4])*epsilon)+(Clones0.0*epsilon))
      All.Pairs$"Po1,o2(1,0)"<- C* ((All.Pairs[,5]-All.Pairs[,4])+((All.Pairs[,5]-All.Pairs[,4])*epsilon)+(Clones0.0*epsilon))
      All.Pairs$"Po1,o2(0,0)"<- C * Clones0.0
      All.Pairs$"Po1"<- Po1
      All.Pairs$"Po2"<- Po2
      All.Pairs<- All.Pairs[c(3,4,1,2,6,7,8,9,10,11)]
      return(All.Pairs)
    }
    
    #Build Probability Matrix 
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

#Fucntion to calculate the MI tot
    
# hist_MI_tot <- function(results){
#   #hash out the steps below to stop generating histograms
# 
#   if (ncol(results)>2) { # if simulations were performed
#     Simulated_MItot<-numeric(ncol(results)-2)
#     for (i in 3:ncol(results)){
#       Simulated_MItot[i-2]<-sum(results[,i])
#     }
#   #   pdf("MI_distribution.pdf")
#   # hist(Simulated_MItot, main=expression('Histogram of Simulated I'['tot']),
#   #      xlab=expression('Simulated I'['tot']),
#   #      xlim=c((min(Simulated_MItot)-min(Simulated_MItot)*.1), sum(results[,2])),
#   #      col="grey")
#   # abline(v=sum(results[,2]), col="red", lwd=2)
#   # dev.off()
#    }
# return(Simulated_MItot)
# }
    

    

    ####Outputs####
    
Results <- Run_Sampling_for_MI(test, 10) 

# miTOT <-  hist_MI_tot(Results)
# hist_Results <- hist_MI_tot(Results)

#Write Out 
write.table(Results, "Pairwise_Mutual_Information.tsv", sep="\t")
# write(Simulated_MItot, "Simulated_Total_I_vector.tsv", sep="\t")



###Parallel GI Pairwise_MI ### 



#### Parallel_GI_pairwise_MI #### 
##############################

# DF<- read.delim(user_input[1], header=TRUE, stringsAsFactors=FALSE)
DF <- Results

simulations<- (ncol(DF)-2)
min_pval<- 1/simulations

Results_df<- DF[1:2]
Results_df$"mean simulated MI"<-NA
Results_df$"unadjusted P value"<- NA


for (i in 1:nrow(DF)){
  obsval<- as.numeric(DF[i,2])
  sim_Vector<- as.numeric(DF[i, 3:ncol(DF)])
  Results_df[i,3] <- mean(sim_Vector)
  x<- length(which(sim_Vector >= obsval))
  if (x == 0){
    Results_df[i,4]<- paste("<", min_pval)
  } 
  else {
    Results_df[i,4]<- x/simulations
  }
}


Results_df<- Results_df[order(Results_df[,4]), ]

write.table(Results_df, "Pairwise_MI_p_values.tsv", sep="\t")



