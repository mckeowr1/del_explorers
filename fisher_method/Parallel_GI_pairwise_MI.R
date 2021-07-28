#!/usr/bin/env Rscript


user_input=commandArgs(trailingOnly = TRUE)

if (length(user_input)<1){
  stop("Pairwise_Mutual_Information.tsv must be provided as an argument", call=FALSE)
} 

DF<- read.delim(user_input[1], header=TRUE, stringsAsFactors=FALSE)
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
