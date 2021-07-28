#!/usr/bin/env Rscript

#run by calling Rscript --vanilla Parallel_GI_Itot.R raw_data.tsv desired_number_simulations(supplied as integer)
#if no simulations are desired enter an argument of 0 for simulations

#raw input data should be tab seperated (tsv) file with clones as columns and gene names in rows in format below: 
#Ex. 
#   clone1name  clone2name  clone3name ...
#   genex       geney       genez
#   geney       genex       genex

#output will be written as .tsv files and .pdf figures. 

#required packages: 
#plyr
#pracma
######

#takes trailing arguments 
user_input=commandArgs(trailingOnly = TRUE)
 
if (length(user_input)< 2){
  stop("Raw data and number of desired replicates must be supplied as arguments", call=FALSE)
} 

source("scripts/Construct_Pairs.R")
source("scripts/Probability_Sample.R")
source("scripts/Generate_Pairs_Dataframe.R")
source("scripts/Build_Probability_Matrix.R")
source("scripts/Run_sampling_for_Mutual_Information.R")

library("plyr")
library("pracma")

DF<- read.delim(user_input[1], header=TRUE, stringsAsFactors=FALSE)
Results<-  Run_Sampling_for_MI(DF, user_input[2])

#hash out the steps below to stop generating histograms 
# 
# if (ncol(Results)>2) { # if simulations were performed
#   Simulated_MItot<-numeric(ncol(Results)-2)
#   for (i in 3:ncol(Results)){
#     Simulated_MItot[i-2]<-sum(Results[,i])
#   }
#   pdf("MI_distribution.pdf")
# hist(Simulated_MItot, main=expression('Histogram of Simulated I'['tot']), 
#      xlab=expression('Simulated I'['tot']), 
#      xlim=c((min(Simulated_MItot)-min(Simulated_MItot)*.1), sum(Results[,2])), 
#      col="grey")
# abline(v=sum(Results[,2]), col="red", lwd=2)
# dev.off()
# }

write.table(Results, "Pairwise_Mutual_Information.tsv", sep="\t")
write(Simulated_MItot, "Simulated_Total_I_vector.tsv", sep="\t")

