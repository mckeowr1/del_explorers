library(genetics) 
library(tidyverse)
library(rebus)


#Create a List of Matrices you want to loop through

file_names <- list.files(pattern = ".Genotype_Matrix.tsv")




for (f in file_names){ 
  #Read In the files
  gm <- read.table(glue::glue("/projects/b1059/projects/Ryan/dauer/LD/{f}"), header = T) 
  
  #Store name for writing out file  
  name <- f %>% str_extract(
    one_or_more(char_class(ALNUM)) %R% "_" %R% one_or_more(char_class(ALNUM)))
  
  #Prep for LD calc 
  if ( nrow(gm) > 1 ) {
    
    #Concat CHROM & POS column to make SNP_ID
    gm <- data.frame(snp_id = paste(gm$CHROM, gm$POS,
                                    sep = "_"), data.frame(gm)[, 5:ncol(gm)])
    #Create a blank list
    sn <- list()
    
    #Convert Binary GT to arbitrary Value
    for (i in 1:nrow(gm)) {
      sn[[i]] <- genetics::genotype(as.character(gsub(1, "T/T",
                                                      gsub(-1, "A/A", gm[i, 4:ncol(gm)]))))
    }
    
    
  }
  
  test <- data.frame(sn) #Convert back to a DF
  colnames(test) <- (gm$snp_id) 
  ldcalc <- t(genetics::LD(test)[["R^2"]]) #Return R^2 - could return other values see LD()
  diag(ldcalc) <- 1
  
  write.table(ldcalc, glue::glue("{name}_LD_between_paired_alleles_corr.tsv"), quote=F, row.names = T, col.names = NA, sep="\t")
}


