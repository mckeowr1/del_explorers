library(genetics)
library(tidyverse)
#library(glue)
gm <- read.table("del_explorers/dauer_pathway/LD/Genotype_Matrix_t.tsv", header = T) %>%  
  mutate(ID = str_c(.$CHROM, "*", .$POS))

test_df <- gm[2,1:20]


# 
# test <- genetics::genotype(select(gm, 5:544))
# 
# check <- genotype(as.character(gm[18]))

#Convert each cell to genotype object

# genotype(as.character(select(gm, 5:544))) 
# 
# select(gm, contains(".")) #Select just the columns with alleles


#Convert to -1 and 1 

Ref <- test_df[1,3]

str_detect(gm[1,5] , Ref )

#Works but based on idex 
replace(test_df, 5:10, 1 )


#How to iterate over each row for each column

for(i in nrow(test_df)){
  ref <- test_df[i,3]
  alt <- test_df[i,4]
  print(test_df[i,])

  #Will break with missing genotypes
  for(c in 1:ncol(test_df[i,])) { 
    print(test_df[i, c])
      if(test_df[i, c] == glue("{ref}/{ref}")) {
      test_df[i, c] = 1 } 
      if(test_df[i,c] == glue("{alt}/{alt}")) {
        test_df[i, c] = -1 }
    
     #  replace(test_df, c , +1)} else {
     #    replace(test_df, c, -1)}
     # 
     #   
    
  
  #If we see the ref then replace
  
}

}

str_replace_all(test_df, Ref, '1')

print(test_df)


sn <- list()

for (i in 1:nrow(gm)){ 
  
  
  sn[[i]]<- genotype(as.character(select(gm, contains("."))))
  

  }

test <- data.frame(sn)
colnames(test) <- (gm$ID) #Name the column with the ID of the snps
ldcalc <- t(genetics::LD(test)[[4]])^2
diag(ldcalc) <- 1





# for (i in 1:nrow(gm)) {
#   #Convert each column to genotype object
#   sn[[i]] <- genetics::genotype(select(gm, 5:545[,1]))
# }

# test <- data.frame(sn)
# colnames(test) <- gm$ID
