library(ggplot2)
library(tidyverse)
library(data.table)
library(plotly)

corr <- data.table::fread("/Users/ryanmckeown/Documents/andersen_lab/data_exploration/del_explorers/dauer_pathway/LD/LD_between_paired_alleles_corr.tsv")  



#Data prep
plotable <- corr %>% 
  as_tibble() %>% 
  gather(key = "Y", value = "Z", -1) %>%  
  mutate(Z = round(Z, digits = 4)) #%>%  
  rename("SNP1" = "V1", 
         "SNP2" = "Y", 
         "LD_COR" = "Z")

#Viz 
ggplot(plotable, aes(V1, Y, fill = Z)) + 
  geom_tile() + 
  geom_text(aes(label = Z))

#Overlay LD on paired variants graph
allele_pair<- data.table::fread("/Users/ryanmckeown/Documents/andersen_lab/data_exploration/del_explorers/dauer_pathway/LD/daf_akt_allele_paire_freq.tsv")

#Replace NA's in plotable so that I can join regardless order of SNP1 and SNP2. All combinations in real data will have a value

for (i in 1:nrow(plotable)){
   if(is.na(plotable[i,3])) {
    SNP1 <- as.character(plotable[i,1])
    SNP2 <- as.character(plotable[i,2])
    val_df <- plotable %>% filter(V1 == SNP2 & Y == SNP1) 
    val <- as.numeric(val_df$Z)
    plotable[i, 3] <- val


    }

}




test_join <- left_join(plotable, allele_pair , by = c("V1" = "SNP1", "Y" = "SNP2"))

#Check Test_Join By Plotting with allele freq as fill 
#Should Recapitulate 

ggplot(test_join, aes(V1, Y, fill = Freq)) + 
  geom_tile() + 
  geom_text(aes(label = Z))



                               

