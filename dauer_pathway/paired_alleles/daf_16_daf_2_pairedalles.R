library(data.table)
library(reshape2)
library(tidyverse)

#Impact Classifications
high_impact <- c("stop_gain", "start_lost", "splice_donor", "splice_acceptor", "frameshift", "stop_lost")
inframe_altering <- c( "inframe_deletion", "inframe_insertion",  "missense&inframe_altering")
moderate_impact <- c("missense", "moderate_deletion", "moderate_insertion", "missense&inframe_altering", "splice_region")
low_impact <- c("start_retained", "synonymous", "splice_region", "stop_retained")



flat_file <- data.table::fread("/projects/b1059/projects/Ryan/WI.20210121.strain-annotation.bcsq.tsv")
dauer_pathway <- data.table::fread("/projects/b1059/projects/Ryan/dauer/insulin_genes.tsv")
d_genes <- dauer_pathway$ens_gene

#Filter for dauer genes
dauer_var <- dplyr::filter(flat_file, WORMBASE_ID %in% d_genes)

#Exapand Strain Column
dauer_var_s <- dauer_var %>% tidyr::separate_rows(Strains, sep = ",")


#Filter for "high" impact
duaer_goi_var_hi <- dauer_var_s %>%
  filter(CONSEQUENCE == "missense" & BLOSUM <= -2 |
           CONSEQUENCE %in% inframe_altering |
           CONSEQUENCE %in% high_impact )


#group by strains
gs_dauer_hi <- duaer_goi_var_hi %>%
  group_by(Strains) %>%
  summarise(genes_aff = (GENE))

#Create Contingency Table
geneaff_by_strain <- table(gs_dauer_hi$Strains,
                           gs_dauer_hi$genes_aff) %>% as.data.frame.matrix


#Filter Contingency Table for >1 Daf-2 and >1 Akt-2 variant
paired_strains_cont <- geneaff_by_strain %>%
  filter(.$`daf-18` >= 1 & .$`akt-2` >= 1 )

paired_strains <- paired_strains_cont %>%
  rownames_to_column() %>%
  .$rowname

#Filter High Var dauer file for strains with >1 Daf-2 and >1 Akt-2
paired_dauer_haps <- duaer_goi_var_hi %>% filter(Strains %in% paired_strains) %>%  
  filter(GENE == "daf-18" | GENE == "akt-2" )

prep_paired_dauer_haps <- paired_dauer_haps %>%
  mutate(CHROM.POS = str_c(.$CHROM, "_", .$POS)) #%>%
# mutate(REF.ALT = str_c(.$REF, ":", .$ALT)) %>%  
# mutate(CHROM.POS.R.A = str_c(.$CHROM.POS, ":", .$REF.ALT))

#Group by strains
gs_paired_dauer_haps <- prep_paired_dauer_haps %>%
  group_by(Strains) %>%
  summarise(var_present = (CHROM.POS))

w <- dcast(gs_paired_dauer_haps, var_present~Strains)
x <- as.matrix(w[,-1]) #Convert to a matrx 
x[is.na(x)] <- 0
x <- apply(x, 2,  function(x) as.numeric(x > 0)) 
v_f <- x %*% t(x)
diag(v_f) <- 0  #Set Diagonal to Zero
dimnames(v_f) <- list(w[, 1], w[,1])
v_f

longData_f <- melt(v_f)
longData_f<-longData_f[longData_f$value!=0,]

ggplot(longData_f, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill=value)) + geom_text(aes(label = value))


#Print out Paired Alleles


#Daf16 & Daf2 

longData_f %>% rename("Freq" = "value", "SNP1" = "Var1", "SNP2" = "Var2" ) %>%  
  write_tsv("daf18_akt2_allele_paired_freq.tsv")



daf18_daf2_allles <- longData_f %>% 
  rename("Freq" = "value", "SNP1" = "Var1", "SNP2" = "Var2" ) %>%
  select(SNP1) %>%  
  separate(SNP1, into = c("CHROM" , "POS"), sep = "_") %>% 
  unique() %>%  
  write_tsv("daf18_akt2_alleles", col_names = FALSE)








