library(dplyr)
library(tidyr)
library(stringr)
library(data.table)

dauer_pathway <- data.table::fread("/Users/ryanmckeown/documents/andersen_lab/data_exploration/del_explorers/insulin_genes.tsv")
d_genes <- dauer_pathway$ens_gene


#Load Flat File
flat_file <- data.table::fread("/Users/ryanmckeown/documents/andersen_lab/data_exploration/del_explorers/WI.20210121.strain-annotation.bcsq.tsv")

#Test File - Not really great for stats test
#test <- dplyr::sample_n(flat_file, 1000) #Need a much larger sample size to make test work 

#Impact Classifications
#Used for filtering
high_impact <- c("stop_gain", "start_lost", "splice_donor", "splice_acceptor", "frameshift", "stop_lost")
inframe_altering <- c( "inframe_deletion", "inframe_insertion",  "missense&inframe_altering")
moderate_impact <- c("missense", "moderate_deletion", "moderate_insertion", "missense&inframe_altering", "splice_region")
low_impact <- c("start_retained", "synonymous", "splice_region", "stop_retained")



rnas <- c("antisense", # List of Non-coding biotypes
          "IG_C_gene",
          "IG_D_gene",
          "IG_J_gene",
          "IG_LV_gene",
          "IG_V_gene",
          "lincRNA",
          "macro_lncRNA",
          "miRNA",
          "misc_RNA",
          "Mt_rRNA",
          "Mt_tRNA",
          "polymorphic_pseudogene",
          "ribozyme",
          "rRNA",
          "sRNA",
          "scRNA",
          "scaRNA",
          "sense_intronic",
          "sense_overlapping",
          "snRNA",
          "snoRNA",
          "TR_C_gene",
          "TR_D_gene",
          "TR_J_gene",
          "TR_V_gene") 


#Get Flat File Ready
prepped_ff <- test %>% dplyr::filter(!BIOTYPE %in% rnas) %>% #Filter out Non-coding biotypes
  dplyr::filter(!stringr::str_detect(CONSEQUENCE, "@")) %>%  #Filter out linker variants 
  tidyr::separate_rows(Strains, sep = ",") #Split up all the strains

#Create a list of strains to analyze 
#Strains w/ >1 Dauer Variant
dauer_var <- dplyr::filter(flat_file, WORMBASE_ID %in% d_genes) %>% 
  tidyr::separate_rows(Strains, sep = ",")

#Filtering Strains w/  high-"er" impact variants & Get Format Ready for Contingency table
  #Intermediate step to get to Paired Strains list ? 
duaer_goi_var_hi <- dauer_var %>%
  dplyr::filter(CONSEQUENCE == "missense" & BLOSUM <= -2 |
           CONSEQUENCE %in% inframe_altering |
           CONSEQUENCE %in% high_impact ) %>% 
  dplyr::group_by(Strains) %>% 
  dplyr::summarise(genes_aff = (GENE)) #Format for Contingency table

#Make Contingency table, #Var in each gene by strain
geneaff_by_strain <- table(duaer_goi_var_hi$Strains,
                           duaer_goi_var_hi$genes_aff) %>% as.data.frame.matrix

#Filter for strains with paired variants in dauer genes 
paired_strains_cont <- geneaff_by_strain %>% #Create contingency table with strains and # of genes w var
  dplyr::filter(.$`akt-2` >= 1 & .$`daf-2` >= 1 ) #Filter for stains with daf-2 , akt-2 pariing

#Create List of Strains to iterate through 
paired_strains <- paired_strains_cont %>% 
  tibble::rownames_to_column() %>%
  .$rowname #Select the strain names


#List of Strains To go through 
  #Strains w/ 
    # - High-"er" Variation in Duaer genes 
    # - Strains with variation in daf and akt2 

#Call it stains list for some reason
strains_list <- paired_strains


### T-test function to compare mean var per Duaer and Non Dauer gene ### 
#######################################################################

#Requires a strain list created with whatever parameters you want

#Set up for loop
dauer_strain_mean_var_by_gene = NULL #Blank Vector for all variants 
del_dauer_strain_mean_var_by_gene = NULL #Blank Vector for del variants 
del_dauer_t_test = NULL #T-test vector
pb = txtProgressBar(min = 0, max = length(strains_list), initial = 0) #Fun progress bar
stepi = 0




for (strain in strains_list){ 
  filtered_strains <- test %>% dplyr::filter(Strains == strain) #Filter for strain in list

#Mean Variants per Dauer Gene
  dauer_genes <- filtered_strains %>% 
    dplyr::filter(WORMBASE_ID %in% d_genes) %>%
    dplyr::group_by(CHROM, POS, REF, ALT, WORMBASE_ID) %>% 
    tidyr::nest() %>%  
    dplyr::group_by(WORMBASE_ID) %>% 
    dplyr::summarise(n = n())
  
  mean_dauer <- mean(dauer_genes$n)
  std_dauer <- sd(dauer_genes$n)
  
#Mean Del Variants per Dauer Gene  
  del_dauer <- filtered_strains %>% 
    dplyr::filter(WORMBASE_ID %in% d_genes) %>%
    dplyr::filter(CONSEQUENCE == "missense" & BLOSUM <= -2 |
             CONSEQUENCE %in% inframe_altering |
             CONSEQUENCE %in% high_impact ) %>% 
    dplyr::group_by(CHROM, POS, REF, ALT, WORMBASE_ID) %>% 
    tidyr::nest() %>%  
    dplyr::group_by(WORMBASE_ID) %>% 
    dplyr::summarise(n = n())

    
  dauer_genes_w_del_var <- nrow(del_dauer) 
  mean_del_dauer <- mean(del_dauer$n) 
  std_del_dauer <- sd(del_dauer$n)

#Mean Varaints per Non-Dauer Genes 
  nd_genes <- filtered_strains %>% 
    dplyr::filter(!WORMBASE_ID %in% d_genes) %>% 
    dplyr::group_by(CHROM, POS, REF, ALT, WORMBASE_ID) %>% 
    tidyr::nest() %>% 
    dplyr::group_by(WORMBASE_ID) %>% 
    dplyr::summarise(n = n())
  
  mean_non_dauer <- mean(nd_genes$n)
  std_non_dauer <- sd(nd_genes$n)

#Mean del Variants per Dauer genes 
  del_nd <- filtered_strains %>% 
    dplyr::filter(!WORMBASE_ID %in% d_genes) %>% 
    dplyr::filter(VARIANT_IMPACT == "HIGH" |
             CONSEQUENCE == "missense" & BLOSUM <= -2) %>% 
    dplyr::group_by(CHROM, POS, REF, ALT, WORMBASE_ID) %>% 
    tidyr::nest() %>% 
    dplyr::group_by(WORMBASE_ID) %>% 
    dplyr::summarise(n = n())
  
  non_dauer_genes_w_del_var <- nrow(del_nd)  
  mean_del_non_dauer <- mean(del_nd$n) #Mean number of del variants per gene
  std_del_non_dauer <- sd(del_nd$n) #Mean number of del variants per gene
  
  ###Run Stats###
  
  #T-test all genes
  test_all <- t.test(dauer_genes$n , nd_genes$n)
  t_stat_all <- test_all$statistic #Extract the test statistic
  t_pval_all <- test_all$p.value #Extract the P-Value
  t_ci_all <- test_all$conf.int #Extract the confidence interval 
  
  
  #All Gene output
  
  dauer_strain_mean_var_by_gene = rbind(dauer_strain_mean_var_by_gene, 
                                        data.frame(strain,
                                                   mean_dauer, 
                                                   std_dauer,
                                                   mean_non_dauer,
                                                   std_non_dauer, 
                                                   t_stat_all, 
                                                   t_pval_all#, 
                                                   #t_ci_all #CI causes duplcate rows, leave out unless 
                                                    #needed
                                        ))
  
  
#If there is not >1 Del dauer var, can't run t-test
    #Shouldnt be a problem with paired strains list, but good check for other inputs
  
  if(nrow(del_dauer) > 1 ){
    #T-test Del genes    
    test_del <- t.test(del_dauer$n , del_nd$n)
    t_stat_del <- test_del$statistic #Extract the test statistic
    t_pval_del <- test_del$p.value #Extract the P-Value
    t_ci_del <- test_del$conf.int #Extract the confidence interval 
    
    del_dauer_t_test = rbind(del_dauer_t_test,
                             data.frame(strain,
                                        t_stat_del, 
                                        t_pval_del,
                                        t_ci_del
                             )) 
    
  }
  
  #Del Gene output 
  del_dauer_strain_mean_var_by_gene = rbind(del_dauer_strain_mean_var_by_gene, 
                                            data.frame(strain,
                                                       dauer_genes_w_del_var,
                                                       mean_del_dauer, 
                                                       std_del_dauer, 
                                                       non_dauer_genes_w_del_var,
                                                       mean_del_non_dauer, 
                                                       std_del_non_dauer))
  
  #Progress Bar
  stepi = stepi + 1 
  setTxtProgressBar(pb,stepi)
  
}



#Write out files 

del_t_test_clean <- del_dauer_t_test %>% select(-t_ci_del)

#Joing Duaer T-test to DF
del_dauer_t <- left_join(del_dauer_strain_mean_var_by_gene, del_t_test_clean, by = "strain")


write.csv(dauer_strain_mean_var_by_gene, "dauer_strain_mean_var_by_gene.csv" )
write.csv(del_dauer_t, "del_dauer_strain_mean_var_by_gene.csv")


#Load Files for stats 
dauer_stats <- data.table::fread("dauer_strain_mean_var_by_gene.csv")
del_dauer_stats<- data.table::fread("del_dauer_strain_mean_var_by_gene.csv")


## Dauer Stats ## - Just Filtering Through the results

higher_d <- dauer_stats %>% filter(mean_dauer > mean_non_dauer) %>% filter (t_pval_all < 0.05)
del_higher_d <- del_dauer_stats %>% filter(mean_del_dauer > mean_del_non_dauer)







# #Troubleshooting 
# 
# #Filter for the strain
# df1 <- prepped_ff %>% filter(Strains =="DL238")
# 
# #Get list of dauer variants
# df2 <- df1 %>% 
#   filter(WORMBASE_ID %in% d_genes) #%>%
# group_by(CHROM, POS, REF, ALT, WORMBASE_ID) %>% 
#   nest() %>%  
#   group_by(WORMBASE_ID) %>% summarise(n = n())
# 
# 
# #Error in del_t-test
# 
# del_df1 <- df1 %>% 
#   filter(WORMBASE_ID %in% d_genes) #%>%
# filter(CONSEQUENCE == "missense" & BLOSUM <= -2 | 
#          CONSEQUENCE %in% inframe_altering | 
#          CONSEQUENCE %in% high_impact ) %>% 
#   group_by(CHROM, POS, REF, ALT, WORMBASE_ID) %>% 
#   nest() %>%  
#   group_by(WORMBASE_ID) %>% summarise(n = n())
# 
# del_df2 <- df1 %>% 
#   filter(!WORMBASE_ID %in% d_genes) %>% 
#   filter(VARIANT_IMPACT == "HIGH" |
#            CONSEQUENCE == "missense" & BLOSUM <= -2) %>% 
#   group_by(CHROM, POS, REF, ALT, WORMBASE_ID) %>% 
#   nest() %>% 
#   group_by(WORMBASE_ID) %>% summarise(n = n())
# 
# 
# t.test(del_df1$n, del_df2$n)
# 
# 
# akt_var <- flat_file %>% filter(GENE == "akt-2")
# 
# akt_var_prep <- prepped_ff %>% filter(GENE == "akt-2")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #### Test Code #### 
# #Collapse by transcript? 
# 
# #Group by everything, collapse dataframe
# collapsed <- filtered %>% group_by(CHROM, POS, REF, ALT, WORMBASE_ID) %>% nest()
# 
# 
# #Group by Genes
# genes <- collapsed %>% group_by(WORMBASE_ID) %>% summarise(n = n())
# 
# #Calculate mean 
# 
# 


