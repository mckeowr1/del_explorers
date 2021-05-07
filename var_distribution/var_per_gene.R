library(dplyr)
library(data.table)

dauer_pathway <- data.table::fread("insulin_genes.tsv")

d_genes <- dauer_pathway$ens_gene


flat_file <- data.table::fread("WI.20210121.strain.annotation-GB-gene-impact-divergent-tempfix.tsv")

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

nd_genes <- flat_file %>%
  filter(!WORMBASE_ID %in% dauer_var) %>%
  filter(BIOTYPE %in% rnas) %>% #Filter out Non-coding biotypes
  filter(!str_detect(CONSEQUENCE, "@")) %>% #Filter out linker variants
  separate_rows(Strains, sep = ",")

d_genes <- flat_file %>%
  filter(WORMBASE_ID %in% dauer_var) %>%
  filter(BIOTYPE %in% rnas) %>% #Filter out Non-coding biotypes
  filter(!str_detect(CONSEQUENCE, "@")) %>% #Filter out linker variants
  separate_rows(Strains, sep = ",") #Make each strain a row


# Strains_list <- unique(filtered$Strains) #Now must be uniqe to stats being run
Strains_list <- c("AB1", "N2", "CB4856")

var_per_nd_gene <- NULL #Empty vector to store for results
pb = txtProgressBar(min = 0, max = length(Strains_list), initial = 0)
stepi = 0


#Need to have this set up for deleterious variants too
for(strain in Strains_list) {
  filtered_f <- nd_genes %>% filter(Strains == strain)
  collpased_f <- filtered_f %>% group_by(CHROM, POS, REF, ALT, WORMBASE_ID) %>% nest()
  genes_f <- collpased_f %>% group_by(WORMBASE_ID) %>% summarise(n = n())
  mean <- mean(genes_f$n)
  std <- sd(genes_f$n)
  var_per_gene = rbind(var_per_gene, data.frame(strain, mean, std))
  stepi = stepi + 1
  setTxtProgressBar(pb,stepi)
}

write.csv(var_per_nd_gene, "var_per_nd_gene.csv")

print("Not dauer Mean calculated")





for(strain in Strains_list){
  filtered_f <- d_genes %>% filter(Strains == strain)
  collpased_f <- filtered_f %>% group_by(CHROM, POS, REF, ALT, WORMBASE_ID) %>% nest()
  genes_f <- collpased_f %>% group_by(WORMBASE_ID) %>% summarise(n = n())
  mean <- mean(genes_f$n)
  std <- sd(genes_f$n)
  var_per_gene = rbind(var_per_gene, data.frame(strain, mean, std))
  stepi = stepi + 1
  setTxtProgressBar(pb,stepi)
}


write.csv(var_per_d_gene, "var_per_d_gene.csv")

print("Dauer Mean Calculated")



#Filter Linkers and non-coding var from flat file
prepped_ff <- flat_file %>% filter(!BIOTYPE %in% rnas) %>% #Filter out Non-coding biotypes
  filter(!str_detect(CONSEQUENCE, "@")) #Filter out linker variants


#Comparing the mean occurs with all the data

var_pergene_dauervar_strains = NULL

#Strains list should just include the strains that have dauer var

dauer_var <- prepped_ff %>% (WORMBASE_ID %in% d_genes) %>%


  strains_list <-

  for (strain in strains_list){
    filtered_strains <- prepped_ff %>% filter(Strains == strain) #FF needs strains separated first
    #Grab daure and non dauer genes & Make Gene Groupings

    dauer_genes <- filtered_strains %>%
      filter(WORMBASE_ID %in% d_genes) %>%
      group_by(CHROM, POS, REF, ALT, WORMBASE_ID) %>%
      nest() %>%
      group_by(WORMBASE_ID) %>% summarise(n = n())

    mean_dauer <- mean(dauer_genes$n)
    std_dauer <- sd(dauer_genes$n)


    nd_genes <- filtered_strains %>%
      filter(!WORMBASE_ID %in% d_genes) %>%
      group_by(CHROM, POS, REF, ALT, WORMBASE_ID) %>%
      nest() %>%
      group_by(WORMBASE_ID) %>% summarise(n = n())

    mean_non_dauer <- mean(nd_genes$n)
    std_dauer <- sd(nd_genes$n)


    #T-test
    test <- t.test(dauer_genes$n , nd_genes$n)
    t_stat <- test$statistic #Extract the test statistic
    t_pval <- test$p.value #Extract the P-Value
    t_ci <- test$conf.int #Extract the confidence interval

    by_strain = rbind(by_strain, mean_dauer , std_dauer,
                      mean_non_dauer, std_non_dauer, t_stat,
                      t_pval, t_ci)
    stepi = stepi + 1
    setTxtProgressBar(pb,stepi)

  }





#### Test Code ####
#Collapse by transcript?

#Group by everything, collapse dataframe
collapsed <- filtered %>% group_by(CHROM, POS, REF, ALT, WORMBASE_ID) %>% nest()


#Group by Genes
genes <- collapsed %>% group_by(WORMBASE_ID) %>% summarise(n = n())

#Calculate mean 
