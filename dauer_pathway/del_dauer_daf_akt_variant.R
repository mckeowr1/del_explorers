###Analyzes Variants Present in Daf-2 and Akt-2 #######

######IN-All Data Prep######
#Load Flat file
flat_file <- data.table::fread("WI.20210121.strain.annotation-GB-gene-impact-divergent-tempfix.tsv")

#Consequence Vectors
high_impact <- c("stop_gain", "start_lost", "splice_donor", "splice_acceptor", "frameshift", "stop_lost")
inframe_altering <- c( "inframe_deletion", "inframe_insertion",  "missense&inframe_altering")
moderate_impact <- c("missense", "moderate_deletion", "moderate_insertion", "missense&inframe_altering", "splice_region")
low_impact <- c("start_retained", "synonymous", "splice_region", "stop_retained")


#Load list of dauer genes
dauer_pathway <- data.table::fread("insulin_genes.tsv")

#Create List of dauer genes
d_genes <- dauer_pathway$ens_gene

#Subset Dauer Genes from flat file
dauer_var <- filter(flat_file, WORMBASE_ID %in% d_genes)

#Spread By Strains
dauer_var_s <- dauer_var %>% separate_rows(Strains, sep = ",")

####Create a Contingency table: strains by gene affected####

#Filter for "deleterious" variants
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

#Contingency Table of Paired Strains
paired_strains_cont <- geneaff_by_strain %>%
  filter(.$`akt-2` >= 1 & .$`daf-2` >= 1 )


#List of Paired Strains

paired_strains <- paired_strains_cont %>%
  rownames_to_column() %>%
  .$rowname


####Variant Analysis #####

#Filter daur_var file for just paired strains
paired_strains_dauer_var <- dauer_var_s %>% filter(Strains %in% paired_strains) #Filter for variants in paired strains

#Filter for variants in akt-2 and daf-2
goi <- c("WBGene00000103", "WBGene00000898") #Create list with akt and daf WBgeneID
paired_strains_goi_var <- paired_strains_dauer_var %>% filter(WORMBASE_ID %in% goi) #Filter Paired Strains Dauer Var DF

#Filter For 'Deleterious' Variants
paired_strains_goi_hi_var <- paired_strains_goi_var %>%
  filter(CONSEQUENCE == "missense" & BLOSUM <= -2 |
           CONSEQUENCE %in% inframe_altering |
           CONSEQUENCE %in% high_impact )

## Akt-2 Analysis ##
akt2_hi <- paired_strains_goi_hi_var %>%
  filter(WORMBASE_ID == "WBGene00000103") %>%
  mutate(CHROM.POS = str_c(.$CHROM, ":", .$POS)) %>%
  mutate(REF.ALT = str_c(.$REF, ":", .$ALT))

akt2_hi_by_var <- akt2_hi %>%
  group_by(CHROM.POS, REF.ALT) %>% #Group By Position and Allele
  nest() %>%
  mutate(Cons = sapply(data, prin_con)) %>%  #Print the predicted consequence for each haplotype
  mutate(TRANSCRIPT_DISTRIBUTION = sapply(data, transcript_check)) %>% #if a transcript has one transcript with more variants
  mutate(TRANSCRIPTS = sapply(data, prin_transcripts)) %>% #Print the transcripts to a string
  mutate(Strains_AFFECTED = sapply(data, prin_strains)) %>% #Print Strains Affected
  mutate(DIVERGENT = sapply(data, divergent_check)) %>% #If divergent region occurs in any Isotype at that Pos
  mutate(Num_Strains = sapply(data, count_strains)) #Number of Paired Strains Affected by the variant

## Daf-2 Analysis ##

daf2_hi <- paired_strains_goi_hi_var %>%
  filter(WORMBASE_ID == "WBGene00000898") %>%
  mutate(CHROM.POS = str_c(.$CHROM, ":", .$POS)) %>%
  mutate(REF.ALT = str_c(.$REF, ":", .$ALT))


daf_hi_by_var <- daf2_hi %>%
  group_by(CHROM.POS, REF.ALT) %>%
  nest() %>%
  mutate(Cons = sapply(data, prin_con)) %>%  #Print the predicted consequence for each haplotype
  mutate(TRANSCRIPT_DISTRIBUTION = sapply(data, transcript_check)) %>%
  mutate(TRANSCRIPTS = sapply(data, prin_transcripts)) %>% #Check if a transcript has one transcript with more variants
  mutate(Strains_AFFECTED = sapply(data, prin_strains)) %>%
  mutate(DIVERGENT = sapply(data, divergent_check)) %>%
  mutate(Num_Strains = sapply(data, count_strains))


## Variant Analysis Functions ##

#Print the transcript isoforms
prin_transcripts <- function(df){ toString(unique(df$TRANSCRIPT))}
#Print the consequences for that haplotype
prin_con <- function(df){ toString(unique(df$CONSEQUENCE))}
#Function to pull strains with that variant
prin_strains <- function(df){toString(unique((df$Strains)))}
#Function to check if the variant is in divergent region
divergent_check <- function(df){unique(df$DIVERGENT)}
#Function to count the number of strains
count_strains <- function(df){length(unique(df$Strains))}
