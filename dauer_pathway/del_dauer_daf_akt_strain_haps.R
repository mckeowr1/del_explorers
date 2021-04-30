#### Plot Strains x By Dauer Alleles ###
#### Plot Strains x Daf-2 & Akt-2 Alleles ###


####Subset data from flat file####

#Load Flat file
flat_file <- data.table::fread("WI.20210121.strain.annotation-GB-gene-impact-divergent-tempfix.tsv")

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


####Subset Strains with Variants in Daf and Akt-2####

#Filter Contingency Table for >1 Daf-2 and >1 Akt-2 variant
paired_strains_cont <- geneaff_by_strain %>%
  filter(.$`akt-2` >= 1 & .$`daf-2` >= 1 )

#Create a list of strains that have variants in daf and akt

paired_strains <- paired_strains_cont %>%
  rownames_to_column() %>%
  .$rowname

#Filter High Var dauer file for strains with >1 Daf-2 and >1 Akt-2
paired_dauer_haps <- duaer_goi_var_hi %>% filter(Strains %in% paired_strains)

####Create a contingency table: strains by variant (strains w Daf-2 and Akt-2 Var)####

#Prep new columns for the contingency table

prep_paired_dauer_haps <- paired_dauer_haps %>%
  mutate(CHROM.POS = str_c(.$CHROM, ":", .$POS)) %>%
  mutate(REF.ALT = str_c(.$REF, ":", .$ALT)) %>%
  mutate(CHROM.POS.R.A = str_c(.$CHROM.POS, ":", .$REF.ALT))

#Group by strains
gs_paired_dauer_haps <- prep_paired_dauer_haps %>%
  group_by(Strains) %>%
  summarise(var_present = (CHROM.POS.R.A))

#Create Contingency Table
strain_by_var <- table(gs_paired_dauer_haps$Strains,
                       gs_paired_dauer_haps$var_present) %>% as.data.frame.matrix()

#Convert to logical affected or unaffected
logical <- function(x){ ifelse( x > 0, 1 , 0) }

strain_by_var_log <- mutate_all(strain_by_var, logical) #%>%


####Plot Dauer Haplotype for strains >1 Daf-2 and >1 Akt-2 (All Dauer genes and just Daf & Akt-2)####

#Prep for ploting
prep_plot_strain_by_var <- strain_by_var_log %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(Column, Value, -rowname)

#Dauer Haplotype(Includes not just Daf-2 and Akt-2)
ggplot(prep_plot_strain_by_var, aes(x = Column, y = rowname, fill = Value)) +
  geom_tile() +
  # scale_fill_gradientn(name = "", colors = terrain.colors(3)) +
  scale_x_discrete(name = "") +
  scale_y_discrete(name = "") +
  theme(axis.text.x = element_text(angle = 90, size = 8),
        axis.text.y = element_text(size = 5)
  )

###Dauer Haplotype with just Daf-2 and Akt-2

#Subset Data for plotting

##Requires akt hi hap and daf hi hap
#Get list of Akt-hi-haplotypes & Daf-hi-haplotypes
akt_var <- akt_hi_hap %>%
  mutate(CHROM.POS.R.A = str_c(CHROM.POS, ":", REF.ALT)) %>%
  .$CHROM.POS.R.A

daf_var <- daf2_hi_hap %>%
  mutate(CHROM.POS.R.A = str_c(CHROM.POS, ":", REF.ALT)) %>%
  .$CHROM.POS.R.A

daf_akt <- c(akt_var, daf_var)

#Filter Strain by variant Logical to include columns in daf_var or akt_var

akt_daf_strain_log <- strain_by_var_log %>%
  as.data.frame() %>% select(c(daf_akt))

#Prep Data for plotting
prep_ak_daf_strain_log <- akt_daf_strain_log %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(Column, Value, -rowname)

ggplot(prep_ak_daf_strain_log, aes(x = Column, y = rowname, fill = Value)) +
  geom_tile() +
  # scale_fill_gradientn(name = "", colors = terrain.colors(3)) +
  scale_x_discrete(name = "") +
  scale_y_discrete(name = "") +
  theme(axis.text.x = element_text(angle = 90, size = 8),
        axis.text.y = element_text(size = 5)) +
  geom_vline(aes(xintercept = "X:14147803:T:C"), col = 'red', position = position_dodge(width = 10))


########## Strain with just daf-2 #####

unpaired_daf_paired_strains_cont <- geneaff_by_strain %>%
  filter(.$`akt-2` == 0 & .$`daf-2` >= 1 )

#Data Prep
unpaired_daf <- unpaired_daf_paired_strains_cont %>% mutate_all( logical) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(Column, Value, -rowname)


#Plot unpaired_daf
ggplot(unpaired_daf, aes(x = Column, y = rowname, fill = Value)) +
  geom_tile() +
  # scale_fill_gradientn(name = "", colors = terrain.colors(3)) +
  scale_x_discrete(name = "") +
  scale_y_discrete(name = "") +
  theme(axis.text.x = element_text(angle = 90, size = 8),
        axis.text.y = element_text(size = 5))


unpaired2_paired_strains_cont <- geneaff_by_strain %>%
  filter(.$`akt-2` >= 1 & .$`daf-2` >= 0 )
