library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(gplots)
library(graphics)
library(reshape2)
library(igraph)

flat_file <- data.table::fread("WI.20210121.strain.annotation-GB-gene-impact-divergent-tempfix.tsv")
dauer_pathway <- data.table::fread("insulin_genes.tsv")
#Make a Vector of Dauer Pathway Genes
d_genes <- dauer_pathway$ens_gene

#Filter For Dauer Variants
dauer_var <- filter(flat_file, WORMBASE_ID %in% d_genes)

#Function to test CHI-Squared Distribution of Variants in genes
#Takes Flat File Filtered by genes

###Load function to do analysis - Outputs Chi Squared Contingency Table, Pairs Table, Pairs Plot

```{r}
var_distribution <- function(df, consequences){
#Data Transformations
  filtered <- df %>%
  filter(CONSEQUENCE %in% consequences) #Filter by input consequence

  strain_split <- filtered %>%
  separate_rows(Strains, sep = ",") #Parse strain rows

  strain_group <- strain_split %>% #List strain in 1 col and gene_aff in 1 col
    group_by(Strains) %>%
    summarise(genes_aff = (GENE))

  contingency <- table(strain_group$Strains, #1 row = 1 strain
                       strain_group$genes_aff) %>% #1 col = 1 gene
                       as.data.frame.matrix()
#Chi-Squared Tests

  chi_table <- cbind(contingency, t(apply(contingency, 1, function(x) {
    ch <- chisq.test(x)
    c(unname(ch$statistic), ch$p.value)})))
  #colnames(df1)[15:16] <- c('x-squared', 'p-value') #The number of columns will change
  #Add Some formatting?

#Chi-Squared Formatting



#Looking at Pairings

logical <- function(x){ ifelse( x > 0, 1 , 0) }
logical_cont <- mutate_all(contingency, logical) # Convert to Affected/Unaffected logic
logical_cont$genes_affected <- rowSums(logical_cont) #Sum number of genes affected in each strain
pairs <- filter(logical_cont, genes_affected > 1) #Subset strains w/ >1 gene affected



#Heatmap of pairing by strains
pairs_plot <- pairs %>%
  as.data.frame() %>%
  select(-genes_affected) %>% #Remove sum column
  rownames_to_column() %>% #Convert matrix rows to column names
  gather(Column, Value, -rowname) #Older version of pivot longer, makes matrix plotable

plot <- ggplot(pairs_plot, aes(x = rowname, #Strain
                               y = Column, #Gene
                               fill = Value)) + #Affected or unaffected
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90))


#Co-Occurence Matrix

melted_data <- melt(strain_group) #wide to long reshaping, may not be needed
w <- dcast(melted_data , genes_aff~Strains) #Long to wide reshaping, strains as cols, genes as row
x <- as.matrix(w[,-1]) #Convert to a matrx w/o 1st column
x[is.na(x)] <- 0 #Replace any NA with zeros
x <- apply(x, 2,  function(x) as.numeric(x > 0)) #Convert to logical data - Is this a silly step?
v_f <- x %*% t(x) #multiply matrix by itself transposed
diag(v_f) <- 0  #Set Diagonal to Zero
dimnames(v_f) <- list(w[, 1], w[,1]) #Lable the matrix
v_f


longData_f <- melt(v_f)
longData_f<-longData_f[longData_f$value!=0,]

v_plot <- ggplot(longData_f, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill=value)) + geom_text(aes(label = value))





# melted_data <- melt(strain_group) #what does melt do?
# w <- dcast(melted_data, genes_aff~Strains) # Make data go wider by gene's affected
# x <- as.matrix(w[,-1]) #Convert to a matrx
# x[is.na(x)] <- 0
# x <- apply(x, 2,  function(x) as.numeric(x > 0))
# co_table <- x %*% t(x)
# diag(co_table) <- 0  #Set Diagonal to Zero
# dimnames(co_table) <- list(w[, 1], w[,1])
#
# co_table <- co_table %>% as.data.frame.matrix()


#Plot Co-occurence matrix

# longData <- melt(co_table)
# longData<-longData[longData$value!=0,]
#
# co_plot <- ggplot(longData, aes(x = Var2, y = Var1)) +
#   geom_tile(aes(fill=value)) + geom_text(aes(label = value))


#List for outputs
out <- list()
  out$chi_table <- chi_table #Return Contingency table w/ Chi-squared values
  out$paired_variants <- pairs #Return Strain by gene table w/ >1 gene affected
  out$plot <- plot #Retirn Strain by gene table above plotted
  out$co_plot <- v_f #Return Co-occurence matrix as DF
  out$co_occurence <- v_plot #Return plot of co-occurence matrix


  return(out)

}





chi_high <- var_distribution(dauer_var, high_impact)

high_impact <- c("stop_gain", "start_lost", "splice_donor", "splice_acceptor")
moderate_impact <- c( "inframe_deletion", "inframe_insertion",  "missense&inframe_altering")
high_moderate_impact <- c(high_impact , moderate_impact)


#Why does Chi_high moderate only return variants
chi_high <- var_distribution(dauer_var, high_impact)
chi_moderate <-var_distribution(dauer_var, moderate_impact)
chi_high_moderate <- var_distribution(dauer_var, high_moderate_impact)

test_out<- chi_moderate$paired_variants



#Sanity Check - Correct Number of genes
high_vars <- filter(dauer_var, CONSEQUENCE %in% high_impact)

print(unique(high_vars$GENE))
