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
  filtered <- df %>% filter(CONSEQUENCE %in% consequences)

  strain_split <- filtered %>% separate_rows(Strains, sep = ",")

  strain_group <- strain_split %>% group_by(Strains) %>%
    summarise(genes_aff = (GENE))

  contingency <- table(strain_group$Strains,
                       strain_group$genes_aff) %>% as.data.frame.matrix()
#Chi-Squared Tests

  chi_table <- cbind(contingency, t(apply(contingency, 1, function(x) {
    ch <- chisq.test(x)
    c(unname(ch$statistic), ch$p.value)})))
  #colnames(df1)[15:16] <- c('x-squared', 'p-value') #The number of columns will change
  #Add Some formatting?

#Chi-Squared Formatting



#Looking at Pairings

logical <- function(x){ ifelse( x > 0, 1 , 0) }
logical_cont <- mutate_all(contingency, logical) #%>%
logical_cont$genes_affected <- rowSums(logical_cont)
pairs <- filter(logical_cont, genes_affected > 1)

#pairs <-  mutate(sum = rowSums(logical_cont))
 #%>% filter(sum > 1)

#Heatmap of pairing by strains
pairs_plot <- pairs %>% as.data.frame() %>% select(-genes_affected) %>%
  rownames_to_column() %>%
  gather(Column, Value, -rowname)

plot <- ggplot(pairs_plot, aes(x = rowname, y = Column, fill = Value)) +
  geom_tile() +
  # scale_fill_gradientn(name = "", colors = terrain.colors(3)) +
  scale_x_discrete(name = "") +
  scale_y_discrete(name = "")

#Co-Occurence Matrix

melted_data <- melt(strain_group) #what does melt do?
w <- dcast(melted_data , genes_aff~Strains) # Make data go wider by gene's affected
x <- as.matrix(w[,-1]) #Convert to a matrx
x[is.na(x)] <- 0
x <- apply(x, 2,  function(x) as.numeric(x > 0))
v_f <- x %*% t(x)
diag(v_f) <- 0  #Set Diagonal to Zero
dimnames(v_f) <- list(w[, 1], w[,1])
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
  out$chi_table <- chi_table
  out$paired_variants <- pairs
  out$plot <- plot
  out$co_plot <- v_f
  out$co_occurence <- v_plot


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
