library(genetics)
library(tidyverse)

gm <- read.table("del_explorers/dauer_pathway/LD/Genotype_Matrix_t.tsv", header = T)

if ( nrow(gm) > 1 ) {

  gm <- data.frame(snp_id = paste(gm$CHROM, gm$POS,
                                       sep = "_"), data.frame(gm)[, 5:ncol(gm)])

  sn <- list()

  for (i in 1:nrow(gm)) {
    sn[[i]] <- genetics::genotype(as.character(gsub(1, "T/T",
                                                    gsub(-1, "A/A", gm[i, 4:ncol(gm)]))))
  }

  test <- data.frame(sn)
  colnames(test) <- (gm$snp_id)
  ldcalc <- t(genetics::LD(test)[[4]])^2
  diag(ldcalc) <- 1

  write.table(ldcalc, "LD_between_paired_alleles.tsv", quote=F, row.names = T, col.names = NA, sep="\t")
}
