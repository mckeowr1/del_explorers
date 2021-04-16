library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(stringr)
library(purrr)
library(hablar)


flat_file <- data.table::fread("WI.20210121.strain-annotation.bcsq.20210401.tsv")
test_data <- head(flat_file, n = 1000)



impact_numeric <- function(x) {
  recode(x,
         "missense" =3,
         "synonymous"=1,
         "stop_lost"= 4,
         "stop_gained"=4,
         "inframe_deletion"=3,
         "inframe_insertion"=3,
         "frameshift"=4,
         "splice_acceptor"=4,
         "splice_donor"=4,
         "start_lost"=4,
         "splice_region"=3,
         "stop_retained"=1,
         "5_prime_utr"=2,
         "3_prime_utr"=2,
         "non_coding"=2,
         "intron"=2,
         "intergenic"=2,
         "inframe_altering"=3,
         "coding_sequence"=2,
         "feature_elongation"=2,
         "start_retained"=1 ,
         "*missense" =3,
         "*synonymous"=1,
         "*stop_lost"= 4,
         "*stop_gained"=4,
         "*inframe_deletion"=3,
         "*inframe_insertion"=3,
         "*frameshift"=4,
         "*splice_acceptor"=4,
         "*splice_donor"=4,
         "*start_lost"=4,
         "*splice_region"=3,
         "*stop_retained"=1,
         "*5_prime_utr"=2,
         "*3_prime_utr"=2,
         "*non_coding"=2,
         "*intron"=2,
         "*intergenic"=2,
         "*inframe_altering"=3,
         "*coding_sequence"=2,
         "*feature_elongation"=2,
         "*start_retained"=1)

}  #convert consequence to numeric

#Collapse By Transcripts
col_tran <- flat_file %>% group_by(CHROM, POS, REF, ALT, CONSEQUENCE, WORMBASE_ID,
                                   BIOTYPE, STRAND, AMINO_ACID_CHANGE, DNA_CHANGE, Strains,
                                   BLOSUM, Grantham, Percent_Protein, GENE, VARIANT_IMPACT, DIVERGENT) %>% nest()


col_tran2 <- col_tran %>% as.data.frame()

multi_con_select <- function(df){
  check_multi <- df %>%
    dplyr::mutate(multi_con = stringr::str_detect(.$CONSEQUENCE, pattern = fixed("&")))

  single_con <- check_multi %>% filter(multi_con == FALSE) %>% as_tibble()

  multi_con <- check_multi %>% filter(multi_con == TRUE) %>% as_tibble()

  sep_con <- multi_con %>%
    tidyr::separate(CONSEQUENCE, sep = "&", into = c("CON1", "CON2", "CON3"), remove= FALSE) %>%
    dplyr::mutate(CON1n = sapply(.$CON1, impact_numeric)) %>% #Recode impact to numeric 4 - highest 1- lowest
    dplyr::mutate(CON2n = sapply(.$CON2, impact_numeric)) %>% #Recode impact to numeric 4 - highest 1- lowest
    dplyr::mutate(CON3n = sapply(.$CON3, impact_numeric)) %>% #Recode impact to numeric 4 - highest 1- lowest
    dplyr::mutate(max_con = pmax(.$CON1n, .$CON2n, .$CON3n, na.rm = TRUE)) #Print the highest consequence

  select_high <- sep_con %>%
    mutate(CONSEQUENCE = (CON1n == max_con,.$CON1,
                                ifelse(CON2n == max_con, .$CON2,
                                       ifelse(CON3n == max_con, .$CON3,
                                              NA ))))
  ref_multi_con <- select_high %>%
    select (-c("CON1", "CON2", "CON3", "CON1n", 'CON2n', "CON3n", "multi_con"))

  ref_full <- bind_rows(single_con, ref_multi_con) %>% select(-c("max_con", "multi_con"))
  return(ref_full)
}

test_f1 <- multi_con_select(col_tran2) %>% mutate(CHROM.POS = str_c(.$CHROM, ":", .$POS))


select_high_con <- function(df){
  #Find_rows that occur more than once
  chrom_pos_freq <- data.table(table(df$CHROM.POS)) %>% rename(CHROM_POS = V1) %>% filter(N > 1 )
  two_row_names <- c(chrom_pos_freq$CHROM_POS)

  single_var <- df %>% filter(!CHROM.POS %in% two_row_names)

  two_or_more <- df %>% filter(CHROM.POS %in% two_row_names) %>%
    mutate(CON_VAL = sapply(.$CONSEQUENCE, impact_numeric)) %>%
    hablar::convert(num(CON_VAL)) %>% as_tibble()

  ply <- two_or_more %>% group_by(CHROM.POS) %>%
    top_n(1, abs(CON_VAL)) %>%
    group_by(CHROM.POS) %>%
    sample_n(1) %>%
    select(-CON_VAL)

  processed <- rbind(single_var, ply)

  return(processed)
}


high_selected <- select_high_con(test_f1)

graphable <- high_selected %>% select(-data)


data.table::fwrite(graphable, "pre_processed.tsv")


#####Graphing##############



test_plot <- graphable %>% mutate(bin = round(POS/1E6)) %>% filter(CHROM != "MtDNA") %>%
  filter(!str_detect(.$CONSEQUENCE, "@")) %>%  #Filter out linker variants
  filter(!str_detect(.$CONSEQUENCE, "\\*"))
#Filter out variants after a stop codon


impact <- c("frameshift",
            "inframe_deletion",
            "inframe_insertion",
            "missense",
            "splice_acceptor",
            "splice_donor",
            "start_lost",
            "stop_gained",
            "stop_lost")

test_plot_1 <- test_plot %>%
  mutate(IMPACT = ifelse(.$CONSEQUENCE %in% impact, TRUE, FALSE))

  ggplot(test_plot_1, aes(x = bin) ) +
    geom_bar(aes(fill = IMPACT)) +
    theme_bw() +
    scale_fill_manual(values = c("#707074", "#000000")) +
    coord_cartesian(ylim = c(0 , 75000))+
    facet_grid(. ~ CHROM, scales = "free_x", space = "free_x") +
    theme(legend.position = "none",
          axis.title.x = element_text(face = "bold", size = 20),
          axis.title.y = element_text(face = "bold", size = 20),
          strip.text.x = element_text(face = "bold", size = 20),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 18)

                                                                  )+
    labs(x = "Genomic Position", y = "Number of Variants")

#Export with with = 7 Size = 3.5
