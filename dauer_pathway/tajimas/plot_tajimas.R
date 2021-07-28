library(tidyverse)
library(ggplot2)
library(cowplot)

chr3 <- data.table::fread("/Users/ryanmckeown/Documents/andersen_lab/data_exploration/tajima/chromIII_td.csv")
chr3_td <- chr3 %>% select(chrom, start, end, td)
chr3_theta <- chr3 %>% select(chrom, start, end, theta)
chr3_pi <- chr3 %>% select(chrom, start, end, pi)

chrX <- data.table::fread("/Users/ryanmckeown/Documents/andersen_lab/data_exploration/tajima/chromX_td.csv")
chrX_td <- chrX %>% select(chrom, start, end, td)
chrX_theta <- chrX %>% select(chrom, start, end, theta)
chrX_pi <- chrX %>% select(chrom, start, end, pi)


chrX_td %>%
  ggplot()+
  aes(x=((start+end)/2), y=td)+
  geom_point()+
  ylab("Tajima's D")+
  xlab("Genomic position (Mb)") + 
  geom_vline(xintercept = 14147700, colour="red")+
  geom_vline(xintercept = 14152070, colour = "red")+ 
  coord_cartesian(xlim = c(13800000, 14400000 ))





chr3_td %>%
  ggplot()+
  aes(x=((start+end)/2), y=td)+
  geom_point()+
  ylab("Tajima's D")+
  xlab("Genomic position (Mb)") + 
  geom_vline(xintercept = 2994514, colour="red")+
  geom_vline(xintercept = 3040846, colour = "red")+ 
  coord_cartesian(xlim = c(2300000, 3700000 ))

chr3_theta %>%
  ggplot()+
  aes(x=((start+end)/2), y=theta)+
  geom_point()+
  ylab("Theta")+
  xlab("Genomic position (Mb)") + 
  geom_vline(xintercept = 2994514, colour="red")+
  geom_vline(xintercept = 3040846, colour = "red")+ 
  coord_cartesian(xlim = c(2993500, 3050000 ))

chr3_pi %>%
  ggplot()+
  aes(x=((start+end)/2), y=pi)+
  geom_point()+
  ylab("Pi")+
  xlab("Genomic position (Mb)") + 
  geom_vline(xintercept = 2994514, colour="red")+
  geom_vline(xintercept = 3040846, colour = "red")+ 
  coord_cartesian(xlim = c(2993500, 3050000 ))
  


  # geom_vline(xintercept = 16069238/1e6, colour="red")+
  # geom_vline(xintercept = 16071318/1e6, colour = "red")+
  # geom_rect(aes(xmin=16069238/1e6, xmax =16071318/1e6, ymax = Inf, ymin = -Inf), color = "red",fill="red", alpha = 0.01)+
  # geom_vline(xintercept = 16071325/1e6, colour="blue", linetype = "dashed")+
  # geom_vline(xintercept = 16072738/1e6, colour="blue")+
  # geom_rect(aes(xmin=16071325/1e6, xmax =16072738/1e6, ymax = Inf, ymin = -Inf), colour = "blue",fill="blue", alpha = 0.01)+
  # cowplot::theme_cowplot(12)+
  # ylim(-2.2,0.2)+
  # theme(legend.position = "none")