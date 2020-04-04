# calculate complexity, richness and in peruvian


# Japan -------------------------------------------------------------------


library(dplyr)
setwd("~/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/japan/GI.1")

japan1 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/japan/GI.1/1_japan_GI.1_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
japan1 <- subset(japan1, subset=TOTAL_DP>=400)
japan1$Sn <- with(japan1, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
japan1 %>% summarise(
  richness = nrow(japan1),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(japan1$ALT_FREQ))



japan2 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/japan/GI.1/2_japan_GI.1_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
japan2 <- subset(japan2, subset=TOTAL_DP>=400)
japan2$Sn <- with(japan2, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
japan2 %>% summarise(
  richness = nrow(japan2),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(japan2$ALT_FREQ))



japan3 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/japan/GI.1/3_japan_GI.1_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
japan3 <- subset(japan3, subset=TOTAL_DP>=400)
japan3$Sn <- with(japan3, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
japan3 %>% summarise(
  richness = nrow(japan3),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(japan3$ALT_FREQ))


japan4 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/japan/GI.1/4_japan_GI.1_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
japan4 <- subset(japan4, subset=TOTAL_DP>=400)
japan4$Sn <- with(japan4, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
japan4 %>% summarise(
  richness = nrow(japan4),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(japan4$ALT_FREQ))
japan5 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/japan/GI.1/5_japan_GI.1_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
japan5 <- subset(japan5, subset=TOTAL_DP>=400)
japan5$Sn <- with(japan5, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
japan5 %>% summarise(
  richness = nrow(japan5),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(japan5$ALT_FREQ))



japan6 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/japan/GI.1/6_japan_GI.1_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
japan6 <- subset(japan6, subset=TOTAL_DP>=400)
japan6$Sn <- with(japan6, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
japan6 %>% summarise(
  richness = nrow(japan6),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(japan6$ALT_FREQ))


japan7 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/japan/GI.1/7_japan_GI.1_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
japan7 <- subset(japan7, subset=TOTAL_DP>=400)
japan7$Sn <- with(japan7, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
japan7 %>% summarise(
  richness = nrow(japan7),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(japan7$ALT_FREQ))

japan8 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/japan/GI.1/8_japan_GI.1_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
japan8 <- subset(japan8, subset=TOTAL_DP>=400)
japan8$Sn <- with(japan8, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
japan8 %>% summarise(
  richness = nrow(japan8),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(japan8$ALT_FREQ))

japan9 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/japan/GI.1/9_japan_GI.1_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
japan9 <- subset(japan9, subset=TOTAL_DP>=400)
japan9$Sn <- with(japan9, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
japan9 %>% summarise(
  richness = nrow(japan9),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(japan9$ALT_FREQ))

japan10 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/japan/GI.1/10_japan_GI.1_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
japan10 <- subset(japan10, subset=TOTAL_DP>=400)
japan10$Sn <- with(japan10, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
japan10 %>% summarise(
  richness = nrow(japan10),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(japan10$ALT_FREQ))


japan11 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/japan/GI.1/11_japan_GI.1_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
japan11 <- subset(japan11, subset=TOTAL_DP>=400)
japan11$Sn <- with(japan11, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
japan11 %>% summarise(
  richness = nrow(japan11),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(japan11$ALT_FREQ))


japan21 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/japan/GI.2/21_japan_GI.2_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
japan21 <- subset(japan21, subset=TOTAL_DP>=400)
japan21$Sn <- with(japan21, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
japan21 %>% summarise(
  richness = nrow(japan21),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(japan21$ALT_FREQ))


japan22 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/japan/GI.2/22_japan_GI.2_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
japan22 <- subset(japan22, subset=TOTAL_DP>=400)
japan22$Sn <- with(japan22, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
japan22 %>% summarise(
  richness = nrow(japan22),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(japan22$ALT_FREQ))

japan23 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/japan/GI.2/23_japan_GI.2_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
japan23 <- subset(japan23, subset=TOTAL_DP>=400)
japan23$Sn <- with(japan23, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
japan23 %>% summarise(
  richness = nrow(japan23),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(japan23$ALT_FREQ))



japan24 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/japan/GI.2/24_japan_GI.2_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
japan24 <- subset(japan24, subset=TOTAL_DP>=400)
japan24$Sn <- with(japan24, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
japan24 %>% summarise(
  richness = nrow(japan24),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(japan24$ALT_FREQ))


japan25 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/japan/GI.2/25_japan_GI.2_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
japan25 <- subset(japan25, subset=TOTAL_DP>=400)
japan25$Sn <- with(japan25, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
japan25 %>% summarise(
  richness = nrow(japan25),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(japan25$ALT_FREQ))




# plotting of sapovirus genetic diversity
library(readxl)
data <-  read_xlsx("C:/Users/Virology Tohoku/Dropbox/phil_data/PRIMAL_sample_conc_ms_saka_kte.xlsx",
                  sheet="diversity_metrics")
head(data)
data_peru <- data %>% filter(country=="peru")
data_japan <- data %>% filter(country=="japan")
# PERU --------------------------------------------------------------------


pr<- ggboxplot(data_peru, x = "genotype", y = "richness",
          color = "genotype", palette = "lancet",xlab = "Genotype",  ylab= "Richness",
          add = "jitter", bxp.errorbar = T)+ 
  stat_compare_means( aes(label = ..p.signif..), 
                      label.x = 1.5 ,label.y = 30 )+ font("xlab", size=15)+ font("ylab", size=15) + font("xy.text", size=13)

pd<- ggboxplot(data_peru, x = "genotype", y = "distance",
               color = "genotype", palette = "lancet",xlab = "Genotype",  ylab= "distance",
               add = "jitter", bxp.errorbar = T)+ 
  stat_compare_means( aes(label = ..p.signif..), 
                      label.x = 1.5 ,label.y = 5 )+ font("xlab", size=15)+ font("ylab", size=15) + font("xy.text", size=13)

pc = ggboxplot(data_peru, x = "genotype", y = "complexity",
                    color = "genotype", palette = "lancet",xlab = "Genotype",  ylab= "Sn",
                    add = "jitter", bxp.errorbar = T) + 
  stat_compare_means( aes(label = ..p.signif..), 
                      label.x = 1.5 ,label.y = 0.6 )+ font("xlab", size=15)+ font("ylab", size=15) + font("xy.text", size=13)
# japan data
library(ggplot2)
library(ggpubr)



jr<- ggboxplot(data_japan, x = "genotype", y = "richness",
               color = "genotype", palette = "lancet",
               xlab = "Genotype",  ylab= "Richness", add = "jitter", bxp.errorbar = T) + 
  stat_compare_means( aes(label = ..p.signif..), 
                      label.x = 1.5 ,label.y = 200 )+ font("xlab", size=15)+ font("ylab", size=15) + font("xy.text", size=13)

jd<- ggboxplot(data_japan, x = "genotype", y = "distance",
               color = "genotype", palette = "lancet",xlab = "Genotype",  ylab= "Distance",
               add = "jitter", bxp.errorbar = T)+ 
  stat_compare_means( aes(label = ..p.signif..), 
                      label.x = 1.5 ,label.y = 200 )+ font("xlab", size=15)+ font("ylab", size=15) + font("xy.text", size=13)

jc = ggboxplot(data_japan, x = "genotype", y = "complexity",
               color = "genotype", palette = "lancet",xlab = "Genotype",  ylab= "Sn",
               add = "jitter", bxp.errorbar = T) + 
  stat_compare_means( aes(label = ..p.signif..), 
                      label.x = 1.5 ,label.y = 0.25 )+ font("xlab", size=15)+ font("ylab", size=15) + font("xy.text", size=13)

library(patchwork)
(jr|jd|jc) / (pr|pd|pc) + plot_layout(guides="collect") + 
  plot_annotation()

(pr|pd|pc) / 
  (jr|jd|jc) + plot_layout(guides="collect")
setwd("C:/Users/Virology Tohoku/Dropbox/manuscript-plos-pathogens/supplements")

ggsave(filename = "sapovirus genetic diversity_ok.jpeg", 
       width = 40, height = 25, dpi = 350, units = "cm")
