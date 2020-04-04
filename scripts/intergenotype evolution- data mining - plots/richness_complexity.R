# calculate complexity, richness and in peruvian


# peru -------------------------------------------------------------------


library(dplyr)


peru1 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/peru/1_peru_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
peru1 <- subset(peru1, subset=TOTAL_DP>=400)
peru1$Sn <- with(peru1, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
peru1 %>% summarise(
  richness = nrow(peru1),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(peru1$ALT_FREQ))



peru2 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/peru/2_peru_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
peru2 <- subset(peru2, subset=TOTAL_DP>=400)
peru2$Sn <- with(peru2, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
peru2 %>% summarise(
  richness = nrow(peru2),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(peru2$ALT_FREQ))



peru3 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/peru/3_peru_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
peru3 <- subset(peru3, subset=TOTAL_DP>=400)
peru3$Sn <- with(peru3, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
peru3 %>% summarise(
  richness = nrow(peru3),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(peru3$ALT_FREQ))


peru4 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/peru/4_peru_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
peru4 <- subset(peru4, subset=TOTAL_DP>=400)
peru4$Sn <- with(peru4, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
peru4 %>% summarise(
  richness = nrow(peru4),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(peru4$ALT_FREQ))

peru5 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/peru/5_peru_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
peru5 <- subset(peru5, subset=TOTAL_DP>=400)
peru5$Sn <- with(peru5, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
peru5 %>% summarise(
  richness = nrow(peru5),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(peru5$ALT_FREQ))



peru9 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/peru/9_peru_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
peru9 <- subset(peru9, subset=TOTAL_DP>=400)
peru9$Sn <- with(peru9, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
peru9 %>% summarise(
  richness = nrow(peru9),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(peru9$ALT_FREQ))


peru11 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/peru/11_peru_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
peru11 <- subset(peru11, subset=TOTAL_DP>=400)
peru11$Sn <- with(peru11, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
peru11 %>% summarise(
  richness = nrow(peru11),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(peru11$ALT_FREQ))


peru14 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/peru/14_peru_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
peru14 <- subset(peru14, subset=TOTAL_DP>=400)
peru14$Sn <- with(peru14, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
peru14 %>% summarise(
  richness = nrow(peru14),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(peru14$ALT_FREQ))

peru15 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/peru/15_peru_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
peru15 <- subset(peru15, subset=TOTAL_DP>=400)
peru15$Sn <- with(peru15, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
peru15 %>% summarise(
  richness = nrow(peru15),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(peru15$ALT_FREQ))


peru16 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/peru/16_peru_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
peru16 <- subset(peru16, subset=TOTAL_DP>=400)
peru16$Sn <- with(peru16, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
peru16 %>% summarise(
  richness = nrow(peru16),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(peru16$ALT_FREQ))


peru17 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/peru/17_peru_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
peru17 <- subset(peru17, subset=TOTAL_DP>=400)
peru17$Sn <- with(peru17, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
peru17 %>% summarise(
  richness = nrow(peru17),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(peru17$ALT_FREQ))


peru18 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/peru/18_peru_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
peru18 <- subset(peru18, subset=TOTAL_DP>=400)
peru18$Sn <- with(peru18, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
peru18 %>% summarise(
  richness = nrow(peru18),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(peru18$ALT_FREQ))

peru19 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/peru/19_peru_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
peru19 <- subset(peru19, subset=TOTAL_DP>=400)
peru19$Sn <- with(peru19, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
peru19 %>% summarise(
  richness = nrow(peru19),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(peru19$ALT_FREQ))



peru20 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/peru/20_peru_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
peru20 <- subset(peru20, subset=TOTAL_DP>=400)
peru20$Sn <- with(peru20, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
peru20 %>% summarise(
  richness = nrow(peru20),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(peru20$ALT_FREQ))


peru21 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/peru/21_peru_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
peru21 <- subset(peru21, subset=TOTAL_DP>=400)
peru21$Sn <- with(peru21, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
peru21 %>% summarise(
  richness = nrow(peru21),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(peru21$ALT_FREQ))

peru22 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/peru/22_peru_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
peru22 <- subset(peru22, subset=TOTAL_DP>=400)
peru22$Sn <- with(peru22, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
peru22 %>% summarise(
  richness = nrow(peru22),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(peru22$ALT_FREQ))



peru23 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/peru/23_peru_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
peru23 <- subset(peru23, subset=TOTAL_DP>=400)
peru23$Sn <- with(peru23, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
peru23 %>% summarise(
  richness = nrow(peru23),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(peru23$ALT_FREQ))



peru25 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/peru/25_peru_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
peru25 <- subset(peru25, subset=TOTAL_DP>=400)
peru25$Sn <- with(peru25, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
peru25 %>% summarise(
  richness = nrow(peru25),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(peru25$ALT_FREQ))


peru26 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/peru/26_peru_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
peru26 <- subset(peru26, subset=TOTAL_DP>=400)
peru26$Sn <- with(peru26, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
peru26 %>% summarise(
  richness = nrow(peru26),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(peru26$ALT_FREQ))


peru27 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/peru/27_peru_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
peru27 <- subset(peru27, subset=TOTAL_DP>=400)
peru27$Sn <- with(peru27, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
peru27 %>% summarise(
  richness = nrow(peru27),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(peru27$ALT_FREQ))



peru28 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/peru/28_peru_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
peru28 <- subset(peru28, subset=TOTAL_DP>=400)
peru28$Sn <- with(peru28, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
peru28 %>% summarise(
  richness = nrow(peru28),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(peru28$ALT_FREQ))

peru30 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/peru/30_peru_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
peru30 <- subset(peru30, subset=TOTAL_DP>=400)
peru30$Sn <- with(peru30, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
peru30 %>% summarise(
  richness = nrow(peru30),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(peru30$ALT_FREQ))



peru31 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/peru/31_peru_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
peru31 <- subset(peru31, subset=TOTAL_DP>=400)
peru31$Sn <- with(peru31, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
peru31 %>% summarise(
  richness = nrow(peru31),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(peru31$ALT_FREQ))



peru33 <- 
  read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/peru/33_peru_final.tsv",
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
peru33 <- subset(peru33, subset=TOTAL_DP>=400)
peru33$Sn <- with(peru33, -((ALT_FREQ*log(ALT_FREQ))+(1-(ALT_FREQ))*log(1-(ALT_FREQ)))/log(2))
peru33 %>% summarise(
  richness = nrow(peru33),
  complexity = round(mean(Sn, na.rm=TRUE), 3), 
  distance = sum(peru33$ALT_FREQ))



# plotting of sapovirus genetic diversity
library(readxl)
data <-  read_xlsx("C:/Users/Virology Tohoku/Dropbox/phil_data/PRIMAL_sample_conc_ms_saka_kte.xlsx",
                   sheet="diversity_metrics")
head(data)
data_peru <- data %>% filter(country=="peru")
data_japan <- data %>% filter(country=="japan")


# plots -------------------------------------------------------------------


## plots


pr<- ggboxplot(data_peru, x = "genotype", y = "richness",
               color = "genotype", palette = "lancet",xlab = "Genotype",  ylab= "Richness",
               add = "jitter", bxp.errorbar = T)+ 
  stat_compare_means( aes(label = ..p.signif..), 
                      label.x = 1.5 ,label.y = 35 )+ font("xlab", size=15)+ font("ylab", size=15) + font("xy.text", size=13)

pd<- ggboxplot(data_peru, x = "genotype", y = "distance",
               color = "genotype", palette = "lancet",xlab = "Genotype",  ylab= "distance",
               add = "jitter", bxp.errorbar = T)+ 
  stat_compare_means( aes(label = ..p.signif..), 
                      label.x = 1.3 ,label.y = 5 )+ font("xlab", size=15)+ font("ylab", size=15) + font("xy.text", size=13)

pc = ggboxplot(data_peru, x = "genotype", y = "complexity",
               color = "genotype", palette = "lancet",xlab = "Genotype",  ylab= "Sn",
               add = "jitter", bxp.errorbar = T) + 
  stat_compare_means( aes(label = ..p.signif..), 
                      label.x = 1.5 ,label.y = 0.7 )+ font("xlab", size=15)+ font("ylab", size=15) + font("xy.text", size=13)
# japan data


library(patchwork)

pr+pd+pc + plot_layout(guides="collect")

setwd("C:/Users/Virology Tohoku/Dropbox/manuscript-plos-pathogens/supplements")

ggsave(filename = "sapovirus_peru_genetic diversity_ok.tiff", 
       width = 22, height = 15, dpi = 350, units = "cm")

