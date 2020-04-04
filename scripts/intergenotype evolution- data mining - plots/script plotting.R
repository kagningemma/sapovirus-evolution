

library(readr)
library(dplyr)
library(ggplot2)


# SP031X ------------------------------------------------------------------


# SP031X
setwd("~/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/aminoacid_peru/SP031X")
peru_2 <-read_tsv("2_peru_final.tsv") %>% filter (TOTAL_DP>=400) %>%  
  dplyr::mutate(Amino_acid = ifelse(S_or_NS=="S", "Synonymous", "Non-synonymous")) 

peru_5 <-read_tsv("5_peru_final.tsv") %>% filter (TOTAL_DP>=400) %>%  
  dplyr::mutate(Amino_acid = ifelse(S_or_NS=="S", "Synonymous", "Non-synonymous")) 

require(gridExtra)
 # plot1
g1=ggplot(peru_2, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=Amino_acid), size=3) + 
  xlim(1, 7400) + 
  ggtitle("iSNPs at day 4 compared to Day1 in SP031X-SaV GI.1" ) + 
  xlab("Position in sapovirus genome") + 
  ylab("iSNP frequency")+
  scale_shape_manual(values = c(3,16))+
  labs("Amino acid mutation")

# plot2
g2=ggplot(peru_5, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=Amino_acid), size=3) + 
  xlim(1, 7400) + 
  ggtitle("iSNPs at day 27 compared to Day1 in SP031X-SaV GI.1" ) + 
  xlab("Position in sapovirus genome") + 
  ylab("iSNP frequency")+
  scale_shape_manual(values = c(3,16))+
  labs("Amino acid mutation")
grid.arrange(g1, g2, nrow=2)

tiff(filename= "S_NS_iSNPs in SP031X.tiff", width = 1024, height = 768, res=300)


# SP051X ------------------------------------------------------------------



setwd("~/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/aminoacid_peru/SP051X")
peru_11 <-read_tsv("11_peru_final.tsv") %>% filter (TOTAL_DP>=400) %>%  
  dplyr::mutate(Amino_acid = ifelse(S_or_NS=="S", "Synonymous", "Non-synonymous")) 

# plot1
g3= ggplot(peru_11, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=Amino_acid), size=3) + 
  ggtitle("iSNPs at day 8 compared to Day1 in SP051X-savGI.1" ) + 
  xlab("Position in sapovirus genome") + 
  ylab("iSNP frequency")+
  theme_bw()+
  scale_shape_manual(values = c(3,16))+
  labs("Amino acid mutation")
g3

t2.rect1 <- data.frame(xmin=0, xmax=800, ymin=-Inf, ymax=+Inf)

g3 + 
  geom_rect(data=t2.rect1, aes(xmin = 0, xmax = 450, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.1, inherit.aes = FALSE)

tiff(filename= "S_NS_iSNPs in SP051X.tiff", width = 1024, height = 768, res=300)

RColorBrewer::brewer.pal()

# SP088X-GI.1 --------------------------------------------------------------


setwd("~/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/aminoacid_peru/SP088X")

peru_22 <-read_tsv("22_peru_final.tsv") %>% filter (TOTAL_DP>=400) %>%  
  dplyr::mutate(Amino_acid = ifelse(S_or_NS=="S", "Synonymous", "Non-synonymous")) 

peru_25 <-read_tsv("25_peru_final.tsv") %>% filter (TOTAL_DP>=400) %>%  
  dplyr::mutate(Amino_acid = ifelse(S_or_NS=="S", "Synonymous", "Non-synonymous")) 

peru_26 <-read_tsv("26_peru_final.tsv") %>% filter (TOTAL_DP>=400) %>%  
  dplyr::mutate(Amino_acid = ifelse(S_or_NS=="S", "Synonymous", "Non-synonymous")) 

require(gridExtra)
# plot1
g4=ggplot(peru_22, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=Amino_acid), size=3) + 
  xlim(1, 7400) + 
  ggtitle("iSNPs at day 9 compared to Day1 in SP088X-SaV GI.1" ) + 
  xlab("Position in sapovirus genome") + 
  ylab("iSNP frequency")+
  scale_shape_manual(values = c(3,16))+
  labs("Amino acid mutation")

# plot2
g5=ggplot(peru_25, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=Amino_acid), size=3) + 
  xlim(1, 7400) + 
  ggtitle("iSNPs at day 15 compared to Day1 in SP088X-SaV GI.1" ) + 
  xlab("Position in sapovirus genome") + 
  ylab("iSNP frequency")+
  scale_shape_manual(values = c(3,16))+
  labs("Amino acid mutation")

# plot2
g6=ggplot(peru_26, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=Amino_acid), size=3) + 
  xlim(1, 7400) + 
  ggtitle("iSNPs at day 21 compared to Day1 in SP088X-SaV GI.1" ) + 
  xlab("Position in sapovirus genome") + 
  ylab("iSNP frequency")+
  scale_shape_manual(values = c(3,16))+
  labs("Amino acid mutation")
grid.arrange(g4, g5, g6, nrow=3)

tiff(filename= "S_SN_iSNPs in SP088X.tiff", width = 1024, height = 768, res=300)



# GI.2-SP0121X ------------------------------------------------------------

setwd("~/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/aminoacid_peru/SP0121X")
peru_33 <-read_tsv("33_peru_final.tsv") %>% filter (TOTAL_DP>=400) %>%  
  dplyr::mutate(Amino_acid = ifelse(S_or_NS=="S", "Synonymous", "Non-synonymous")) 

require(gridExtra)
# plot1
g7=ggplot(peru_33, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=Amino_acid), size=3) + 
  xlim(1, 7400) + 
  ggtitle("iSNPs at day 9 compared to Day1 in SP0121X-SaV GI.2" ) + 
  xlab("Position in sapovirus genome") + 
  ylab("iSNP frequency")+
  scale_shape_manual(values = c(3,16))+
  labs("Amino acid mutation")
g7
tiff(filename= "GI.2.S_SN_iSNPs in SP0121X.tiff", width = 1024, height = 768, res=300)


# GI.2-SP223X -------------------------------------------------------------


setwd("~/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/aminoacid_peru/SP223X")


peru_15 <-read_tsv("15_peru_final.tsv") %>% filter (TOTAL_DP>=400) %>%  
  dplyr::mutate(Amino_acid = ifelse(S_or_NS=="S", "Synonymous", "Non-synonymous")) 
peru_16 <-read_tsv("16_peru_final.tsv") %>% filter (TOTAL_DP>=400) %>%  
  dplyr::mutate(Amino_acid = ifelse(S_or_NS=="S", "Synonymous", "Non-synonymous")) 
peru_17 <-read_tsv("17_peru_final.tsv") %>% filter (TOTAL_DP>=400) %>%  
  dplyr::mutate(Amino_acid = ifelse(S_or_NS=="S", "Synonymous", "Non-synonymous")) 
peru_18 <-read_tsv("18_peru_final.tsv") %>% filter (TOTAL_DP>=400) %>%  
  dplyr::mutate(Amino_acid = ifelse(S_or_NS=="S", "Synonymous", "Non-synonymous")) 
require(gridExtra)
# plot1
g8=ggplot(peru_15, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=Amino_acid), size=3) + 
  xlim(1, 7400) + 
  ggtitle("iSNPs at day 9 compared to Day1 in SP0121X-SaV GI.2" ) + 
  xlab("Position in sapovirus genome") + 
  ylab("iSNP frequency")+
  scale_shape_manual(values = c(3,16))+
  labs("Amino acid mutation")
g9=ggplot(peru_16, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=Amino_acid), size=3) + 
  xlim(1, 7400) + 
  ggtitle("iSNPs at day 9 compared to Day1 in SP0121X-SaV GI.2" ) + 
  xlab("Position in sapovirus genome") + 
  ylab("iSNP frequency")+
  scale_shape_manual(values = c(3,16))+
  labs("Amino acid mutation")

g10=ggplot(peru_17, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=Amino_acid), size=3) + 
  xlim(1, 7400) + 
  ggtitle("iSNPs at day 9 compared to Day1 in SP0121X-SaV GI.2" ) + 
  xlab("Position in sapovirus genome") + 
  ylab("iSNP frequency")+
  scale_shape_manual(values = c(3,16))+
  labs("Amino acid mutation")

g11=ggplot(peru_18, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=Amino_acid), size=3) + 
  xlim(1, 7400) + 
  ggtitle("iSNPs at day 9 compared to Day1 in SP0121X-SaV GI.2" ) + 
  xlab("Position in sapovirus genome") + 
  ylab("iSNP frequency")+
  scale_shape_manual(values = c(3,16))+
  labs("Amino acid mutation")

grid.arrange(g8, g9, g10,g11, nrow=4)

tiff(filename= "GI.2.S_SN_iSNPs in SP0223X.tiff", width = 1024, height = 768, res=300)


# GI.2 - SP265X -----------------------------------------------------------


setwd("~/MyRdirectory/SAPO/Sapovirus evolution over time/called iSNV files/aminoacid_peru/SP265X")


peru_21 <-read_tsv("21_peru_final.tsv") %>% filter (TOTAL_DP>=400) %>%  
  dplyr::mutate(Amino_acid = ifelse(S_or_NS=="S", "Synonymous", "Non-synonymous")) 
peru_23 <-read_tsv("23_peru_final.tsv") %>% filter (TOTAL_DP>=400) %>%  
  dplyr::mutate(Amino_acid = ifelse(S_or_NS=="S", "Synonymous", "Non-synonymous")) 
peru_27 <-read_tsv("27_peru_final.tsv") %>% filter (TOTAL_DP>=400) %>%  
  dplyr::mutate(Amino_acid = ifelse(S_or_NS=="S", "Synonymous", "Non-synonymous")) 

# plot1
g12=ggplot(peru_21, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=Amino_acid), size=3) + 
  xlim(1, 7400) + 
  ggtitle("iSNPs at day 5 compared to Day1 in SP265X-SaV GI.2" ) + 
  xlab("Position in sapovirus genome") + 
  ylab("iSNP frequency")+
  scale_shape_manual(values = c(3,16))+
  labs("Amino acid mutation")
g13=ggplot(peru_23, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=Amino_acid), size=3) + 
  xlim(1, 7400) + 
  ggtitle("iSNPs at day 12 compared to Day1 in SP265X-SaV GI.2" ) + 
  xlab("Position in sapovirus genome") + 
  ylab("iSNP frequency")+
  scale_shape_manual(values = c(3,16))+
  labs("Amino acid mutation")

g14=ggplot(peru_27, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=Amino_acid), size=3) + 
  xlim(1, 7400) + 
  ggtitle("iSNPs at day 23 compared to Day1 in SP265X-SaV GI.2" ) + 
  xlab("Position in sapovirus genome") + 
  ylab("iSNP frequency")+
  scale_shape_manual(values = c(3,16))+
  labs("Amino acid mutation")



grid.arrange(g12, g13, g14, nrow=3)

tiff(filename= "GI.2.S_SN_iSNPs in SP0265X.tiff", width = 1024, height = 768, res=300)




# all samples -------------------------------------------------------------


grid.arrange(g1, g2, g3,g4, g5, g6,
             g7, g8, g9,
             g10, g11, g12, g13, g14, nrow=14)

