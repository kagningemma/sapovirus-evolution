library(readr)
library(dplyr)
library(ggplot2)
require(gridExtra)
library(grid)

# SP031X ------------------------------------------------------------------

# SP031X
setwd("C:/Users/viro10/Desktop/called_iSNV_files/aminoacid_peru/SP031X")
peru_2 <-read_tsv("2_peru_final.tsv") %>% filter (TOTAL_DP>=400) %>%  
  dplyr::mutate(Amino_acid = ifelse(S_or_NS=="S", "Synonymous", "Non-synonymous")) 

peru_5 <-read_tsv("5_peru_final.tsv") %>% filter (TOTAL_DP>=400) %>%  
  dplyr::mutate(Amino_acid = ifelse(S_or_NS=="S", "Synonymous", "Non-synonymous")) 


# plot1
g1 = ggplot(peru_2, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=Amino_acid), size=3) + 
  ggtitle("Day 4 vs Day1 (peru_2)" ) + 
  xlab("Sapovirus genome position") + 
  ylab("iSNP frequency")+
  theme_bw()+
  scale_shape_manual(values = c(3,16))

# to create rectangle 1 of g1
g1 <- g1 + geom_rect(data=peru_2, aes(xmin = 0, xmax = 600, ymin = -Inf, ymax = +Inf), 
               fill="#636363", alpha=0.05, inherit.aes = FALSE) +
 # to create rectangle2 of g1
geom_rect(data=peru_2, aes(xmin = 1300, xmax = 2250, ymin = -Inf, ymax = +Inf), 
               fill="#636363", alpha=0.05, inherit.aes = FALSE) +

# to create rectangle3 of g1
geom_rect(data=peru_2, aes(xmin = 2500, xmax = 2540, ymin = -Inf, ymax = +Inf), 
               fill="#636363", alpha=0.05, inherit.aes = FALSE) +
# to create rectangle4 of g1
geom_rect(data=peru_2, aes(xmin = 3200, xmax = 3450, ymin = -Inf, ymax = +Inf), 
          fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  
geom_rect(data=peru_2, aes(xmin = 3950, xmax = 4300, ymin = -Inf, ymax = +Inf), 
          fill="#636363", alpha=0.05, inherit.aes = FALSE) +
geom_rect(data=peru_2, aes(xmin = 4630, xmax = 4640, ymin = -Inf, ymax = +Inf), 
          fill="#636363", alpha=0.05, inherit.aes = FALSE) +
geom_rect(data=peru_2, aes(xmin = 4900, xmax = 4950, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
geom_rect(data=peru_2, aes(xmin = 5720, xmax = 5750, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
geom_rect(data=peru_2, aes(xmin = 6100, xmax = 6160, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE)+
geom_rect(data=peru_2, aes(xmin = 6400, xmax = 6460, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
geom_rect(data=peru_2, aes(xmin = 6800, xmax = 6850, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
geom_rect(data=peru_2, aes(xmin = 7200, xmax = 7250, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE)




# plot2
g2= ggplot(peru_5, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=Amino_acid), size=3) + 
  ggtitle("Day 27 vs Day1 (peru_5)" ) + 
  theme(plot.title = element_text(size = 1)) +
  xlab("Sapovirus genome position") + 
  ylab("iSNP frequency")+
  theme_bw ()
  
  g2



# to create rectangle 1 of g2
g2 <- g2 + geom_rect(data=peru_5, aes(xmin = 0, xmax = 450, ymin = -Inf, ymax = +Inf), 
               fill="#636363", alpha=0.04, inherit.aes = FALSE) +
geom_rect(data=peru_5, aes(xmin = 1400, xmax = 2250, ymin = -Inf, ymax = +Inf), 
                 fill="#636363", alpha=0.04, inherit.aes = FALSE) +
geom_rect(data=peru_5, aes(xmin = 2500, xmax = 2560, ymin = -Inf, ymax = +Inf), 
                 fill="#636363", alpha=0.04, inherit.aes = FALSE) +
geom_rect(data=peru_5, aes(xmin = 2950, xmax = 2990, ymin = -Inf, ymax = +Inf), 
                 fill="#636363", alpha=0.04, inherit.aes = FALSE) +
geom_rect(data=peru_5, aes(xmin = 3200, xmax = 3250, ymin = -Inf, ymax = +Inf), 
                 fill="#636363", alpha=0.04, inherit.aes = FALSE) +
geom_rect(data=peru_5, aes(xmin = 3700, xmax = 3750, ymin = -Inf, ymax = +Inf), 
                 fill="#636363", alpha=0.04, inherit.aes = FALSE) +
geom_rect(data=peru_5, aes(xmin = 3900, xmax = 4150, ymin = -Inf, ymax = +Inf), 
                 fill="#636363", alpha=0.04, inherit.aes = FALSE) +
geom_rect(data=peru_5, aes(xmin = 4400, xmax = 4410, ymin = -Inf, ymax = +Inf), 
                 fill="#636363", alpha=0.04, inherit.aes = FALSE) +
geom_rect(data=peru_5, aes(xmin = 4630, xmax = 4650, ymin = -Inf, ymax = +Inf), 
                 fill="#636363", alpha=0.04, inherit.aes = FALSE) +
geom_rect(data=peru_5, aes(xmin = 4800, xmax = 4850, ymin = -Inf, ymax = +Inf), 
                 fill="#636363", alpha=0.04, inherit.aes = FALSE) +
geom_rect(data=peru_5, aes(xmin = 5500, xmax = 5600, ymin = -Inf, ymax = +Inf), 
                 fill="#636363", alpha=0.04, inherit.aes = FALSE) +
geom_rect(data=peru_5, aes(xmin = 5750, xmax = 5760, ymin = -Inf, ymax = +Inf), 
                 fill="#636363", alpha=0.04, inherit.aes = FALSE) +
geom_rect(data=peru_5, aes(xmin = 5780, xmax = 5790, ymin = -Inf, ymax = +Inf), 
                 fill="#636363", alpha=0.04, inherit.aes = FALSE) +
geom_rect(data=peru_5, aes(xmin = 6050, xmax = 6150, ymin = -Inf, ymax = +Inf), 
                 fill="#636363", alpha=0.04, inherit.aes = FALSE) +
geom_rect(data=peru_5, aes(xmin = 6300, xmax = 6350, ymin = -Inf, ymax = +Inf), 
                 fill="#636363", alpha=0.04, inherit.aes = FALSE) +
geom_rect(data=peru_5, aes(xmin = 6800, xmax = 6840, ymin = -Inf, ymax = +Inf), 
                 fill="#636363", alpha=0.04, inherit.aes = FALSE) +
geom_rect(data=peru_5, aes(xmin = 7200, xmax = 7250, ymin = -Inf, ymax = +Inf), 
                 fill="#636363", alpha=0.04, inherit.aes = FALSE) 
g2


g1g2 = grid.arrange(
  g1,
  g2,
  nrow = 2,
  top = textGrob("sapovirus GI.1-SP031X",
    gp = gpar(fontface = 6, fontsize = 20)),
  bottom = textGrob(
    "",
    gp = gpar(fontface = 6, fontsize = 20),
    hjust = 1,
    x = 1
  )
)

ggsave(filename = "results/sapovirus GI.1-SP031X.tiff",
       plot = g1g2, width = 25, height = 15, dpi = 300, units = "cm")

# SP051X ------------------------------------------------------------------



setwd("C:/Users/viro10/Desktop/called_iSNV_files/aminoacid_peru/SP051X")
peru_11 <-read_tsv("11_peru_final.tsv") %>% filter (TOTAL_DP>=400) %>%  
  dplyr::mutate(Amino_acid = ifelse(S_or_NS=="S", "Synonymous", "Non-synonymous")) 

# plot1
g3 = ggplot(peru_11, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=Amino_acid), size=3) + 
  ggtitle("Day 8 vs Day1 (peru_11)" ) + 
  xlab("Sapovirus genome position") + 
  ylab("iSNP frequency")+
  theme_bw() + 
  scale_color_manual(values = c("Synonymous" = "#00BFC4", "Non-synonymous" = "red"))
 # to plot g3
g3
g3 <- g3 + geom_rect(data=peru_11, aes(xmin = 0, xmax = 450, ymin = -Inf, ymax = +Inf), 
               fill="#636363", alpha=0.05, inherit.aes = FALSE) +
geom_rect(data=peru_11, aes(xmin = 1150, xmax = 1300, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
geom_rect(data=peru_11, aes(xmin = 2850, xmax = 2890, ymin = -Inf, ymax = +Inf), 
          fill="#636363", alpha=0.05, inherit.aes = FALSE) +
geom_rect(data=peru_11, aes(xmin = 3100, xmax = 3160, ymin = -Inf, ymax = +Inf), 
          fill="#636363", alpha=0.05, inherit.aes = FALSE) +
geom_rect(data=peru_11, aes(xmin = 3500, xmax = 3560, ymin = -Inf, ymax = +Inf), 
          fill="#636363", alpha=0.05, inherit.aes = FALSE) +
geom_rect(data=peru_11, aes(xmin = 6950, xmax = 7000, ymin = -Inf, ymax = +Inf), 
          fill="#636363", alpha=0.05, inherit.aes = FALSE) +
geom_rect(data=peru_11, aes(xmin = 7200, xmax = 7250, ymin = -Inf, ymax = +Inf), 
          fill="#636363", alpha=0.05, inherit.aes = FALSE) 

g3
ggsave(filename = "results/sapovirus GI.1-SP051X.png",
       plot = g3, width = 25, height = 15, dpi = 300, units = "cm")

# SP088X-GI.1 --------------------------------------------------------------


setwd("C:/Users/viro10/Desktop/called_iSNV_files/aminoacid_peru/SP088X")

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
  ggtitle("Day 9 vs Day1(peru_22)" ) + 
  theme_bw()+
  xlab("Sapovirus genome position") + 
  ylab("iSNP frequency")+
  scale_color_manual(values = c("Synonymous" = "#00BFC4", "Non-synonymous" = "red"))

g4 <- g4 + geom_rect(data=peru_22, aes(xmin = 0, xmax = 650, ymin = -Inf, ymax = +Inf), 
               fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_22, aes(xmin = 850, xmax = 890, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_22, aes(xmin = 1200, xmax = 1400, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_22, aes(xmin = 3510, xmax = 3530, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_22, aes(xmin = 6950, xmax = 7000, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_22, aes(xmin = 7200, xmax = 7250, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) 

# plot2
g5=ggplot(peru_25, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=Amino_acid), size=3) + 
  ggtitle("Day 15 vs Day1(peru_25)" ) + 
  xlab("Sapovirus genome position") + 
  ylab("iSNP frequency")+
  theme_bw()+
  scale_color_manual(values = c("Synonymous" = "#00BFC4", "Non-synonymous" = "red"))
  g5
g5 <- g5 + geom_rect(data=peru_25, aes(xmin = 0, xmax = 650, ymin = -Inf, ymax = +Inf), 
               fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_25, aes(xmin = 900, xmax = 950, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_25, aes(xmin = 1200, xmax = 1400, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_25, aes(xmin = 3550, xmax = 3600, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_25, aes(xmin = 6850, xmax = 7000, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) 
# plot2

g6=ggplot(peru_26, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=Amino_acid), size=3) + 
  ggtitle("Day 21 vs Day1(peru_26)" ) + 
  xlab("Sapovirus genome position") + 
  ylab("iSNP frequency")+
  theme_bw()+
  scale_shape_manual(values = c(3,16))+
  labs("Amino acid mutation")

g6 <- g6 + geom_rect(data=peru_26, aes(xmin = 0, xmax = 650, ymin = -Inf, ymax = +Inf), 
               fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_26, aes(xmin = 960, xmax = 990, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_26, aes(xmin = 1200, xmax = 1400, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_26, aes(xmin = 3550, xmax = 3600, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_26, aes(xmin = 6850, xmax = 7000, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) 

# move to the next person
g4g5g6 = grid.arrange(g4, g5, g6, nrow=3, 
             top = textGrob("sapovirus GI.1-SP088X",
                            gp = gpar(fontface = 6, fontsize = 20)))


ggsave(filename = "results/sapovirus GI.1-SP088X.tiff",
       plot = g1g2, width = 25, height = 15, dpi = 300, units = "cm")

# GI.2-SP0121X ------------------------------------------------------------

setwd("C:/Users/viro10/Desktop/called_iSNV_files/aminoacid_peru/SP0121X")
peru_33 <-read_tsv("33_peru_final.tsv") %>% filter (TOTAL_DP>=400) %>%  
  dplyr::mutate(Amino_acid = ifelse(S_or_NS=="S", "Synonymous", "Non-synonymous")) 

require(gridExtra)
# plot1
g7=ggplot(peru_33, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=Amino_acid), size=3)+ 
  ggtitle("Day 12 vs Day1(peru_33)" ) + 
  xlab("Sapovirus genome position") + 
  ylab("iSNP frequency")+
  theme_bw()


g7= g7 + geom_rect(data=peru_33, aes(xmin = 1, xmax = 40, ymin = -Inf, ymax = +Inf), 
               fill="#636363", alpha=0.01, inherit.aes = FALSE) +
  geom_rect(data=peru_33, aes(xmin = 300, xmax = 360, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.03, inherit.aes = FALSE) +
  geom_rect(data=peru_33, aes(xmin = 600, xmax = 650, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.03, inherit.aes = FALSE) +
  geom_rect(data=peru_33, aes(xmin = 1550, xmax = 1570, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.03, inherit.aes = FALSE) +
  geom_rect(data=peru_33, aes(xmin = 1800, xmax = 2200, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.03, inherit.aes = FALSE) +
  geom_rect(data=peru_33, aes(xmin = 2850, xmax = 2890, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.03, inherit.aes = FALSE) +
  geom_rect(data=peru_33, aes(xmin = 3200, xmax = 3300, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.03, inherit.aes = FALSE) +
  geom_rect(data=peru_33, aes(xmin = 4400, xmax = 4650, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.03, inherit.aes = FALSE) +
  geom_rect(data=peru_33, aes(xmin = 4900, xmax = 4950, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.03, inherit.aes = FALSE) +
  geom_rect(data=peru_33, aes(xmin = 5300, xmax = 5360, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.03, inherit.aes = FALSE) +
  geom_rect(data=peru_33, aes(xmin = 5940, xmax = 5965, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.03, inherit.aes = FALSE) +
  geom_rect(data=peru_33, aes(xmin = 6200, xmax = 6300, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.03, inherit.aes = FALSE) +
geom_rect(data=peru_33, aes(xmin = 6500, xmax = 6600, ymin = -Inf, ymax = +Inf), 
          fill="#636363", alpha=0.03, inherit.aes = FALSE) +
geom_rect(data=peru_33, aes(xmin = 6750, xmax = 6780, ymin = -Inf, ymax = +Inf), 
          fill="#636363", alpha=0.05, inherit.aes = FALSE) +
geom_rect(data=peru_33, aes(xmin = 7010, xmax = 7160, ymin = -Inf, ymax = +Inf), 
          fill="#636363", alpha=0.05, inherit.aes = FALSE) 


g7

ggsave(filename = "results/sapovirus GI.2-SP121X.png",
       plot = g1g2, width = 25, height = 15, dpi = 300, units = "cm")

# GI.2-SP223X -------------------------------------------------------------


setwd("C:/Users/viro10/Desktop/called_iSNV_files/aminoacid_peru/SP223X")


peru_15 <-read_tsv("15_peru_final.tsv") %>% filter (TOTAL_DP>=400) %>%  
  dplyr::mutate(Amino_acid = ifelse(S_or_NS=="S", "Synonymous", "Non-synonymous")) 
peru_16 <-read_tsv("16_peru_final.tsv") %>% filter (TOTAL_DP>=400) %>%  
  dplyr::mutate(Amino_acid = ifelse(S_or_NS=="S", "Synonymous", "Non-synonymous")) 
peru_17 <-read_tsv("17_peru_final.tsv") %>% filter (TOTAL_DP>=400) %>%  
  dplyr::mutate(Amino_acid = ifelse(S_or_NS=="S", "Synonymous", "Non-synonymous")) 
peru_18 <-read_tsv("18_peru_final.tsv") %>% filter (TOTAL_DP>=400) %>%  
  dplyr::mutate(Amino_acid = ifelse(S_or_NS=="S", "Synonymous", "Non-synonymous")) 

# plot1
g8=ggplot(peru_15, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=Amino_acid), size=3) + 
  ggtitle("Day 5 vs Day1(peru_15)" ) + 
  xlab("Sapovirus genome position") + 
  ylab("iSNP frequency")+
  theme_bw()

g8= g8 + geom_rect(data=peru_15, aes(xmin = 0, xmax = 70, ymin = -Inf, ymax = +Inf), 
               fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_15, aes(xmin = 250, xmax = 300, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_15, aes(xmin = 550, xmax = 600, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_15, aes(xmin = 1440, xmax = 2380, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_15, aes(xmin = 2600, xmax = 2650, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_15, aes(xmin = 2900, xmax = 2950, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_15, aes(xmin = 4100, xmax = 4150, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_15, aes(xmin = 4300, xmax = 4700, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_15, aes(xmin = 5250, xmax = 5300, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_15, aes(xmin = 6200, xmax = 6500, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_15, aes(xmin = 6700, xmax = 6800, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_15, aes(xmin = 7050, xmax = 7100, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE)

g9=ggplot(peru_16, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=Amino_acid), size=3) + 
  ggtitle("Day 12 vs Day1(peru_16)" ) + 
  xlab("Sapovirus genome position") + 
  ylab("iSNP frequency")+
  theme_bw()+ 
  geom_label( 
    data= peru_16 %>% filter(stop=="stop"), # Filter data first
    aes(label=stop),
    label.padding = unit(0.05, "lines"), # Rectangle size around label
    label.size = 0.05)


g9= g9 + geom_rect(data=peru_16, aes(xmin = 0, xmax = 70, ymin = -Inf, ymax = +Inf), 
               fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_16, aes(xmin = 250, xmax = 300, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_16, aes(xmin = 550, xmax = 600, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_16, aes(xmin = 1440, xmax = 2500, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_16, aes(xmin = 2900, xmax = 2950, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_16, aes(xmin = 3200, xmax = 3300, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_16, aes(xmin = 4100, xmax = 4150, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_16, aes(xmin = 4300, xmax = 4700, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_16, aes(xmin = 5250, xmax = 5300, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_16, aes(xmin = 7050, xmax = 7100, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE)

g10=ggplot(peru_17, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=Amino_acid), size=3) + 
  ggtitle("Day 19 vs Day1(peru_17)" ) + 
  xlab("Sapovirus genome position") + 
  ylab("iSNP frequency")+
  theme_bw()

g10= g10 + geom_rect(data=peru_17, aes(xmin = 0, xmax = 70, ymin = -Inf, ymax = +Inf), 
               fill="#636363", alpha=0.05, inherit.aes = FALSE) +
geom_rect(data=peru_17, aes(xmin = 250, xmax = 300, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
geom_rect(data=peru_17, aes(xmin = 550, xmax = 600, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
geom_rect(data=peru_17, aes(xmin = 1400, xmax = 1750, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
geom_rect(data=peru_17, aes(xmin = 2050, xmax = 2100, ymin = -Inf, ymax = +Inf), 
          fill="#636363", alpha=0.05, inherit.aes = FALSE) +
geom_rect(data=peru_17, aes(xmin = 2400, xmax = 2500, ymin = -Inf, ymax = +Inf), 
          fill="#636363", alpha=0.05, inherit.aes = FALSE) +
geom_rect(data=peru_17, aes(xmin = 2530, xmax = 2570, ymin = -Inf, ymax = +Inf), 
          fill="#636363", alpha=0.05, inherit.aes = FALSE) +
geom_rect(data=peru_17, aes(xmin = 2850, xmax = 2900, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_17, aes(xmin = 3100, xmax = 3200, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_17, aes(xmin = 4080, xmax = 4130, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_17, aes(xmin = 4300, xmax = 4700, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
geom_rect(data=peru_17, aes(xmin = 6500, xmax = 6550, ymin = -Inf, ymax = +Inf), 
          fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_17, aes(xmin = 7050, xmax = 7100, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) 

 g11=ggplot(peru_18, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=Amino_acid), size=3) + 
  ggtitle("Day 26 vs Day1(peru_18)" ) + 
  xlab("Sapovirus genome position") + 
  ylab("iSNP frequency")+
  theme_bw() + 
  geom_label( 
    data= peru_18 %>% filter(stop=="stop"), # Filter data first
    aes(label=stop),
    label.padding = unit(0.05, "lines"), # Rectangle size around label
    label.size = 0.05)
  g11
  
  g11= g11 + geom_rect(data=peru_18, aes(xmin = 0, xmax = 70, ymin = -Inf, ymax = +Inf), 
                fill="#636363", alpha=0.03, inherit.aes = FALSE) +
  geom_rect(data=peru_18, aes(xmin = 270, xmax = 310, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.03, inherit.aes = FALSE) +
  geom_rect(data=peru_18, aes(xmin = 560, xmax = 600, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.03, inherit.aes = FALSE) +
  geom_rect(data=peru_18, aes(xmin = 1400, xmax = 2200, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.03, inherit.aes = FALSE) +
  geom_rect(data=peru_18, aes(xmin = 2900, xmax = 2950, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.03, inherit.aes = FALSE) +
  geom_rect(data=peru_18, aes(xmin = 3100, xmax = 3200, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.03, inherit.aes = FALSE) +
  geom_rect(data=peru_18, aes(xmin = 4070, xmax = 4100, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.03, inherit.aes = FALSE) +
  geom_rect(data=peru_18, aes(xmin = 4300, xmax = 4700, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.03, inherit.aes = FALSE) +
  geom_rect(data=peru_18, aes(xmin = 5100, xmax = 5130, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.03, inherit.aes = FALSE) +
  geom_rect(data=peru_18, aes(xmin = 7050, xmax = 7100, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.03, inherit.aes = FALSE) 

g8g9g10g11 = grid.arrange(g8, g9, g10,g11, nrow=4, 
             top = textGrob("sapovirus GI.2-SP223X", gp = gpar(fontface = 6, fontsize = 20)))

ggsave(filename = "results/sapovirus GI.2-SP223X.png",
       plot = g1g2, width = 25, height = 15, dpi = 300, units = "cm")




# GI.2 - SP265X -----------------------------------------------------------


setwd("C:/Users/viro10/Desktop/called_iSNV_files/aminoacid_peru/SP265X")


peru_21 <-read_tsv("21_peru_final.tsv") %>% filter (TOTAL_DP>=400) %>%  
  dplyr::mutate(Amino_acid = ifelse(S_or_NS=="S", "Synonymous", "Non-synonymous")) 
peru_23 <-read_tsv("23_peru_final.tsv") %>% filter (TOTAL_DP>=400) %>%  
  dplyr::mutate(Amino_acid = ifelse(S_or_NS=="S", "Synonymous", "Non-synonymous")) 
peru_27 <-read_tsv("27_peru_final.tsv") %>% filter (TOTAL_DP>=400) %>%  
  dplyr::mutate(Amino_acid = ifelse(S_or_NS=="S", "Synonymous", "Non-synonymous")) 

# plot1
g12= ggplot(peru_21, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=Amino_acid), size=3)+ 
  ggtitle("Day 5 vs Day1(peru_21)" ) + 
  xlab("Sapovirus genome position") + 
  ylab("iSNP frequency")+
  theme_bw()

g12= g12 + geom_rect(data=peru_21, aes(xmin = 0, xmax = 60, ymin = -Inf, ymax = +Inf), 
                fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_21, aes(xmin = 320, xmax = 350, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_21, aes(xmin = 600, xmax = 660, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_21, aes(xmin = 1420, xmax = 1800, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_21, aes(xmin = 2100, xmax = 2200, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_21, aes(xmin = 2400, xmax = 2450, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_21, aes(xmin = 2900, xmax = 2950, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_21, aes(xmin = 4550, xmax = 4700, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_21, aes(xmin = 5500, xmax = 5630, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_21, aes(xmin = 6500, xmax = 6800, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) 

g13=ggplot(peru_23, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=Amino_acid), size=3) + 
  ggtitle("Day 12 vs Day1(peru_23)" ) + 
  xlab("Sapovirus genome position") + 
  ylab("iSNP frequency")+
  theme_bw()+
  geom_label( 
    data= peru_23 %>% filter(stop=="stop"), # Filter data first
    aes(label=stop),
    label.padding = unit(0.05, "lines"), # Rectangle size around label
    label.size = 0.01)

g13= g13 + geom_rect(data=peru_23, aes(xmin = 0, xmax = 60, ymin = -Inf, ymax = +Inf), 
                fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_23, aes(xmin = 320, xmax = 350, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_23, aes(xmin = 600, xmax = 660, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_23, aes(xmin = 1420, xmax = 1800, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_23, aes(xmin = 2100, xmax = 2200, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_23, aes(xmin = 2300, xmax = 2500, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_23, aes(xmin = 2900, xmax = 2950, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_23, aes(xmin = 4550, xmax = 4650, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_23, aes(xmin = 5300, xmax = 5400, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_23, aes(xmin = 5500, xmax = 5550, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_23, aes(xmin = 6510, xmax = 6520, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_23, aes(xmin = 7150, xmax = 7200, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) 

g14=ggplot(peru_27, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=Amino_acid), size=3) + 
  ggtitle("Day 20 vs Day1(peru_27)" ) + 
  xlab("Sapovirus genome position") + 
  ylab("iSNP frequency")+
  theme_bw()

g14= g14 + geom_rect(data=peru_27, aes(xmin = 0, xmax = 60, ymin = -Inf, ymax = +Inf), 
                fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_27, aes(xmin = 320, xmax = 350, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_27, aes(xmin = 600, xmax = 660, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_27, aes(xmin = 1400, xmax = 2400, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_27, aes(xmin = 600, xmax = 660, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_23, aes(xmin = 2900, xmax = 2950, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_23, aes(xmin = 4550, xmax = 4650, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
  geom_rect(data=peru_23, aes(xmin = 6590, xmax = 6700, ymin = -Inf, ymax = +Inf), 
            fill="#636363", alpha=0.05, inherit.aes = FALSE) 

g12g13g14 <- grid.arrange(g12, g13, g14, nrow=3, 
             top = textGrob("sapovirus GI.2-SP265X", gp = gpar(fontface = 6, fontsize = 20)))

ggsave(filename = "results/sapovirus GI.2-SP265X.png",
       plot = g1g2, width = 25, height = 15, dpi = 300, units = "cm")

# all GI.2 samples -------------------------------------------------------------



g8g9g10g11g12g13g14 <- grid.arrange(g8,g9,g10,g11,g12, g13, g14, ncols=2, 
                          top = textGrob("sapovirus GI.2", gp = gpar(fontface = 6, fontsize = 20)))

ggsave(filename = "results/sapovirus GI.2-SP265X.png",
       plot = g1g2, width = 25, height = 15, dpi = 300, units = "cm")



library("cowplot")
plot_grid(g1, g2, g3 + rremove("x.text"), 
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)






grid.arrange(g1, g2, g3,g4, g5, g6,
             g7, g8, g9,
             g10, g11, g12, g13, g14, nrow=14)


plist <- list(g1, g2, g3,g4, g5, g6,
              g7, g8, g9,
              g10, g11, g12, g13, g14)
cowplot::plot_grid(plotlist = plist, ncol = 3)


library(reshape2)
plotData <- melt()




