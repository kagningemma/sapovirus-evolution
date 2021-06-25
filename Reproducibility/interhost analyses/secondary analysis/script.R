
# following the secondary mapping 
# to generate the count of snps and NS mutation
library(Rcmdr)
library(dplyr)
library(readr)
library(tidyr)
files =  list.files()
setwd("C:/Users/Virology Tohoku/Google Drive/Manuscripts/My own/final edits/tsv files_secondary_20210622")


# 1-import data -----------------------------------------------------------


jap1 <- read_tsv("1_japan_GI.1_final.tsv") %>% 
  filter(TOTAL_DP>= 400 & ALT_FREQ >=0.03 & !is.na(ALT_CODON) & !is.na(ALT_CODON)) %>% 
  select(-c(REGION,REF_DP, ALT_QUAL, ALT_DP, REF_RV, ALT_RV , ALT_RV, REF_CODON, PASS) ) %>% 
  mutate(mut = if_else(REF_AA==ALT_AA, "syn", "non_syn"),
         genotype= "GI.1", outbreak= "2013_k9", sample="s63")
table(jap1$mut)

jap2 <- read_tsv("2_japan_GI.1_final.tsv") %>% 
  filter(TOTAL_DP>= 400 & ALT_FREQ >=0.03 & !is.na(ALT_CODON) ) %>% 
  select(-c(REGION,REF_DP, ALT_QUAL, ALT_DP, REF_RV, GFF_FEATURE, ALT_RV , ALT_RV, REF_CODON, PASS) ) %>% 
  mutate(mut = if_else(REF_AA==ALT_AA, "syn", "non_syn"),genotype= "GI.1", 
         outbreak= "2013_k9", sample="s64")
table(jap2$mut)

jap3 <- read_tsv("3_japan_GI.1_final.tsv") %>% 
  filter(TOTAL_DP>= 400 & ALT_FREQ >=0.03 & !is.na(ALT_CODON) ) %>% 
  select(-c(REGION,REF_DP, ALT_QUAL, ALT_DP, REF_RV, GFF_FEATURE, ALT_RV , ALT_RV, REF_CODON, PASS) ) %>% 
  mutate(mut = if_else(REF_AA==ALT_AA, "syn", "non_syn"), genotype= "GI.1", 
         outbreak= "2013_k9", sample="s65")
table(jap3$mut)

jap4 <- read_tsv("4_japan_GI.1_final.tsv") %>% 
  filter(TOTAL_DP>= 400 & ALT_FREQ >=0.03 & !is.na(ALT_CODON) ) %>% 
  select(-c(REGION,REF_DP, ALT_QUAL, ALT_DP, REF_RV, GFF_FEATURE, ALT_RV , ALT_RV, REF_CODON, PASS) ) %>% 
  mutate(mut = if_else(REF_AA==ALT_AA, "syn", "non_syn"), genotype= "GI.1", 
         outbreak= "2015_k5", sample="s90")
table(jap4$mut)

jap5 <- read_tsv("5_japan_GI.1_final.tsv") %>% 
  filter(TOTAL_DP>= 400 & ALT_FREQ >=0.03 & !is.na(ALT_CODON) ) %>% 
  select(-c(REGION,REF_DP, ALT_QUAL, ALT_DP, REF_RV, GFF_FEATURE, ALT_RV , ALT_RV, REF_CODON, PASS) ) %>% 
  mutate(mut = if_else(REF_AA==ALT_AA, "syn", "non_syn"), genotype= "GI.1", 
         outbreak= "2015_k5", sample="s91")
table(jap5$mut)


jap6 <- read_tsv("6_japan_GI.1_final.tsv") %>% 
  filter(TOTAL_DP>= 400 & ALT_FREQ >=0.03 & !is.na(ALT_CODON) ) %>% 
  select(-c(REGION,REF_DP, ALT_QUAL, ALT_DP, REF_RV, GFF_FEATURE, ALT_RV , ALT_RV, REF_CODON, PASS) ) %>% 
  mutate(mut = if_else(REF_AA==ALT_AA, "syn", "non_syn"), genotype= "GI.1", 
         outbreak= "2015_k5", sample="s92")
table(jap6$mut)


jap7 <- read_tsv("7_japan_GI.1_final.tsv") %>% 
  filter(TOTAL_DP>= 400 & ALT_FREQ >=0.03 & !is.na(ALT_CODON) ) %>% 
  select(-c(REGION,REF_DP, ALT_QUAL, ALT_DP, REF_RV, GFF_FEATURE, ALT_RV , ALT_RV, REF_CODON, PASS) ) %>% 
  mutate(mut = if_else(REF_AA==ALT_AA, "syn", "non_syn"), genotype= "GI.1",
         outbreak= "2015_k5", sample="s93")
table(jap7$mut)

jap8 <- read_tsv("8_japan_GI.1_final.tsv") %>% 
  filter(TOTAL_DP>= 400 & ALT_FREQ >=0.03 & !is.na(ALT_CODON) ) %>% 
  select(-c(REGION,REF_DP, ALT_QUAL, ALT_DP, REF_RV, GFF_FEATURE, ALT_RV , ALT_RV, REF_CODON, PASS) ) %>% 
  mutate(mut = if_else(REF_AA==ALT_AA, "syn", "non_syn"), genotype= "GI.1", 
         outbreak= "2015_k5", sample="s94")
table(jap8$mut)

#jap9 <- read_tsv("9_japan_GI.1_final.tsv") %>% 
#  filter(TOTAL_DP>= 400 & ALT_FREQ >=0.03 & !is.na(ALT_CODON) ) %>% 
 # select(-c(REGION,REF_DP, ALT_QUAL, ALT_DP, REF_RV, GFF_FEATURE, ALT_RV , ALT_RV, REF_CODON, PASS) ) %>% 
  #mutate(mut = if_else(REF_AA==ALT_AA, "syn", "non_syn"), sample="s82")
#table(jap9$mut)



jap21 <- read_tsv("21_japan_GI.2_final.tsv") %>% 
  filter(TOTAL_DP>= 400 & ALT_FREQ >=0.03 & !is.na(ALT_CODON) ) %>% 
  select(-c(REGION,REF_DP, ALT_QUAL, ALT_DP, REF_RV, GFF_FEATURE, ALT_RV , ALT_RV, REF_CODON, PASS) ) %>% 
  mutate(mut = if_else(REF_AA==ALT_AA, "syn", "non_syn"), genotype= "GI.2", 
         outbreak= "2013_k4",sample="s66")
table(jap21$mut)

jap22 <- read_tsv("22_japan_GI.2_final.tsv") %>% 
  filter(TOTAL_DP>= 400 & ALT_FREQ >=0.03 & !is.na(ALT_CODON) ) %>% 
  select(-c(REGION,REF_DP, ALT_QUAL, ALT_DP, REF_RV, GFF_FEATURE, ALT_RV , ALT_RV, REF_CODON, PASS) ) %>% 
  mutate(mut = if_else(REF_AA==ALT_AA, "syn", "non_syn"), genotype= "GI.2", 
         outbreak= "2013_k4", sample="s67")
table(jap22$mut)

jap23 <- read_tsv("23_japan_GI.2_final.tsv") %>% 
  filter(TOTAL_DP>= 400 & ALT_FREQ >=0.03 & !is.na(ALT_CODON) ) %>% 
  select(-c(REGION,REF_DP, ALT_QUAL, ALT_DP, REF_RV, GFF_FEATURE, ALT_RV , ALT_RV, REF_CODON, PASS) ) %>% 
  mutate(mut = if_else(REF_AA==ALT_AA, "syn", "non_syn"), genotype= "GI.2", 
         outbreak= "2013_k4", sample="s68")
table(jap23$mut)

jap24 <- read_tsv("24_japan_GI.2_final.tsv") %>% 
  filter(TOTAL_DP>= 400 & ALT_FREQ >=0.03 & !is.na(ALT_CODON) ) %>% 
  select(-c(REGION,REF_DP, ALT_QUAL, ALT_DP, REF_RV, GFF_FEATURE, ALT_RV , ALT_RV, REF_CODON, PASS) ) %>% 
  mutate(mut = if_else(REF_AA==ALT_AA, "syn", "non_syn"), genotype= "GI.2", 
         outbreak= "2013_k4", sample="s69")
table(jap24$mut)

jap25 <- read_tsv("25_japan_GI.2_final.tsv") %>% 
  filter(TOTAL_DP>= 400 & ALT_FREQ >=0.03 & !is.na(ALT_CODON) ) %>% 
  select(-c(REGION,REF_DP, ALT_QUAL, ALT_DP, REF_RV, GFF_FEATURE, ALT_RV , ALT_RV, REF_CODON, PASS) ) %>% 
  mutate(mut = if_else(REF_AA==ALT_AA, "syn", "non_syn"), genotype= "GI.2", 
         outbreak= "2013_k4", sample="s70")
table(jap25$mut)



# 2- make dataframe for supplements ---------------------------------------


# table of mutations compare with outbreak consensus in miyagi. 
library(plyr)
secondary_SNPS_table <- join_all(list(jap1,jap2,jap3, jap4,jap5, jap6,
              jap7, jap8,jap21, jap22,jap23, jap24,
              jap25), type='full')

# export this table for supplements file 
write.csv(secondary_SNPS_table, 
          "outbreak_secondary_SNPS_table.csv")



# 3- Multiple plots for newfigure 3 ---------------------------------------------


require(gridExtra)
require(ggplot2)
# plot1
j4=ggplot(jap4, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=mut), size=7) + 
  ggtitle("s90" ) + 
  theme_bw()+
  xlab("Sapovirus genome position") + 
    ylab("iSNP frequency")+ scale_color_manual(values = c("syn" = "#00BFC4", "non_syn" = "red"))+
  #scale_color_manual(values = c("Synonymous" = "#00BFC4", "Non-synonymous" = "red"))+
  theme(legend.position = "none",panel.background = element_rect (fill = '#E7E7E7')) + 
  labs(x=NULL, y=NULL) + expand_limits(x = c(1, 7300))+ 
  theme(plot.title=element_text(size=20, vjust=0.5, family="Times New Roman"),  
        axis.text.x = element_text(size=14, vjust=0.5, family="Times New Roman"), 
        axis.text.y = element_text(size=14, vjust=0.5, family="Times New Roman"), 
        axis.title.x=element_text(vjust=-0.25, size=40, family="Times New Roman"),
        axis.title.y=element_text(vjust=-0.25, size=40, family="Times New Roman"), 
        legend.text= element_text(size=23, family="Times New Roman"), 
        legend.title = element_text(size=23, family="Times New Roman"))

#j4 <- j4 + geom_rect(data=jap4, aes(xmin = 0, xmax = 650, ymin = -Inf, ymax = +Inf), 
#                    fill="#636363", alpha=0.05, inherit.aes = FALSE) +
#  geom_rect(data=jap4, aes(xmin = 850, xmax = 890, ymin = -Inf, ymax = +Inf), 
#            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
#  geom_rect(data=jap4, aes(xmin = 1200, xmax = 1400, ymin = -Inf, ymax = +Inf), 
#            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
#  geom_rect(data=jap4, aes(xmin = 3510, xmax = 3530, ymin = -Inf, ymax = +Inf), 
#            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
#  geom_rect(data=jap4, aes(xmin = 6950, xmax = 7000, ymin = -Inf, ymax = +Inf), 
#            fill="#636363", alpha=0.05, inherit.aes = FALSE) +
#  geom_rect(data=jap4, aes(xmin = 7200, xmax = 7250, ymin = -Inf, ymax = +Inf), 
#            fill="#636363", alpha=0.05, inherit.aes = FALSE) 

j5=ggplot(jap5, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=mut), size=7) + 
  ggtitle("s91" ) + 
  theme_bw()+
  xlab("Sapovirus genome position") + 
    ylab("iSNP frequency")+ scale_color_manual(values = c("syn" = "#00BFC4", "non_syn" = "red"))+
  #scale_color_manual(values = c("Synonymous" = "#00BFC4", "Non-synonymous" = "red"))+
  theme(legend.position = "none",panel.background = element_rect (fill = '#E7E7E7')) + 
  labs(x=NULL, y=NULL)+ expand_limits(x = c(1, 7300))+ 
  theme(plot.title=element_text(size=20, vjust=0.5, family="Times New Roman"),  
        axis.text.x = element_text(size=14, vjust=0.5, family="Times New Roman"), 
        axis.text.y = element_text(size=14, vjust=0.5, family="Times New Roman"), 
        axis.title.x=element_text(vjust=-0.25, size=40, family="Times New Roman"),
        axis.title.y=element_text(vjust=-0.25, size=40, family="Times New Roman"), 
        legend.text= element_text(size=23, family="Times New Roman"), 
        legend.title = element_text(size=23, family="Times New Roman"))

j6=ggplot(jap6, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=mut), size=7) + 
  ggtitle("s92" ) + 
  theme_bw()+ 
  scale_color_manual(values = c("syn" = "#00BFC4", "non_syn" = "red"))+
  #scale_color_manual(values = c("Synonymous" = "#00BFC4", "Non-synonymous" = "red"))+
  theme(legend.position = "none",panel.background = element_rect (fill = '#E7E7E7')) + 
  labs(x=NULL, y="iSNP frequency")+ expand_limits(x = c(1, 7300))+ 
  theme(plot.title=element_text(size=20, vjust=0.5, family="Times New Roman"),  
        axis.text.x = element_text(size=14, vjust=0.5, family="Times New Roman"), 
        axis.text.y = element_text(size=14, vjust=0.5, family="Times New Roman"), 
        axis.title.x=element_text(vjust=-0.25, size=40, family="Times New Roman"),
        axis.title.y=element_text(vjust=-0.25, size=20, family="Times New Roman"), 
        legend.text= element_text(size=23, family="Times New Roman"), 
        legend.title = element_text(size=23, family="Times New Roman"))

j7=ggplot(jap7, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=mut), size=7) + 
  ggtitle("s93" ) + 
  theme_bw()+
  xlab("Sapovirus genome position") + 
    ylab("iSNP frequency")+ scale_color_manual(values = c("syn" = "#00BFC4", "non_syn" = "red"))+
  #scale_color_manual(values = c("Synonymous" = "#00BFC4", "Non-synonymous" = "red"))+
  theme(legend.position = "right", 
        panel.background = element_rect (fill = '#E7E7E7')) + 
  labs(x=NULL, y=NULL)+ expand_limits(x = c(1, 7300))+ 
  theme(plot.title=element_text(size=20, vjust=0.5, family="Times New Roman"),  
        axis.text.x = element_text(size=14, vjust=0.5, family="Times New Roman"), 
        axis.text.y = element_text(size=14, vjust=0.5, family="Times New Roman"), 
        axis.title.x=element_text(vjust=-0.25, size=40, family="Times New Roman"),
        axis.title.y=element_text(vjust=-0.25, size=40, family="Times New Roman"), 
        legend.text= element_text(size=23, family="Times New Roman"), 
        legend.title = element_text(size=23, family="Times New Roman"))

## initial plotting figur 3

j8=ggplot(jap8, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=mut), size=7) + 
  ggtitle("s94" ) + 
  theme_bw()+ labs(x="Sapovirus genome position", y=NULL) + 
  scale_color_manual(values = c("syn" = "#00BFC4", "non_syn" = "red"))+
  #scale_color_manual(values = c("Synonymous" = "#00BFC4", "Non-synonymous" = "red"))+
  theme(legend.position = "none", 
        panel.background = element_rect (fill = '#E7E7E7'), 
        axis.title.x = element_text(size = 20))+ expand_limits(x = c(1, 7300))+ 
  theme(plot.title=element_text(size=20, vjust=0.5, family="Times New Roman"),  
        axis.text.x = element_text(size=14, vjust=0.5, family="Times New Roman"), 
        axis.text.y = element_text(size=14, vjust=0.5, family="Times New Roman"), 
        axis.title.x=element_text(vjust=-0.25, size=40, family="Times New Roman"),
        axis.title.y=element_text(vjust=-0.25, size=40, family="Times New Roman"), 
        legend.text= element_text(size=23, family="Times New Roman"), 
        legend.title = element_text(size=23, family="Times New Roman"))

j21=ggplot(jap21, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=mut), size=7) + 
  ggtitle("s66" ) + 
  theme_bw()+
  xlab("Sapovirus genome position") + 
    ylab("iSNP frequency")+ scale_color_manual(values = c("syn" = "#00BFC4", "non_syn" = "red"))+
  #scale_color_manual(values = c("Synonymous" = "#00BFC4", "Non-synonymous" = "red"))+
  theme(legend.position = "none", 
        panel.background = element_rect (fill = '#E7E7E7')) + 
  labs(x=NULL, y=NULL)+ expand_limits(x = c(1, 7300))+ 
  theme(plot.title=element_text(size=20, vjust=0.5, family="Times New Roman"),  
        axis.text.x = element_text(size=14, vjust=0.5, family="Times New Roman"), 
        axis.text.y = element_text(size=14, vjust=0.5, family="Times New Roman"), 
        axis.title.x=element_text(vjust=-0.25, size=40, family="Times New Roman"),
        axis.title.y=element_text(vjust=-0.25, size=40, family="Times New Roman"), 
        legend.text= element_text(size=23, family="Times New Roman"), 
        legend.title = element_text(size=23, family="Times New Roman"))

j22=ggplot(jap22, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=mut), size=7) + 
  ggtitle("s67" ) + 
  theme_bw()+ 
  xlab("Sapovirus genome position") + 
    ylab("iSNP frequency")+ scale_color_manual(values = c("syn" = "#00BFC4", "non_syn" = "red"))+
  #scale_color_manual(values = c("Synonymous" = "#00BFC4", "Non-synonymous" = "red"))+
  theme(legend.position = "none",
        panel.background = element_rect (fill = '#E7E7E7')) + 
  labs(x=NULL, y=NULL)+ expand_limits(x = c(1, 7300))+ 
  theme(plot.title=element_text(size=20, vjust=0.5, family="Times New Roman"),  
        axis.text.x = element_text(size=14, vjust=0.5, family="Times New Roman"), 
        axis.text.y = element_text(size=14, vjust=0.5, family="Times New Roman"), 
        axis.title.x=element_text(vjust=-0.25, size=40, family="Times New Roman"),
        axis.title.y=element_text(vjust=-0.25, size=40, family="Times New Roman"), 
        legend.text= element_text(size=23, family="Times New Roman"), 
        legend.title = element_text(size=23, family="Times New Roman"))

j23=ggplot(jap23, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=mut), size=7) + 
  ggtitle("s68" ) + 
  theme_bw()+
  scale_color_manual(values = c("syn" = "#00BFC4", "non_syn" = "red"))+
  #scale_color_manual(values = c("Synonymous" = "#00BFC4", "Non-synonymous" = "red"))+
  theme(legend.position = "none",panel.background = element_rect (fill = '#E7E7E7')) + 
  labs(x=NULL, y="iSNP frequency")+ expand_limits(x = c(1, 7300))+ 
  theme(plot.title=element_text(size=20, vjust=0.5, family="Times New Roman"),  
        axis.text.x = element_text(size=14, vjust=0.5, family="Times New Roman"), 
        axis.text.y = element_text(size=14, vjust=0.5, family="Times New Roman"), 
        axis.title.x=element_text(vjust=-0.25, size=40, family="Times New Roman"),
        axis.title.y=element_text(vjust=-0.25, size=40, family="Times New Roman"), 
        legend.text= element_text(size=23, family="Times New Roman"), 
        legend.title = element_text(size=23, family="Times New Roman"))

j24=ggplot(jap24, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=mut), size=7) + 
  ggtitle("s69" ) + 
  theme_bw()+
  xlab("Sapovirus genome position") + 
    ylab("iSNP frequency")+ scale_color_manual(values = c("syn" = "#00BFC4", "non_syn" = "red"))+
  #scale_color_manual(values = c("Synonymous" = "#00BFC4", "Non-synonymous" = "red"))+
  theme(legend.position = "none", panel.background = element_rect (fill = '#E7E7E7')) + 
  labs(x=NULL, y=NULL) + expand_limits(x = c(1, 7300))+ 
  theme(plot.title=element_text(size=20, vjust=0.5, family="Times New Roman"),  
        axis.text.x = element_text(size=14, vjust=0.5, family="Times New Roman"), 
        axis.text.y = element_text(size=14, vjust=0.5, family="Times New Roman"), 
        axis.title.x=element_text(vjust=-0.25, size=40, family="Times New Roman"),
        axis.title.y=element_text(vjust=-0.25, size=40, family="Times New Roman"), 
        legend.text= element_text(size=23, family="Times New Roman"), 
        legend.title = element_text(size=23, family="Times New Roman"))
 
 
j25= ggplot(jap25, aes(x=POS, y = ALT_FREQ))+ 
  geom_point(aes(color=mut), size=7) + ggtitle("s70" ) + 
  theme_bw()+ labs(x="Sapovirus genome position", y=NULL)+
  xlab("Sapovirus genome position") + 
  scale_color_manual(values = c("syn" = "#00BFC4", "non_syn" = "red"))+
  theme(legend.position = "right", 
        panel.background = element_rect (fill = '#E7E7E7'), 
        axis.title.x = element_text(size = 20)) + expand_limits(x = c(1, 7300)) + 
  theme(plot.title=element_text(size=20, vjust=0.5, family="Times New Roman"),  
        axis.text.x = element_text(size=14, vjust=0.5, family="Times New Roman"), 
        axis.text.y = element_text(size=14, vjust=0.5, family="Times New Roman"), 
        axis.title.x=element_text(vjust=-0.25, size=40, family="Times New Roman"),
        axis.title.y=element_text(vjust=-0.25, size=40, family="Times New Roman"), 
        legend.text= element_text(size=23, family="Times New Roman"), 
        legend.title = element_text(size=23, family="Times New Roman"))







# plot GI.1
library(patchwork)
GI.1 <- j4/j5/j6/j7/j8
GI.1 + 
  plot_annotation(title = 'Mutations compared with 2015 sapovirus GI.1 outbreak reference',
                  theme = theme(plot.title = element_text(size = 18)))& 
  theme(text = element_text("serif"))

ggsave(filename = "GI.1 outbreak SNPS remapped.jpg", 
       width = 25, height = 15, dpi = 350, units = "cm")


# plot GI.2 

GI.2 <- j21/j22/j23/j24/j25
GI.2 + 
  plot_annotation(title = 'Mutations compared with 2014 sapovirus GI.2 outbreak reference',
                  theme = theme(plot.title = element_text(size = 18)))& 
  theme(text = element_text("serif"))

ggsave(filename = "GI.2 outbreak SNPS remapped.jpg", 
       width = 25, height = 15, dpi = 350, units = "cm")

GI.1 | GI.2 + plot_layout(guides = "collect") 


ggsave(filename = "outbreak SNPS remapped.jpg", 
       width = 50, height = 25, dpi = 350, units = "cm")





data <- 
  readXL("/home/viro102/Documents/UGENE_Data/outbreak_interhost_SNPs.xlsx", 
         rownames=FALSE, header=TRUE, na="", sheet="Sheet1", stringsAsFactors=TRUE)
require("ggplot2")
.df <- data.frame(x = data$genotype, y = data$Ratio.NS.S)
.plot <- ggplot(data = .df, aes(x = factor(x), y = y)) + 
  stat_boxplot(geom = "errorbar", width = 0.5) + 
  geom_boxplot(outlier.colour = "transparent") + 
  geom_jitter(colour = "black", width = 0.1, height = 0) + 
  xlab("Genotype") + 
  ylab("ratio non-synonymoys/synonymous") + 
  theme_bw(base_size = 16, base_family = "sans")
print(.plot)
# rm(.df, .plot)to remove all plots and dataframes
