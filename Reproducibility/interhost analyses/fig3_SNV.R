setwd("C:/Users/Virology Tohoku/Google Drive/Manuscripts/My own/final edits/new outbreak data/tsv files")
## data cleaning 
data <- readr::read_tsv("1_japan_GI.1_final.tsv")
data <- data %>% filter(TOTAL_DP >= 400 & ALT_FREQ >= 0.03) %>% 
  select(-ALT_DP,-ALT_RV, -ALT_QUAL,-REF_DP,-REF_RV,-PASS,-REGION, -ALT, -REF, -REF_QUAL)%>% 
  mutate(mut = ifelse(REF_AA==ALT_AA, "syn", "non_syn"))
table(data$mut) # write the number in the file "outbreak_interhost_SNPs.xlsx"


# 1-group all outbreaks together -------------------------------------------------------

# to plot figure on interhost
setwd("C:/Users/Virology Tohoku/Google Drive/Manuscripts/My own/final edits/new outbreak data")


data <-  readxl::read_xlsx("outbreak_interhost_SNPs.xlsx")
library(ggpubr)
library(dplyr)
data$total.SNV = data$Synonymous + data$`non synonymous`
Jap.ratio.SNVs <- ggboxplot(data, x = "genotype", y = "Ratio NS/S",
                color = "genotype", palette = "lancet",
                xlab = "Genotype", ylab=("Ratio of non-synonymous/synonymous SNVs"),
                add = "jitter", bxp.errorbar = T, 
                ggtheme = theme_bw()) + 
  theme(text=element_text(size=18,  family="serif"), legend.position = "none") + 
  stat_compare_means(paired = F, label= "p.format" ,
                     label.x = 1.39 , show.legend=F,
                     label.y = 0.69) +
  #stat_compare_means( aes(label = ..p.signif..), label.x = 1.5 ,label.y = 0.25 ) + 
  font("xlab", size=18)+ font("ylab", size=18) + font("xy.text", size=18)+
  # Add a horizontal line segment
  geom_segment(aes(x = 1, y = 0.67, xend = 2, yend = 0.67))+
  # Add vertical line segment
  geom_segment(aes(x = 1, y = 0.67, xend = 1, yend = 0.65))+ 
  geom_segment(aes(x = 2, y = 0.67, xend = 2, yend = 0.65))



Jap.total.SNVs <- ggboxplot(data, x = "genotype", y = "total.SNV",
                      color = "genotype", palette = "lancet", 
                      xlab = "Genotype", ylab=("total SNV count"),
                      add = "jitter", bxp.errorbar = T, 
                      ggtheme = theme_bw()) + 
  theme(text=element_text(size=18,  family="serif"), legend.position = "none") + 
  stat_compare_means(paired = F, label= "p.format" ,
                     label.x = 1.39 , show.legend=F, label.y = 260)+
  #stat_compare_means( aes(label = ..p.signif..), label.x = 1.5 ,label.y = 0.25 ) + 
  font("xlab", size=18)+ font("ylab", size=18) + font("xy.text", size=18)+
  # Add a horizontal line segment
 geom_segment(aes(x = 1, y = 250, xend = 2, yend = 250))+
# Add vertical line segment
  geom_segment(aes(x = 1, y = 250, xend = 1, yend = 242))+ 
  geom_segment(aes(x = 2, y = 250, xend = 2, yend = 242))

## Fig 6

library(patchwork)

Jap.ratio.SNVs + theme(plot.margin = unit(c(0,60,0,0), "pt")) + 
  Jap.total.SNVs + plot_layout() + plot_annotation(tag_levels = "A")
ggsave(filename = "miyagi_outbreak_mutations1_byoutbreak.jpg", 
       width = 25, height = 15, dpi = 350, units = "cm")






# 2.Color by outbreak -------------------------------------------------------


Jap.ratio.SNVs <- ggboxplot(data, x = "genotype", y = "Ratio NS/S",
                            color = "Outbreak", palette = "lancet",
                            xlab = "Genotype", ylab=("Ratio of non-synonymous/synonymous SNVs"),
                            add = "jitter", bxp.errorbar = T, 
                            ggtheme = theme_bw()) + 
  theme(text=element_text(size=18,  family="serif"), legend.position = "right") + 
  stat_compare_means(paired = F, label= "p.format" ,
                     label.x = 1.39 , show.legend=F,
                     label.y = 0.69) +
  #stat_compare_means( aes(label = ..p.signif..), label.x = 1.5 ,label.y = 0.25 ) + 
  font("xlab", size=18)+ font("ylab", size=18) + font("xy.text", size=18)+
  # Add a horizontal line segment
  geom_segment(aes(x = 1, y = 0.67, xend = 2, yend = 0.67))+
  # Add vertical line segment
  geom_segment(aes(x = 1, y = 0.67, xend = 1, yend = 0.65))+ 
  geom_segment(aes(x = 2, y = 0.67, xend = 2, yend = 0.65))


Jap.total.SNVs <- ggboxplot(data, x = "genotype", y = "total.SNV",
                            color = "Outbreak", palette = "lancet", 
                            xlab = "Genotype", ylab=("total SNV count"),
                            add = "jitter", bxp.errorbar = T, 
                            ggtheme = theme_bw()) + 
  theme(text=element_text(size=18,  family="serif"), legend.position = "none") + 
  stat_compare_means(paired = F, label= "p.format" ,
                     label.x = 1.39 , show.legend=F, label.y = 260)+
  #stat_compare_means( aes(label = ..p.signif..), label.x = 1.5 ,label.y = 0.25 ) + 
  font("xlab", size=18)+ font("ylab", size=18) + font("xy.text", size=18)+
  # Add a horizontal line segment
  geom_segment(aes(x = 1, y = 250, xend = 2, yend = 250))+
  # Add vertical line segment
  geom_segment(aes(x = 1, y = 250, xend = 1, yend = 242))+ 
  geom_segment(aes(x = 2, y = 250, xend = 2, yend = 242))

## Fig 6

library(patchwork)

Jap.total.SNVs + theme(plot.margin = unit(c(0,60,0,0), "pt")) + 
  Jap.ratio.SNVs + plot_layout() + plot_annotation(tag_levels = "A")
ggsave(filename = "miyagi_outbreak_mutations1_byoutbreak.jpg", 
       width = 25, height = 15, dpi = 350, units = "cm")


# 2.jitter points by outbreak -------------------------------------------------------
# use this for paper

Jap.ratio.SNVs <- ggstripchart(data, x = "genotype", y = "Ratio NS/S",size = 4,
                            shape = "Outbreak", color = "Outbreak", palette = "lancet",
                            xlab = "Genotype", ylab=("Ratio of non-synonymous/synonymous SNV"),
                            add = "boxplot", add.params = list(color = "black"), 
                            ggtheme = theme_bw()) + 
  theme(text=element_text(size=18,  family="serif"), legend.position = "right") + 
  stat_compare_means(paired = F, label= "p.format" ,
                     label.x = 1.39 , show.legend=F,
                     label.y = 0.69) +
  #stat_compare_means( aes(label = ..p.signif..), label.x = 1.5 ,label.y = 0.25 ) + 
  font("xlab", size=18)+ font("ylab", size=18) + font("xy.text", size=18)+
  # Add a horizontal line segment
  geom_segment(aes(x = 1, y = 0.67, xend = 2, yend = 0.67))+
  # Add vertical line segment
  geom_segment(aes(x = 1, y = 0.67, xend = 1, yend = 0.65))+ 
  geom_segment(aes(x = 2, y = 0.67, xend = 2, yend = 0.65))
 
Jap.total.SNVs <- ggstripchart(data, "genotype", "total.SNV", shape = "Outbreak",
             color = "Outbreak", palette = "lancet", 
             xlab = "Genotype", ylab=("total SNV count"),
             add = "boxplot", bxp.errorbar = T, size = 4,
             add.params = list(color = "black"), 
  ggtheme = theme_bw()) + 
 theme(text= element_text(size=18,  family="serif"), 
       legend.position = "right")+ 
       #legend.background = element_rect(fill=NULL, colour = "black", size=0.5, linetype="solid")) + 
  stat_compare_means(paired = F, label= "p.format" ,
                     label.x = 1.39 , show.legend=F, label.y = 260)+
  #stat_compare_means( aes(label = ..p.signif..), label.x = 1.5 ,label.y = 0.25 ) + 
  font("xlab", size=18)+ font("ylab", size=18) + font("xy.text", size=18)+
  # Add a horizontal line segment
  geom_segment(aes(x = 1, y = 250, xend = 2, yend = 250))+
  # Add vertical line segment
  geom_segment(aes(x = 1, y = 250, xend = 1, yend = 242))+ 
  geom_segment(aes(x = 2, y = 250, xend = 2, yend = 242))

Jap.total.SNVs + theme(plot.margin = unit(c(0,30,0,0), "pt")) + 
  Jap.ratio.SNVs + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")
ggsave(filename = "miyagi_outbreak_strip_OK.jpg", 
       width = 25, height = 15, dpi = 350, units = "cm")









# 3 - Fig3-remapped-good data for paper ---------------------------------------


  #jitter points by outbreak 
# use this for paper


# to plot figure on interhost excuding sample 
setwd("C:/Users/Virology Tohoku/Google Drive/Manuscripts/My own/final edits/new outbreak data")


data <-  readxl::read_xlsx("outbreak_interhost_SNPs.xlsx", sheet = 2)
library(ggpubr)
library(dplyr)
data$total.SNV = data$Synonymous + data$`non synonymous`

Jap.ratio.SNVs <- ggstripchart(data, x = "genotype", y = "Ratio NS/S",size = 4,
                               shape = "Outbreak", color = "Outbreak", palette = "lancet",
                               xlab = "Genotype", ylab=("Ratio of non-synonymous/synonymous SNV"),
                               add = "boxplot", add.params = list(color = "black"), 
                               ggtheme = theme_bw()) + 
  theme(text=element_text(size=18,  family="serif"), legend.position = "right") + 
  stat_compare_means(paired = F, label= "p.format" ,
                     label.x = 1.39 , show.legend=F,
                     label.y = 0.69) +
  #stat_compare_means( aes(label = ..p.signif..), label.x = 1.5 ,label.y = 0.25 ) + 
  font("xlab", size=18)+ font("ylab", size=18) + font("xy.text", size=18)+
  # Add a horizontal line segment
  geom_segment(aes(x = 1, y = 0.67, xend = 2, yend = 0.67))+
  # Add vertical line segment
  geom_segment(aes(x = 1, y = 0.67, xend = 1, yend = 0.65))+ 
  geom_segment(aes(x = 2, y = 0.67, xend = 2, yend = 0.65))

Jap.total.SNVs <- ggstripchart(data, "genotype", "total.SNV", shape = "Outbreak",
                               color = "Outbreak", palette = "lancet", 
                               xlab = "Genotype", ylab=("total SNV count"),
                               add = "boxplot", bxp.errorbar = T, size = 4,
                               add.params = list(color = "black"), 
                               ggtheme = theme_bw()) + 
  theme(text= element_text(size=18,  family="serif"), 
        legend.position = "right")+ 
  #legend.background = element_rect(fill=NULL, colour = "black", size=0.5, linetype="solid")) + 
  stat_compare_means(paired = F, label= "p.format" ,
                     label.x = 1.39 , show.legend=F, label.y = 260)+
  #stat_compare_means( aes(label = ..p.signif..), label.x = 1.5 ,label.y = 0.25 ) + 
  font("xlab", size=18)+ font("ylab", size=18) + font("xy.text", size=18)+
  # Add a horizontal line segment
  geom_segment(aes(x = 1, y = 250, xend = 2, yend = 250))+
  # Add vertical line segment
  geom_segment(aes(x = 1, y = 250, xend = 1, yend = 242))+ 
  geom_segment(aes(x = 2, y = 250, xend = 2, yend = 242))

library(patchwork)
Jap.total.SNVs + theme(plot.margin = unit(c(0,30,0,0), "pt")) + 
  Jap.ratio.SNVs + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")


# LEt us exclude low breath samples


Jap.ratio.SNVs <- ggstripchart(data[data$breadth_400... > 50, ], x = "genotype", y = "Ratio.NS.S" ,size = 4,
                               shape = "Outbreak", color = "Outbreak", palette = "lancet",
                               xlab = "Genotype", ylab=("Ratio of non-synonymous/synonymous SNV"),
                               add = "boxplot", add.params = list(color = "black"), 
                               ggtheme = theme_bw()) + 
  theme(text=element_text(size=18,  family="serif"), legend.position = "right") + 
  stat_compare_means(paired = F, label= "p.format" ,
                     label.x = 1.39 , show.legend=F,
                     label.y = 0.69) +
  #stat_compare_means( aes(label = ..p.signif..), label.x = 1.5 ,label.y = 0.25 ) + 
  font("xlab", size=18)+ font("ylab", size=18) + font("xy.text", size=18)+
  # Add a horizontal line segment
  geom_segment(aes(x = 1, y = 0.67, xend = 2, yend = 0.67))+
  # Add vertical line segment
  geom_segment(aes(x = 1, y = 0.67, xend = 1, yend = 0.65))+ 
  geom_segment(aes(x = 2, y = 0.67, xend = 2, yend = 0.65))

Jap.total.SNVs <- ggstripchart(data[data$breadth_400... > 50, ], "genotype", "total.SNV", shape = "Outbreak",
                               color = "Outbreak", palette = "lancet", 
                               xlab = "Genotype", ylab=("total SNV count"),
                               add = "boxplot", bxp.errorbar = T, size = 4,
                               add.params = list(color = "black"), 
                               ggtheme = theme_bw()) + 
  theme(text= element_text(size=18,  family="serif"), 
        legend.position = "right")+ 
  #legend.background = element_rect(fill=NULL, colour = "black", size=0.5, linetype="solid")) + 
  stat_compare_means(paired = F, label= "p.format" ,
                     label.x = 1.39 , show.legend=F, label.y = 260)+
  #stat_compare_means( aes(label = ..p.signif..), label.x = 1.5 ,label.y = 0.25 ) + 
  font("xlab", size=18)+ font("ylab", size=18) + font("xy.text", size=18)+
  # Add a horizontal line segment
  geom_segment(aes(x = 1, y = 250, xend = 2, yend = 250))+
  # Add vertical line segment
  geom_segment(aes(x = 1, y = 250, xend = 1, yend = 242))+ 
  geom_segment(aes(x = 2, y = 250, xend = 2, yend = 242))

library(patchwork)
Jap.total.SNVs + theme(plot.margin = unit(c(0,30,0,0), "pt")) + 
  Jap.ratio.SNVs + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")


ggsave(filename = "miyagi_outbreak_secondary_excluding low breath.jpg", 
       width = 25, height = 15, dpi = 350, units = "cm")

