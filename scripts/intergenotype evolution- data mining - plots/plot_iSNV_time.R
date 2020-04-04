# Analysis of summary distribution in peru data

library(Rcmdr)
data <- 
  readXL("C:/Users/Virology Tohoku/Dropbox/phil_data/PRIMAL_sample_conc_ms_saka_kte.xlsx",
         rownames=FALSE, header=TRUE, na="", sheet="suppl-peru", 
         stringsAsFactors=TRUE)
summary(data)


require("ggplot2")
.df <- data.frame(x = data$Day.difference, y = data$raw.number.of.SNPs, z = data$genotype)
.plot <- ggplot(data = .df, aes(x = x, y = y, colour = z, shape = z)) + 
  geom_point(size=4) + 
  stat_smooth(aes(fill = z), method = "lm") + 
  xlab("Day difference") + 
  ylab("cumulative iSNVs") + 
  labs(colour = "genotype", shape = "genotype", fill = "genotype") + 
  theme_bw(base_size = 14, base_family = "sans") + 
  theme(legend.position = "right")
print(.plot)



.df <- data.frame(x = data$Day.difference, y = 
                    data$raw.number.of.SNPs, z = data$genotype)
.plot <- ggplot(data = .df, aes(x = x, y = y, colour = z, shape = z)) + geom_point() + stat_smooth(aes(fill = z), method = "lm", se = FALSE) + 
  xlab("day difference") + ylab("cumulative iSNVs") + labs(colour = "genotype", shape = "genotype", fill = "genotype")
+ theme(legend.position = "right")+theme_bw(base_size = 14, base_family = "sans")

print(.plot)
# ggplot_build(.plot)

setwd("C:/Users/Virology Tohoku/Dropbox/manuscript-plos-pathogens/supplements")
ggsave(filename = "iSNVs with time222.tiff", 
       width = 22, height = 13, dpi = 300, units = "cm")



