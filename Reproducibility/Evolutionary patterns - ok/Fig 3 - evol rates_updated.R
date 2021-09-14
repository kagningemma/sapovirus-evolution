

# plotting MCMC estimates - substitution rates ----------------------------


# MCMC plots --------------------------------------------------------------
library(ggpubr)
library(tidyverse)


# polymerase

g1.1polylogdata <- read.table("C:/Users/Virology Tohoku/Google Drive/MyRdirectory/SAPO/Sapovirus evolution over time/my data+tohma_data/polymerase/v.1.8.2/GI.1/GI.1_n29_polym_v1.8.log.txt", header = TRUE)# no row labels required!
g1.2polylogdata <- read.table("C:/Users/Virology Tohoku/Google Drive/MyRdirectory/SAPO/Sapovirus evolution over time/my data+tohma_data/polymerase/v.1.8.2/GI.2/GI.2_n13_polym_trim3150_5170.log.txt", header = TRUE)# no row labels required!

#g1.1caplogdata <- read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/mcmc files/capsid/sav.GI.1RLC2.log", header = TRUE)# no row labels required!
#g1.2capylogdata <- read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/mcmc files/capsid/sav.GI.2_RCE.log", header = TRUE)# no row labels required!

g1.1caplogdata <- read.table("C:/Users/Virology Tohoku/Google Drive/MyRdirectory/SAPO/Sapovirus evolution over time/my data+tohma_data/VP1/BEAST/GI.1 NT/GI.1_VP1_NT_inframe_n60.log.txt", header = TRUE)# no row labels required!
g1.2capylogdata <- read.table("C:/Users/Virology Tohoku/Google Drive/MyRdirectory/SAPO/Sapovirus evolution over time/my data+tohma_data/VP1/BEAST/GI.2 NT/GI.2_VP1_NT_inframe_n56.log.txt", header = TRUE)# no row labels required!




str(g1.2polylogdata) # clockRate
str(g1.1caplogdata) # meanClockRate   / clockRate
str(g1.2capylogdata)# ucedMean
# arranging the dataset for plotting
 
g1.1ratecap <- as.data.frame(g1.1caplogdata$clock.rate)
g1.2ratecap <- as.data.frame(g1.2capylogdata$clock.rate)

library(gdata)
cap.rate <- cbindX(g1.1ratecap, g1.2ratecap)
str(cap.rate)
head(cap.rate)
#cap.rate <- cap.rate %>% dplyr::mutate( 
# GI.1= g1.1ratecap $ meanClockRate, 
#GI.2 = g1.2ratecap $ ucedMean)
colnames(cap.rate)<- c("GI.1", "GI.2")
cap.rate.r <- reshape2::melt(cap.rate, value.name="clock.rate")
cap.rate.r <- tibble::as_tibble(cap.rate.r %>% filter(clock.rate < 1)) # exclude the first value
names(cap.rate.r)[c(1)] <- c("genotypes")
head(cap.rate.r)
str(cap.rate.r)
# width: change box plots width
library(ggpubr)
ggboxplot(cap.rate.r, x = "genotypes", y = "clock.rate", width = 0.8)
library(scales)

blue.bold.italic.16.text <- element_text(color = "black", size = 15)

## axis.text.x for x axis only  | Without outliers

cap <- ggboxplot(cap.rate.r, x = "genotypes", y = "clock.rate", 
                 width = 0.8, color = "genotypes", palette = "lancet", outlier.shape = NA,
                 xlab = "Genotype",  ylab= "substitutions/site/year", bxp.errorbar = T)+ 
  theme_bw() + font("xlab", size=15)+ font("ylab", size=15) + font("xy.text", size=13) + 
  labs(title = "VP1") + theme_bw(base_size = 18, base_family = "serif")+  
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.title.y = element_blank()) + 
  theme(axis.text.x = blue.bold.italic.16.text, axis.text.y = blue.bold.italic.16.text ) + 
  geom_segment(aes(x = 2, y = 0.002, xend = 2, yend = 0.0021)) + 
  geom_segment(aes(x = 1, y = 0.0021, xend = 1, yend = 0.0021)) + 
  geom_segment(aes(x = 1, y = 0.0021, xend = 2, yend = 0.0021)) + 
  annotate("text", x = 1.5, y = 0.00207, parse = TRUE, size = 4,
           label = "'p-value' == 2 %*% 10^{-16}")
  
# without annotations 
cap <- ggboxplot(cap.rate.r, x = "genotypes", y = "clock.rate", 
                 width = 0.8, color = "genotypes", palette = "lancet", outlier.shape = NA,
                 xlab = "Genotype",  ylab= "substitutions/site/year", bxp.errorbar = T)+ 
  theme_bw() + font("xlab", size=15)+ font("ylab", size=15) + font("xy.text", size=13) + 
  labs(title = "VP1") + theme_bw(base_size = 18, base_family = "serif")+  
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.title.y = element_blank()) + 
  theme(axis.text.x = blue.bold.italic.16.text, axis.text.y = blue.bold.italic.16.text ) 

# See ?plotmath for many examples of mathematical expressions, and ?demo(plotmath) for graphical examples of mathematical expressions.

# to determine the pvalue : +   stat_compare_means( aes(label = ..p.format..), label.x = 1 ,label.y = 0.001 ) 

cap


# cap <- cap + stat_compare_means(aes(label = ..p.signif..), label.x = 1.5 ,label.y = 0.0037) 



# arranging poly dataseet
library(dplyr)
str(g1.1polylogdata)
# arranging the dataset for plotting
g1.1ratepoly <- g1.1polylogdata %>% select(clock.rate)
g1.2ratepoly <- g1.2polylogdata%>% select(clock.rate)
poly.rate <- data.frame(g1.1ratepoly, g1.2ratepoly)
colnames(poly.rate)<- c("GI.1", "GI.2")
poly.rate.r <- reshape2::melt(poly.rate, value.name="clock.rate")
poly.rate.r <- tibble::as_tibble(poly.rate.r %>% filter(clock.rate < 1))
names(poly.rate.r)[c(1)] <- c("genotypes")


str(poly.rate.r)
# width: change box plots width
ggboxplot(poly.rate.r, x = "genotypes", y = "clock.rate", width = 0.8)

str(poly.rate.r)



# pol<- pol + stat_compare_means(aes(label = ..p.signif.., size=100), label.x = 1.5 ,label.y = 0.0039) 



pol <- ggboxplot(poly.rate.r, x = "genotypes", y = "clock.rate", 
          width = 0.8, color = "genotypes", palette = "lancet", outlier.shape = NA,
          xlab = "Genotype",  ylab= "Substitutions/site/year", bxp.errorbar = T)+ 
  theme_bw() + font("xlab", size=15)+ font("ylab", size=15) + font("xy.text", size=13) + 
  labs(title = "Polymerase") + theme_bw(base_size = 18, base_family = "serif")+  
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = blue.bold.italic.16.text, axis.text.y = blue.bold.italic.16.text ) + 
  geom_segment(aes(x = 2, y = 0.0053, xend = 2, yend = 0.0055)) + 
  geom_segment(aes(x = 1, y = 0.0053, xend = 1, yend = 0.0055)) + 
  geom_segment(aes(x = 1, y = 0.0055, xend = 2, yend = 0.0055)) + 
  annotate("text", x = 1.5, y = 0.0056, parse = TRUE,  size = 4, label = "'p-value' == 2 %*% 10^{-16}")


# pol + stat_compare_means(aes(label = ..p.format..), label.x = 1.5 ,label.y = 0.0037) 
pol



# without annotations

pol <- ggboxplot(poly.rate.r, x = "genotypes", y = "clock.rate", 
                 width = 0.8, color = "genotypes", palette = "lancet", outlier.shape = NA,
                 xlab = "Genotype",  ylab= "Substitutions/site/year", bxp.errorbar = T)+ 
  theme_bw() + font("xlab", size=15)+ font("ylab", size=15) + font("xy.text", size=13) + 
  labs(title = "Polymerase") + theme_bw(base_size = 18, base_family = "serif")+  
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = blue.bold.italic.16.text, axis.text.y = blue.bold.italic.16.text )

library(patchwork)

pol|cap

ggsave(filename = "sapo.cap.pol.evol.rates.tiff", 
       width = 22, height = 15, dpi = 300, units = "cm")

