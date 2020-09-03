

# plotting MCMC estimates - substitution rates ----------------------------


# MCMC plots --------------------------------------------------------------
library(ggpubr)
library(dplyr)
# polymerase
g1.1polylogdata <- read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/mcmc files/poly/GI.1_polym_strict_n23.log", header = TRUE)# no row labels required!
g1.2polylogdata <- read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/mcmc files/poly/GI.2_polym_strict_aa_n17.log", header = TRUE)# no row labels required!
# capsid
g1.1caplogdata <- read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/mcmc files/capsid/sav.GI.1RLC2.log", header = TRUE)# no row labels required!
g1.2capylogdata <- read.table("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/mcmc files/capsid/sav.GI.2_RCE.log", header = TRUE)# no row labels required!

str(g1.2polylogdata)
str(g1.1caplogdata) # meanClockRate
str(g1.2capylogdata)# ucedMean
# arranging the dataset for plotting

g1.1ratecap <- as.data.frame(g1.1caplogdata$meanClockRate)
g1.2ratecap <- as.data.frame(g1.2capylogdata$ucedMean)

library(gdata)
cap.rate <- cbindX(g1.1ratecap, g1.2ratecap)
str(cap.rate)
head(cap.rate)
#cap.rate <- cap.rate %>% dplyr::mutate( 
# GI.1= g1.1ratecap $ meanClockRate, 
#GI.2 = g1.2ratecap $ ucedMean)
colnames(cap.rate)<- c("GI.1", "GI.2")
cap.rate.r <- reshape2::melt(cap.rate, value.name="clock.rate")
cap.rate.r <- tibble::as_tibble(cap.rate.r %>% filter(clock.rate < 0.0035))
names(cap.rate.r)[c(1)] <- c("genotypes")
head(cap.rate.r)
str(cap.rate.r)
# width: change box plots width
ggboxplot(cap.rate.r, x = "genotypes", y = "clock.rate", width = 0.8)
library(scales)

cap <- ggboxplot(cap.rate.r, x = "genotypes", y = "clock.rate", 
                 width = 0.8, color = "genotypes", palette = "lancet", 
                 xlab = "Genotype",  
                 ylab= "substitutions/site/year", bxp.errorbar = T)+ 
  theme_bw() + font("xlab", size=15)+ font("ylab", size=15) + font("xy.text", size=13)

# cap <- cap + stat_compare_means(aes(label = ..p.signif..), label.x = 1.5 ,label.y = 0.0037) 



# arranging poly dataseet
library(dplyr)
head(g1.1polylogdata)
# arranging the dataset for plotting
g1.1ratepoly <- g1.1polylogdata %>% select(clockRate)
g1.2ratepoly <- g1.2polylogdata%>% select(clockRate)
poly.rate <- data.frame(g1.1ratepoly, g1.2ratepoly)
colnames(poly.rate)<- c("GI.1", "GI.2")
poly.rate.r <- reshape2::melt(poly.rate, value.name="clock.rate")
poly.rate.r <- tibble::as_tibble(poly.rate.r %>% filter(clock.rate < 0.0038))
names(poly.rate.r)[c(1)] <- c("genotypes")


str(poly.rate.r)
# width: change box plots width
ggboxplot(poly.rate.r, x = "genotypes", y = "clock.rate", width = 0.8)
class(df)
str(poly.rate.r)


pol <- ggboxplot(poly.rate.r, x = "genotypes", y = "clock.rate", 
                 width = 0.8, color = "genotypes", palette = "lancet", 
                 xlab = "Genotype",  
                 ylab= "substitutions/site/year", bxp.errorbar = T)+ 
  theme_bw()+ font("xlab", size=15)+ font("ylab", size=15) + font("xy.text", size=13)


pol<- pol + stat_compare_means(aes(label = ..p.signif.., size=100), label.x = 1.5 ,label.y = 0.0039) 



library(patchwork)

pol|cap

ggsave(filename = "sapo.cap.pol.evol.rates.tiff", 
       width = 22, height = 15, dpi = 300, units = "cm")

