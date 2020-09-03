


library(dplyr)
library(readxl)

# Figure 1 --------------------------
## Plots of diversification over time

# for Sapovirus GI ---------------------------------------------------
setwd("~/MyRdirectory/SAPO/Sapovirus evolution over time")
data = read_xlsx("SaVaminoacid_diffbigdata_4plotting.xlsx", sheet=1, col_names = T)

data <- data %>% dplyr::rename(specie.1="Species 1", specie.2= "Species 2")

str_extract(data$specie.1,"[G][I-IV|1-5][.| ][0-9 | I-IV]") # to exttact g1, g2 , g4 g5
table(str_extract(data$specie.1,"[G][I-IV|1-5][.| ][0-9 | I-IV]"))
# how to use stringi https://r4ds.had.co.nz/strings.html
str_extract(data$specie.1[1],"[G][I|1][.| ][0-9 | I-IV]") # to select G1 only
table(str_extract(data$specie.1, "[G][I|1][.| ][0-9 | I-IV]"))

data$gen_spe_1 = str_extract(data$specie.1, "[G][I|1][.| ][0-9 | I-IV]")
data$gen_spe_2 = str_extract(data$specie.2, "[G][I|1][.| ][0-9 | I-IV]")# to extract the genotypes

# table(data$specie.1) %>% table(data$specie.2)

table(data$gen_spe_2)
table(data$gen_spe_2)

# fix the levels 

data$gen_spe_1[data$gen_spe_1=="GI 2"] <- "GI.2"
data$gen_spe_1[data$gen_spe_1=="GI.I"] <- "GI.1"
data$gen_spe_1[data$gen_spe_1=="G1.1"] <- "GI.1"
table(data$gen_spe_1)

data$gen_spe_2[data$gen_spe_2=="GI 2"] <- "GI.2"
data$gen_spe_2[data$gen_spe_2=="GI.I"] <- "GI.1"
data$gen_spe_2[data$gen_spe_2=="G1.1"] <- "GI.1"
table(data$gen_spe_2)

table(is.na(data$gen_spe_2))

#FALSE  TRUE 
#5405 22561 
table(is.na(data$gen_spe_1))
#FALSE  TRUE 
#18667  9299 

data<- data[complete.cases(data), ] # to keep only lines without NAs

data$label<- paste(data$gen_spe_1,data$gen_spe_2, sep = "Vs")
ggplot(data, aes(aa_difference))+ 
  geom_density()+ggtitle(label = "distribution of pairwise distances for Sav G1-aa") 
# some samples have pdist >0.5 (amino acid distance >200), it is abnormal, let us identify those sequences and filter them out 

data%>% dplyr::filter(aa_difference>200) # to see all weird pairwie comparison with
# sapovirus GQ261222 in bangladesh is a recombinant sequence, thus i will remove it.
data <- data[!grepl("GQ261222.1 Sapovirus BD/697 BGD GI 2005", data$specie.1),] 
data <- data[!grepl("GQ261222.1 Sapovirus BD/697 BGD GI 2005", data$specie.2),]


# now i need to compute one column which represent the difference in isolation year 

# to extract the last 4 characters
data$year1 = as.numeric(str_sub(data$specie.1,-4,-1))
data$year2= as.numeric(str_sub(data$specie.2,-4,-1))

data$year_diff<- abs(data$year1-data$year2) # compute the difference
#ggplot(data, aes(year1))+geom_histogram()
#ggplot(data, aes(year2))+geom_histogram()
ggplot(data, aes(year_diff)) + geom_histogram() + ggtitle( label="histogram of isolation years difference of SaV GI")
hist(data$year_diff) # it is skewed to the left. 

# GI sequences were isolated over a period of 34 years. 

# # Gomogenotype analyses
# let us focus on same genotype compasisons
data$label
data2<- data %>% 
  dplyr::filter(label %in% c("GI.1VsGI.1", "GI.2VsGI.2", "GI.3VsGI.3", "GI.4VsGI.4", 
                             "GI.5VsGI.5", "GI.6VsGI.6", "GI.7VsGI.7", "GI.8VsGI.8"))

ggplot(data2, aes(aa_difference, year_diff)) +
  geom_point(aes(colour = label)) + 
  labs(title = "Labelled scatterplot of amino acid differences for genogroup 1")+
  xlab("amino acid difference")+
  ylab("detection year difference")+ 
  theme_bw(base_size = 15)

data2 %>% filter(aa_difference >90) # Chanthanburi GI.2 strain is very different from other GI.2 strains

ggpubr::ggscatter(data2, x="aa_difference",  y= "year_diff", 
                  color="label", palette = "lancet")




# for Sapovirus GI.1 and GI.2 heatmap ---------------------------------------------------

#  make a subset of the dataset containing only GI.1 and GI.2 pairwise distances

data2<- data %>% filter(label %in% c("GI.1VsGI.1", "GI.2VsGI.2" ))

# Heatmap for GI.1 only 
data1.1<- data2 %>% filter(label %in% c("GI.1VsGI.1"))# to make sure the labels are correct
table(data1.1$gen_spe_1,data1.2$gen_spe_2)
# Heatmap for GI.2 only 
data1.2<- data2 %>% filter(label %in% c("GI.2VsGI.2"))# to make sure the labels are correct
table(data1.2$gen_spe_1,data1.2$gen_spe_2)
require(gridExtra)

g1 <- ggplot(data1.1, aes(x=aa_difference, y=year_diff))+ 
  geom_bin2d() + 
  theme_bw(base_size = 15) + 
  xlab("Amino acid differences") + ylab("Detection year difference") + 
  scale_fill_gradient('counts', low = "blue", high = "red")+
  ggtitle("Sapovirus GI.1")+
  theme(plot.title=element_text(size=15, vjust=2, family="Times New Roman"),  
        axis.title.x=element_text(vjust=-1, size=15, family="Times New Roman"),
        axis.title.y=element_text(vjust=-0.25, size=15, family="Times New Roman"), 
        legend.text=element_text(size=15, family="Times New Roman"))

g2 <- ggplot(data1.2, aes(x=aa_difference, y=year_diff))+ 
  geom_bin2d() + 
  theme_bw(base_size = 15) + 
  xlab("Amino acid differences") + ylab("Detection year difference") + 
  scale_fill_gradient('counts', low = "blue", high = "red")+
  ggtitle("Sapovirus GI.2")+
  theme(plot.title=element_text(size=15, vjust=2, family="Times New Roman"),  
        axis.title.x=element_text(vjust=-1, size=15, family="Times New Roman"),
        axis.title.y=element_text(vjust=-0.25, size=15, family="Times New Roman"), 
        legend.text=element_text(size=15, family="Times New Roman"))

library(patchwork)

g1/g2 + plot_layout()+
  plot_annotation(tag_levels = 'A')



ggsave(filename = "sapovirus_gb_genetic_divers.tiff", 
       width = 10, height = 13, dpi = 350, units = "in")















# Computing the correlation coeficient for GI.1 and GI.2
# Correlation coeficient between aa_differnce and isolation year_diff of GI.1 and GI.2

data3<- data2 %>% dplyr::filter(label=="GI.1VsGI.1")
data4<- data2 %>% dplyr::filter(label=="GI.2VsGI.2")

hist(data3$aa_difference) # the distribution of this variable is skewed to the left
hist(data3$year_diff) # the distribution of this variable skewed to the left
qqnorm(data3$aa_difference, main = "Normal Q-Q Plot",
       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
       plot.it = TRUE, datax = FALSE)
qqline(data3$aa_difference, col = 2)

# skewed distribution to the left

cor.test(data3$aa_difference, data3$year_diff, method = "spearman", conf.level = 0.95) 
cor.test(data4$aa_difference, data4$year_diff, method = "spearman", conf.level = 0.95)

# 1) for SaV GI.1, r_SaV GI.2 = 0.36 (p-value < 2.2e-16)
# 2) for Sav GI.2 ,  # r_SaV GI.2 = 0.86 (p-value < 2.2e-16)












# for Sapovirus GII ---------------------------------------------------
setwd("~/MyRdirectory/SAPO/Sapovirus evolution over time")
library(readxl)
data = read_xlsx("SaVaminoacid_diffbigdata_4plotting.xlsx", sheet=1, col_names = T)
data <- data %>% dplyr::rename(specie.1="Species 1", specie.2= "Species 2")
table(data$specie.1) # GII X, GII.3, GII.2, GII.1, GII.6, GII.8, GII.3, GII.5, GII.4, 
table(str_extract(data$specie.1,"GII.."))# match GII and any2 following character
table(str_extract(data$specie.1,"GII.*"))
data2<- data %>% filter(grepl("GII.*", specie.1))  %>%  filter(grepl("GII.*", specie.2)) # to only retain row containing GII.. pattern
table(data2$specie.1) # i will need this if i need to focus on some specific groups
# str_extract(data$specie.1[1],"[G][II|2][.| ][0-9 | I-IV]") # to select G2 only
# table(str_extract(data$specie.1, "[G][I|1][.| ][0-9 | I-IV]"))

data$gen_spe_1 = str_extract(data$specie.1,"GII..")
data$gen_spe_2 = str_extract(data$specie.2,"GII..")# to extract the genotypes

table(data$gen_spe_1)
table(data$gen_spe_2)
which(data$gen_spe_1=="GII X")
which(data$gen_spe_2=="GII X")-> x
X<- data[x, ]
# data <- data %>% filter(data!=12)
table(is.na(data$gen_spe_2))
table(is.na(data$gen_spe_1))

data = data[!is.na(data$gen_spe_1), ] # to keep all rows of the 2nd col with no missing value
data = data[!is.na(data$gen_spe_2), ] # to keep all rows of the 2nd col with no missing value

data<- data[complete.cases(data), ] # to keep only lines without NAs

data <- data[!grepl("GQ261222.1 Sapovirus BD/697 BGD GI 2005", data$specie.1),] 
data <- data[!grepl("GQ261222.1 Sapovirus BD/697 BGD GI 2005", data$specie.2),]


data$gen_spe_1<- ifelse(data$gen_spe_1=="GII ." , "GII.8", data$gen_spe_1) # changing GII . to GII.8 
data$gen_spe_2<- ifelse(data$gen_spe_2=="GII ." , "GII.8", data$gen_spe_2)
glimpse(data)

data$label <- paste(data$gen_spe_1,data$gen_spe_2, sep = "vs")
class(data)

f<- ggplot(data, aes(aa_difference))
f+geom_density()+ggtitle(label = "Distribution of pairwise distances for Sav GII-aa") # bimodal distribution, showing variants, genotypes and genogroup distribution

# Here, everythings seems ok. OKa et al, 2012 reported that mean +- 3sd of each pairwise distance value for strain is 0-0.127, and for genotype is 0.2-0.37. 
# therefore, this distribution is in accordance with oka's observations.


# Let us compute one column which represent the difference in isolation year 

# to extract last 4 characters
data$year1 = as.numeric(str_sub(data$specie.1,-4,-1))
data$year2= as.numeric(str_sub(data$specie.2,-4,-1))

summary(data)
describe(data)
table(data$label)

data$year_diff<- abs(data$year1-data$year2)# compute the difference

ggplot(data, aes(year_diff)) + geom_histogram() + ggtitle(label="histogram of isolation years difference of SaV GII")
hist(data$year_diff)
hist(data$aa_difference)# it is skewed to the left. most samples were isolated during the same period
ggplot(data, aes(aa_difference))+ 
  geom_density()+ggtitle(label = "distribution of # AA differences for Sav GII")


#  scatter plot and label by points

ggplot(data, aes(aa_difference, year_diff)) +
  geom_point(aes(colour = label)) + 
  labs(title = "Evolutionary pattern of sapovirus Genogroup II")




## Genotype-wise analyses

table(data$specie.1)
table(data$gen_spe_1)
data <- data %>% mutate(homo_hetero  =  ifelse(data$gen_spe_1==data$gen_spe_2, "homogenotype", "heterogenotype")) # creating new variable
data %>% select(homo_hetero) %>% table() 

class(data)
ggplot(data, aes(aa_difference, year_diff))+ geom_point(aes(colour=homo_hetero, shape=homo_hetero)) + labs(title = "scatterplot Sav GII-aa diff-Homogenotypes-heterogenotypes")

#ggplotly(p)

# most homogenotypes have a aa diff <50 (99% homology) and there is a little trend if i only focus on variants. 
# below, I will only focus on those with pdist<0.1 (that variants resulting from pairwise comparison)
table(data$label)
library(dplyr)
data2<- data %>% 
  filter(label %in% c("GII.1vsGII.1", "GII.2vsGII.2", "GII.3vsGII.3", "GII.4vsGII.4", "GII.5vsGII.5", 
                      "GII.6vsGII.6", "GII.7vsGII.7", "GII.8vsGII.8", "GII XvsGII X"))

table(data2$label)

# I made same plots as above, and observed a trend 

# I think i will focus only on GI.2VsGI.2 and GI.1VsGI.1

g <- ggplot(data2, aes(x=aa_difference, y=year_diff)) 
g + geom_bin2d() + 
  theme_bw() + labs(title = "Heatmap-pairwise distances for Sav GII-aa among variants of same genotypes") # to make a Heatmap of 2d bin counts
# This heatmap shows no trend specific trend 

# it seems that SaV G1 has a statist evolutionnary pattern 

#  scatter plot  and label by points


#  scatter plot  and label by points
setwd("C:/Users/Virology Tohoku/Dropbox/manuscript-plos-pathogens/supplements/")

# output
ggplot(data2, aes(aa_difference, year_diff)) +
  geom_point(aes(colour = label), size=3) + 
  xlab("amino acid distance") + 
  ylab("Detection year difference") + 
  font("xlab", size=14)+ 
  font("ylab", size=15)+
  theme_bw(base_size = 14)


ggsave(filename = "evolution pattern genogroup II.tiff", 
       width = 8, height = 5, dpi = 300, units = "in")



# For Sapovirus GIV  
setwd("~/MyRdirectory/SAPO/Sapovirus evolution over time")
data = read_xlsx("SaVaminoacid_diffbigdata_4plotting.xlsx", sheet=1, col_names = T) # to use the # aa differences 

head(data)

data <- data %>% dplyr::rename(specie.1="Species 1", specie.2= "Species 2")

table(grepl("AB436383.1_Sapovirus_GIV_Nagano18-4_JPN_2004", data$specie.2))

table(data$specie.1) # GII X, GII.3, GII.2, GII.1, GII.6, GII.8, GII.3, GII.5, GII.4, 


data<- dplyr::filter(data, grepl('GIV', specie.1))
data<- dplyr::filter(data, grepl('GIV', specie.2)) # keep only GVI sequences. 

table(data$specie.1)# "rename "MG012462.1 Sapovirus/Lima1864 PE GIV.12016" 
data$specie.1 <- ifelse(data$specie.1=="MG012462.1 Sapovirus/Lima1864 PE GIV.12016", "MG012462.1 Sapovirus/Lima1864 PE GIV.1 2016", data$specie.1)
data$specie.2 <- ifelse(data$specie.2=="MG012462.1 Sapovirus/Lima1864 PE GIV.12016", "MG012462.1 Sapovirus/Lima1864 PE GIV.1 2016", data$specie.2)

table(data$specie.1)
table(data$specie.2)# "GIV" ,"GIV.1", 

library(stringr)
data$gen_spe_1 <- word(data$specie.1, -2)
data$gen_spe_2 <- word(data$specie.2, -2)
head(data)

data$gen_spe_2[data$gen_spe_2=="GIV"] <- "GIV.1"
data$gen_spe_1[data$gen_spe_1=="GIV"] <- "GIV.1"
table(data$gen_spe_1) # GIV and GIV.1 
table(data$gen_spe_2)

data$label<- paste(data$gen_spe_1,data$gen_spe_2, sep = "vs")
table(data$label)

data <- data[!grepl("AB436383.1_Sapovirus_GI.1_Nagano18-4_JPN_2004", data$specie.2),]
data <- data[!grepl("AB436383.1_Sapovirus_GI.1_Nagano18-4_JPN_2004", data$specie.1),]

class(data)
library(ggplot2)

ggplot(data, aes(aa_difference)) +geom_density()+ggtitle(label = "distribution # aa differences for Sav GIV ")

# Here, everythings seems ok. OKa et al, 2012 reported that mean +- 3sd of each pairwise distance value for strain is 0-0.127, and for genotype is 0.2-0.37. 
# therefore, this distribution is in accordance with oka's observations.
# now i need to compute one column which represent the difference in isolation year 

# to extract last 4 characters
data$year1 = as.numeric(str_sub(data$specie.1,-4,-1))
data$year2= as.numeric(str_sub(data$specie.2,-4,-1))
head(data)
# rename 
summary(data)
describe(data)
table(data$label)

data$year_diff<- abs(data$year1-data$year2)# compute the difference

hist(data$year_diff) # it is skewed to the left. most samples were isolated during the same period/ Error
# correlation coeficient between aa_differnce and isolation year_diff

head(data)
hist(data$aa_difference) # the distribution of this variable is skewed to the left
hist(data$year_diff) # the distribution of this variable skewed to the left
qqnorm(data$aa_difference, main = "Normal Q-Q Plot",
       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
       plot.it = TRUE, datax = FALSE)
qqline(data$aa_difference, col = 2)

# sincethe data is bimodal, i will select a subset
data2<- filter(data, aa_difference>30)
hist(data2$aa_difference)
# because we are sure that at least one of our variables meets the normality assumption, we can proceed to the correlation test
cor.test(data$aa_difference, data$year_diff, method = "spearman", conf.level = 0.95) 

# data:  data$aa_difference and data$year_diff
# S = 70283000, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#  rho 
# 0.6661697 = Value of the correlation coeficient in the original dataset


cor.test(data2$aa_difference, data2$year_diff, method = "spearman", conf.level = 0.95) 
#  rho 
# 0.8178146   = Value of the correlation coeficient in the subset of the dataset
# p-value = 3.999e-12



g <- ggplot(data, aes(x=aa_difference, y=year_diff)) 
g+ geom_bin2d() + 
  theme_bw() + labs(title = "Evolutionary pattern with time for Sav GIV(heatmap)-# aa differences") # to make a Heatmap of 2d bin counts

#  scatter plot and label by points

ggplot(data, aes(aa_difference, year_diff)) +
  geom_point(aes(colour = label)) + labs(title = "Evolutionary pattern with time for Sav GIV") 

# Fig 1 B 
ggplot(data, aes(aa_difference, year_diff)) +
  geom_point(aes(colour = label), size=4) + 
  xlab("amino acid difference") + 
  ylab("detection year difference") + 
  theme_bw(base_size = 14)


ggsave(filename = "evolution pattern genogroup IV.tiff", 
       width = 8, height = 5, dpi = 300, units = "in")


# For sapovirus GV ------------------------------------------------------------------

data = read_xlsx("SaVaminoacid_diffbigdata_4plotting.xlsx", sheet=1, col_names = T) # to use the # aa differences 


data <- data %>% rename(specie.1="Species 1", specie.2= "Species 2")
table(data$specie.1) # GII X, GII.3, GII.2, GII.1, GII.6, GII.8, GII.3, GII.5, GII.4, 
table(str_extract(data$specie.2,"GV.[0-5]")) # match GII and any2 following character

data$gen_spe_1 = str_extract(data$specie.1,"GV.[0-5]")# to extract the genotypes
data$gen_spe_2 = str_extract(data$specie.2,"GV.[0-5]")# to extract the genotypes

# fix the levels 

data$gen_spe_1[data$gen_spe_1=="GV 1"] <- "GV.1"
data$gen_spe_1[data$gen_spe_1=="GV 2"] <- "GV.2"
data$gen_spe_2[data$gen_spe_2=="GV 1"] <- "GV.1"
data$gen_spe_2[data$gen_spe_2=="GV 2"] <- "GV.2"

table(data$gen_spe_1)

table(is.na(data$gen_spe_2))
table(is.na(data$gen_spe_1))

data = data[!is.na(data$gen_spe_1), ] # to keep all rows of the 2nd col with no missing value
data = data[!is.na(data$gen_spe_2), ] # to keep all rows of the 2nd col with no missing value

data<- data[complete.cases(data), ] # to keep only lines without NAs

data$label<- paste(data$gen_spe_1,data$gen_spe_2, sep = "VS")

data$aa_difference
class(data)


f<- ggplot(data, aes(aa_difference))
f+geom_density()+ggtitle(label = "distribution of pairwise distances for Sav GV-aa") # trimodal distribution, showing variants, genotypes and genogroup distribution

# Here, everythings seems ok. OKa et al, 2012 reported that mean +- 3sd of each pairwise distance value for strain is 0-0.127, and for genotype is 0.2-0.37. 
# therefore, this distribution is in accordance with oka's observations.


# Compute one column which represent the difference in isolation year 

# to extract last 4 characters
data$year1 = as.numeric(str_sub(data$specie.1,-4,-1))
data$year2= as.numeric(str_sub(data$specie.2,-4,-1))

summary(data)
describe(data)
table(data$label)

data$year_diff<- abs(data$year1-data$year2)# compute the difference

ggplot(data, aes(year_diff))+geom_histogram()+ggtitle(label="histogram of year of isolation GV")

# Homogenotypes

data2<- data %>% 
  filter(label %in% c("GV.1VSGV.1", "GV.2VSGV.2"))

hist(data$year_diff) # it is skewed to the left. most samples were isolated during the same period/ Error

# it seems that SaV G2 has also has static evolutionnary pattern 

ggplot(data2, aes(aa_difference, year_diff)) +
  geom_point(aes(colour = label), size=3) + 
  xlab("Amino acid difference")+ ylab("Detection year difference")+ 
  theme_bw(base_size = 15)


cor.test(data2$aa_difference, data2$year_diff, method = "spearman", conf.level = 0.95) 
# the correlation is not significantly different from zero

ggsave(filename = "evolution pattern genogroup5.tiff", 
       width = 8, height = 5, dpi = 300, units = "in")












# Figure 2 --------------------------


# Plot of GI.1 and GI.2 heatmaps

# Heatmap for GI.1 only 
data1.1<- data2 %>% filter(label %in% c("GI.1VsGI.1"))# to make sure the labels are correct
table(data1.1$gen_spe_1,data1.2$gen_spe_2)
# Heatmap for GI.2 only 
data1.2<- data2 %>% filter(label %in% c("GI.2VsGI.2"))# to make sure the labels are correct
table(data1.2$gen_spe_1,data1.2$gen_spe_2)
require(gridExtra)

g1 <- ggplot(data1.1, aes(x=aa_difference, y=year_diff))+ 
  geom_bin2d() + 
  theme_bw(base_size = 15) + 
  xlab("Amino acid differences") + ylab("Detection year difference") + 
  scale_fill_gradient('counts', low = "blue", high = "red")+
  ggtitle("Sapovirus GI.1")+
  theme(plot.title=element_text(size=15, vjust=2, family="Times New Roman"),  
        axis.title.x=element_text(vjust=-1, size=15, family="Times New Roman"),
        axis.title.y=element_text(vjust=-0.25, size=15, family="Times New Roman"), 
        legend.text=element_text(size=15, family="Times New Roman"))

g2 <- ggplot(data1.2, aes(x=aa_difference, y=year_diff))+ 
  geom_bin2d() + 
  theme_bw(base_size = 15) + 
  xlab("Amino acid differences") + ylab("Detection year difference") + 
  scale_fill_gradient('counts', low = "blue", high = "red")+
  ggtitle("Sapovirus GI.2")+
  theme(plot.title=element_text(size=15, vjust=2, family="Times New Roman"),  
        axis.title.x=element_text(vjust=-1, size=15, family="Times New Roman"),
        axis.title.y=element_text(vjust=-0.25, size=15, family="Times New Roman"), 
        legend.text=element_text(size=15, family="Times New Roman"))

library(patchwork)

g1/g2 + plot_layout()+
  plot_annotation(tag_levels = 'A')

ggsave(filename = "sapovirus_gb_genetic_divers.tiff", 
       width = 10, height = 13, dpi = 350, units = "in")



## Compute the correlation coeficient between aa_differnce and isolation year_diff

data3<- data2 %>% dplyr::filter(label=="GI.1VsGI.1")
data4<- data2 %>% dplyr::filter(label=="GI.2VsGI.2")

hist(data3$aa_difference) # the distribution of this variable is skewed to the left
hist(data3$year_diff) # the distribution of this variable skewed to the left
qqnorm(data3$aa_difference, main = "Normal Q-Q Plot",
       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
       plot.it = TRUE, datax = FALSE)
qqline(data3$aa_difference, col = 2)
# because we are sure that at least one of our variables meets the normality assumption, we can proceed to the correlation test

cor.test(data3$aa_difference, data3$year_diff, method = "spearman", conf.level = 0.95) 

cor.test(data4$aa_difference, data4$year_diff, method = "spearman", conf.level = 0.95)

# 1) for SaV GI.1, r_SaV GI.2 = 0.36 (p-value < 2.2e-16)
# 2) for Sav GI.2 ,  # r_SaV GI.2 = 0.86 (p-value < 2.2e-16)

