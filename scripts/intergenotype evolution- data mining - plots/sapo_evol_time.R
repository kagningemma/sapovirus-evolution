
setwd("~/MyRdirectory/SAPO/Sapovirus evolution over time")
library(ggplot2)
library(plotly)
library(stringr)
library(readxl)
library(psych)
library(tibble)

data = read_xls("pairwisebigdata_4plottingxls.xls", sheet=1, col_names = T)
# this dataset was obtained by using poisson distribution based pairwise amino acid comparison of the dataset containing all 
# human sapovirus genotypes in the file allSaVdata_for pairwise.fas. An excell file was generated and uploaded here. 
head(data)
glimpse(data)

# For SaV GI --------------------------------------------------------------

# Genogroup wise analyses

library(tidyverse)
library(stringi)


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
table(is.na(data$gen_spe_1))
#FALSE  TRUE 
#5405 22561 

#data <- data %>% filter(specie.1 !="NA") %>% filter(specie.2 !="NA")# not working
#data <- data %>% filter(!is.na(specie.1))# not working

data<- data[complete.cases(data), ] # to keep only lines without NAs

data$label<- paste(data$gen_spe_1,data$gen_spe_2, sep = "Vs")
data$Dist

# there is a group of pdist >0.5, why ? i need to checkk it out

data%>% filter(Dist>0.5) # to see all weird pairwie comparison with
# sapovirus GQ261222 in bangladesh is a recombinant sequence, thus i will remove it.
outlier <- data %>% filter(Dist>0.5)

data <- data[!grepl("GQ261222.1 Sapovirus BD/697 BGD GI 2005", data$specie.1),] 
data <- data[!grepl("GQ261222.1 Sapovirus BD/697 BGD GI 2005", data$specie.2),]


ggplot(data, aes(Dist))+ 
  geom_density()+ggtitle(label = "distribution of pairwise distances for Sav G1-aa") # trimodal distribution, showing variants, genotypes and genogroup distribution


# now i need to compute one column which represent the difference in isolation year 

# to extract last 4 characters
data$year1 = as.numeric(str_sub(data$specie.1,-4,-1))
data$year2= as.numeric(str_sub(data$specie.2,-4,-1))
hist(data$year1)
hist(data$year2)

summary(data)
describe(data)
table(data$label)
data$year1
data$year_diff<- abs(data$year1-data$year2) # compute the difference
ggplot(data, aes(year1))+geom_histogram()
ggplot(data, aes(year2))+geom_histogram()
ggplot(data, aes(year_diff)) + geom_histogram() + ggtitle( label="histogram of isolation years difference of SaV GI")
hist(data$year_diff) # it is skewed to the left. most samples were isolated during the same period/ Error

g <- ggplot(data, aes(x=Dist, y=year_diff)) 
g+ geom_bin2d() + 
  theme_bw() + labs(title = "Heatmap-Distribution of pairwise distances for Sav GI-aa") # to make a Heatmap of 2d bin counts

# it seems that SaV G1 has a statist evolutionnary pattern 

#  scatter plot  and label by points
setwd("C:/Users/Virology Tohoku/Dropbox/manuscript-plos-pathogens/supplements/")

evg1 = ggplot(data, aes(Dist, year_diff)) +
  geom_point(aes(colour = label),size=2) + xlab("Distance") + 
  ylab("Detection year difference") +  theme_bw()

library(patchwork)  
evg1 + plot_annotation(tag_levels = 'A')
ggsave(filename = "evolution pattern gengroup1.tiff", 
       width = 22, height = 13, dpi = 350, units = "cm")

# Genotype-wise analyses

data<- data %>% mutate(homo_hetero  =  ifelse(data$gen_spe_1==data$gen_spe_2, "homogenotype", "heterogenotype")) 
data %>% select(gen_spe_1, gen_spe_2, homo_hetero) %>% filter(homo_hetero=="homogenotype") %>% table() 

class(data)
head(data)
table(data$gen_spe_1)


ggplot(data, aes(Dist, year_diff))+ geom_point(aes(colour=homo_hetero, shape=homo_hetero)) + labs(title = "scatterplot Sav GI-aa-Homogenotypes-heterogenotypes")
ggplotly(p)


# most homogenotypes have a pdist<0.1 (99% homology) and there is a little trend if i only focus on variants. 
 

# below, I will only focus on those with pdist<0.1 (that variants resulting from pairwise comparison)
 
 
 data2 <- data %>% filter(Dist<0.1)


 # GI.1VsGI.1 GI.1VsGI.2 GI.2VsGI.1 GI.2VsGI.2 GI.3VsGI.3 GI.5VsGI.5 GI.6VsGI.6 GI.7VsGI.7 
 # 780         21         19        666         21         10         15          3 
 
 # I made same plots as above, and observed a trend 
 
 # I think i will focus only on GI.2VsGI.2 and GI.1VsGI.1
 
 g <- ggplot(data2, aes(x=Dist, y=year_diff)) 
 g + geom_bin2d() + 
   theme_bw() + labs(title = "Heatmap-pairwise distances for Sav GI-aa among sequences of the same genotypes", 
                     subtitle = "There is a trend but it's not well visible") # to make a Heatmap of 2d bin counts
 
 # it seems that SaV G1 has a statist evolutionnary pattern 
 
 #  scatter plot  and label by points
 
 ggplot(data2, aes(Dist, year_diff)) +
     geom_point(aes(colour = label)) + labs(title = "scatterplot of pairwise distances for Sav GI-aa among sequences of same genotypes")
 # the visible trend is not well visible
 
 ggplot(data2, aes(Dist, year_diff)) +
   geom_point(aes(shape = label)) + labs(title = "scatterplot of pairwise distances for Sav GI-aa among sequences of same genotypes")
 
 
 p= ggplot(data2, aes(Dist, year_diff))+ geom_point(aes(colour=homo_hetero, shape=homo_hetero)) + 
   labs(title = "scatterplot Sav GI-aa- among sequences of same genogytpes")
 
 ggplotly(p)
 
 
#Since i observed a trend between G1.1 and G1.2 homogenotype comparisons, i will focus only on that
 
 data2 <- data2 %>% filter(label %in% c("GI.2VsGI.2", "GI.1VsGI.1"))
 
# I think i will focus only on GI.2VsGI.2 and GI.1VsGI.1
 
g <- ggplot(data2, aes(x=Dist, y=year_diff)) 
g + geom_bin2d() + 
   theme_bw() + labs(title = "Heatmap-pairwise distances for Sav GI-aa among sequences of same genotype") # to make a Heatmap of 2d bin counts

 # same plot as above - try to make plot for calicivirus conference
 # black and white plot

 # we can see that there seem to be an evolutionary pattern, but not very informative
 
 #  scatter plot  and label by points
head(data2)
 ggplot(data2, aes(Dist, year_diff)) +
   geom_point(aes(colour = label)) + labs(title = "scatterplot of pairwise distances for Sav GI.1 and GI.2 genotypes")
 
 # we can see that GI.2 genotype changes over times more than GI.1 genotypes. 
 # based on this, i hypothesize that sapovirus GI.2 has more evolutionary changes over time compared to sapovirus GI.1
 
 ggplot(data2, aes(Dist, year_diff)) +
   geom_point(aes(shape = label)) + labs(title = "scatterplot of pairwise distances for Sav GI.1 and GI.2 genotypes")
 
 
 # same plot as above
 
 
 # I would like to see this heatmap for GI.1 only 
 data1.1<- data2 %>% filter(label %in% c("GI.1VsGI.1"))
 table(data1.2$gen_spe_1,data1.1$gen_spe_2) # to make sure the labels are correct
 g <- ggplot(data1.1, aes(x=Dist, y=year_diff)) 
 g + geom_bin2d() + 
   theme_bw() + labs(title = "Heatmap-pairwise Amino Acid distances between Sapovirus GI.1 sequences", 
                     subtitle = "static evolutionary pattern", caption = "source= genbank sequences") + 
   xlab("poisson adjusted distance") + ylab("Detection year difference")# to make a Heatmap of 2d bin counts
 
 # Heatmap for GI.2 only 
 data1.2<- data2 %>% filter(label %in% c("GI.2VsGI.2"))# to make sure the labels are correct
 table(data1.2$gen_spe_1,data1.2$gen_spe_2)
 g <- ggplot(data1.2, aes(x=Dist, y=year_diff)) 
 g + geom_bin2d() + 
   theme_bw() + labs(title = "Heatmap-pairwise Amino Acid distances between Sapovirus GI.2 sequences", 
                     subtitle = "'Evolving' evolutionary pattern", caption = "source= genbank sequences") + 
   xlab("poisson adjusted distance") + ylab("Detection year difference")# to make a Heatmap of 2d bin counts
 
 # Because of the observation above, i would like to make a separate alignement containing only G1.1 and GII.2. 
 # and make the above plots using x as number of amino acid differences. 
 
 library(stringr)
 library(stringi)
 allfasta=readDNAStringSet(filepath = "C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/allSaV_capsid_data_for pairwise_ok.fas")
 grep(names(allfasta), "GI.1")
 
 stri_extract(names(allfasta))
 stri_match(names(allfasta), "GI.1")
 
 allfasta[grep(c("GI.1"), names(allfasta)), ]
 allfasta[grep(c("GV.2"), names(allfasta)), ]
 
 mtcars[grep("Merc", rownames(mtcars)), ] # to select all car names starting with "merc"

 allfasta2<- allfasta[grep(c("GI.1|GI.2"), names(allfasta)), ] # to select only fasta sequences containing GI.1 or GI.2

 out23a <- tempfile()
 writeXStringSet(allfasta2, out23a)# export the fasta file
 out23a # go to the location and copy paste the temp file into your working directory. rename the extension of file into .fas
 
 allfasta3<- allfasta[grep("GI.2", names(allfasta)), ] 
 allfasta4<- allfasta[grep("GI.1", names(allfasta)), ] 
 g1.1 <- tempfile()
 g1.2<- tempfile()
 writeXStringSet(allfasta4,  g1.1)
 writeXStringSet(allfasta3,  g1.2)
 g1.1 # to see fasta file containing only G1.1 sequences_very fast method, rather than doing manually
 g1.2
 
 # export as fasta file with name "G1.1orGI.2_msa_forplotting.fas" 
 # open the file into MEGA, exclude the bangladesh sequence (GQ261222.1_Sapovirus_BD/697_BGD_GI_2005) make alignment, make protein p-distaces and re-import for plotting
 # align amino acid using muscle, find the best substitution model, (JTT+Gamma with 5 rate cat)
 # i calculated the within group distance between GI.1 and GI.2. 
 # mean number of aa differences is 10.83 for GI.2 and 4.68 for GI.1
 # this confirming my hypothesis that GI.2 is more evolves faster than GI.1
 # now let us plot that using 
 
 
 

# Plotting AA distances for GI.1 and GI.2---------------------------------------------------


 
 data = read_xlsx("sapoGI.2andGI.1numberofaadiff.xlsx", sheet=1, col_names = T)
 
 colnames(data)[1] <- "specie.1"
 colnames(data)[2] <- "specie.2"
 colnames(data)[3] <- "aa_difference"
 
 data <- data[!grepl("Chanthaburi-74 Thailand AY646854 GI 2002", data$specie.1),] # It was mislabelled as GI instead of GII, and there was already a duplicate
 data <- data[!grepl("Chanthaburi-74 Thailand AY646854 GI 2002", data$specie.2),] # to exclude that one
 data <- data[!grepl("GQ261222.1 Sapovirus BD/697 BGD GI 2005", data$specie.2),]
 table(grepl("Chanthaburi-74 Thailand AY646854 GI 2002", data$specie.1))
 
 
 data$gen_spe_1 = str_extract(data$specie.1, "[G][I|1][.| ][0-9 | I-IV]")
 data$gen_spe_2 = str_extract(data$specie.2, "[G][I|1][.| ][0-9 | I-IV]")# to extract the genotypes
 
 
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
 table(is.na(data$gen_spe_1))
 
 
 
 data$label<- paste(data$gen_spe_1,data$gen_spe_2, sep = "Vs")
 data$label2<- paste(data$specie.1,data$specie.2, sep = "AND")
 
 ggplot(data, aes(Dist))+ 
   geom_area(stat="bin")+ggtitle(label = "distribution of #aa diff btw GI.1 andGI.2") # bimodal distribution, showing variants, genotypes and genogroup distribution
 
 ggplot(data, aes(aa_difference))+ 
   geom_area(stat="bin")+ggtitle(label = "distribution of #aa diff btw GI.1 andGI.2") # bimodal distribution, showing variants, genotypes and genogroup distribution
 
 # now i need to compute one column which represent the difference in isolation year 
 
 # to extract last 4 characters
 data$year1 = as.numeric(str_sub(data$specie.1,-4,-1))
 data$year2= as.numeric(str_sub(data$specie.2,-4,-1))
 
 summary(data)
 describe(data)
 table(data$label)
 
 data$year_diff<- abs(data$year1-data$year2)# compute the difference
 ggplot(data, aes(year_diff)) + geom_histogram() + ggtitle( label="histogram of isolation years difference of GI.1 and GI.2")
 hist(data$year_diff) # it is skewed to the left. most samples were isolated during the same period/ Error
 
 g <- ggplot(data, aes(x=aa_difference, y=year_diff)) 
 g+ geom_bin2d() + 
   theme_bw() + labs(title = "Heatmap-Distribution of #aa diff btw GI.1 and GI.2")
 
 # unline those of GI.6 of dr parra, we have only 2 lines, representing inter-cluster comparison and intra-cluster comparisons.
 
 #  scatter plot  and label by points
 
 ggplot(data, aes(aa_difference, year_diff)) +
   geom_point(aes(colour = label)) + labs(title = "Labelled scatterplot of #aa diff btw GI.1 and GI.2")
 
 
 # # Genotype-wise analyses
 # let us focus on same genotype compasisons
 data$label
data2<- data %>% 
  filter(label %in% c("GI.1VsGI.1", "GI.2VsGI.2", "GI.3VsGI.3", "GI.4VsGI.4", 
                       "GI.5VsGI.5", "GI.6VsGI.6", "GI.7VsGI.7", "GI.8VsGI.8"))
 
ggplot(data2, aes(aa_difference, year_diff)) +
  geom_point(aes(colour = label)) + 
  labs(title = "Labelled scatterplot of aa diff btw GI.1 and GI.2")+
  xlab("amino acid difference")+
  ylab("detection year difference")+ 
  theme_bw()

 
 
 data2<- data %>% filter(label %in% c("GI.1VsGI.1", "GI.2VsGI.2" ))
 table(data2$label)
 
 g <- ggplot(data2, aes(x=aa_difference, y=year_diff)) 
 g + geom_bin2d() + 
   theme_bw() + labs(title = "Heatmap-# of aa differences for Sav GI.1 VS GI.2 ") # to make a Heatmap of 2d bin counts
 
 # it seems that SaV G1 has a static evolutionnary pattern 
# Calicivurus conference plot ---------------------------------------------

 
 #  scatter plot  and label by points
data2 %>% count(label) %>% knitr::kable()
head(data2$label)
data2$Genotype <- ifelse(data2$label=="GI.1VsGI.1", "SaV-GI.1", "SaV-GI.2")
 
 ggplot(data2, aes(x=year_diff, y=aa_difference)) + theme_bw() +
   geom_point(aes(shape = Genotype), size=6) + 
   ggtitle("Evolutionary patterns of sapovirus GI.1 and GI.2") + 
   ylab("Number of amino acid difference") + 
   xlab("Detection year difference (years)")+
   scale_shape_manual(values = c(3,16))+
   theme(legend.title = element_blank(), 
         axis.text = element_text(size = 20),
         plot.title = element_text(color="black", size=19),
         axis.title.x = element_text(color="black", size=20),
         axis.title.y = element_text(color="black", size=20),
         legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"),
         legend.justification=c(1,0), 
         legend.position=c(1,0.85), 
         legend.text=element_text(size=20))

  
 

data2 %>% filter(Dist>90) 
 
 # this plot shows that, among all sav GI.1 variatns, amino acid sequence in most of detected variants do not change a lot despites a 30 years span
 # while among GI.2 variants, there is a positive and steep slope describing the relation between aa changes and time.
 # in some GI.2 variants isolated the 1year later, there were more than 100 amino acid differences. 
 
 ggplot(data2, aes(aa_difference, year_diff)) +
   geom_point(aes(shape = label)) + labs(title = "scatterplot of pairwise distances for Sav GI.1 variants VS GI.2 variants")
 
 # it is difficult to choose colour for calicivirus conference colours
# http://research.stowers.org/mcm/efg/R/Color/Chart/
 
 
 # I would like to see this heatmap for GI.1 only 
 data1.1<- data2 %>% filter(label %in% c("GI.1VsGI.1"))
 table(data1.2$gen_spe_1,data1.1$gen_spe_2) # to make sure the labels are correct
 g <- ggplot(data1.1, aes(x=aa_difference, y=year_diff)) 
 g + geom_bin2d() + 
   theme_bw() + labs(title = "Heatmap of pairwise amino acid distances between sapovirus GI.1 sequences", 
                     subtitle = "static evolutionary pattern", caption = "Data source = GenBank") + 
   xlab("number of amino acid difference") + ylab("Detection year difference")# to make a Heatmap of 2d bin counts
 

# plot of GI.1 and GI.2 heatmaps for publication --------------------------

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
 
 setwd("C:/Users/Virology Tohoku/Dropbox/manuscript-plos-pathogens/supplements")
 
 library(patchwork)
 
g1/g2 + plot_layout()+
  plot_annotation(tag_levels = 'A')
 
 setwd("C:/Users/Virology Tohoku/Dropbox/manuscript-plos-pathogens/supplements")
 
 ggsave(filename = "sapovirus_gb_genetic_divers.tiff", 
        width = 10, height = 13, dpi = 350, units = "in")
 



 
 # black and white plot
 
 g + geom_bin2d() + 
   theme_bw() + labs(title = "Heatmap-pairwise amino acid differences among sapovirus GI.2 sequences", 
                     subtitle = "'Evolving' evolutionary pattern", caption = "data source= genbank") + 
   xlab("number of amino acid difference") + ylab("Detection year difference") + 
   scale_fill_gradient('counts', low = "blue", high = "red") # to change the fill color palette
 
 colours() # to see all colours of ggplot
 
 write.csv(data1.2, "sapovirus GI.2 evolutionary pattern with time.csv") #  for future use, ill save the csv file to generate this plot easier
 write.csv(data1.1, "sapovirus GI.1 evolutionary pattern with time.csv")

# Correlation coeficient G1.1-G1.2 ----------------------------------------

 
 ## correlation coeficient between aa_differnce and isolation year_diff
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
 

 
 
 
 
 
 
 
 
 
# for Sapovirus G2 ---------------------------------------------------
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
library(ggplot2)

f<- ggplot(data, aes(aa_difference))
f+geom_density()+ggtitle(label = "Distribution of pairwise distances for Sav GII-aa") # bimodal distribution, showing variants, genotypes and genogroup distribution
# Here, everythings seems ok. OKa et al, 2012 reported that mean +- 3sd of each pairwise distance value for strain is 0-0.127, and for genotype is 0.2-0.37. 
# therefore, this distribution is in accordance with oka's observations.


# now i need to compute one column which represent the difference in isolation year 

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


g <- ggplot(data, aes(x=aa_difference, y=year_diff)) 
g+ geom_bin2d() + 
  theme_bw() + labs(title = "Evolutionary pattern of Sav GII with time (heatmap)-# aa differences") # to make a Heatmap of 2d bin counts

# it seems that SaV G2 has also has static evolutionnary pattern 

#  scatter plot and label by points

ggplot(data, aes(aa_difference, year_diff)) +
  geom_point(aes(colour = label)) + 
  labs(title = "Evolutionary pattern of sapovirus Genogroup II")



# GII_ HOMOGENOTyPES=HETEROGENOTYPEs --------------------------------------


# # Genotype-wise analyses

table(data$specie.1)
table(data$gen_spe_1)
data <- data %>% mutate(homo_hetero  =  ifelse(data$gen_spe_1==data$gen_spe_2, "homogenotype", "heterogenotype")) # creating new variable
data %>% select(homo_hetero) %>% table() 

class(data)
library(plotly)
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
  ylab("Detection year difference") +  font("xlab", size=14)+ font("ylab", size=15)+
  theme_bw(base_size = 14)

 
ggsave(filename = "evolution pattern genogroup II.tiff", 
       width = 8, height = 5, dpi = 300, units = "in")







# the visible trend is not well visible

ggplot(data2, aes(aa_difference, year_diff)) +
  geom_point(aes(shape = label)) + labs(title = "scatterplot of pairwise distances for Sav GI-aa among variants of same genotypes")



p= ggplot(data2, aes(Dist, year_diff))+ geom_point(aes(colour=homo_hetero, shape=homo_hetero)) + labs(title = "scatterplot Sav GI-aa- among variants of same genogytpes")

ggplotly(p)




# Since i observed a trend between G1 and G2 homogenotype comparisons, i will focus only on that

data2<- data2 %>% filter(label %in% c("GI.2VsGI.2", "GI.1VsGI.1"))

# I think i will focus only on GI.2VsGI.2 and GI.1VsGI.1

g <- ggplot(data2, aes(x=Dist, y=year_diff)) 
g + geom_bin2d() + 
  theme_bw() + labs(title = "Heatmap-pairwise distances for Sav GI-aa among variants of same genotypes") # to make a Heatmap of 2d bin counts

# we can see that there seem to be an evolutionary pattern, but not very informative
























# For GIV -----------------------------------------------------------------

library(dplyr)
library(readxl)

data = read_xls("pairwisebigdata_4plottingxls.xls", sheet=1, col_names = T) # to use the pairwise distances, Dist
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
#sentences <- c("Jane saw a cat", "Jane sat down")
#word(sentences, 1) # to extract the 1st word of each sentence / string
#word(sentences, 2)
#word(sentences, -2)

data$gen_spe_1 <- word(data$specie.1, -2)
data$gen_spe_2 <- word(data$specie.2, -2)
head(data)

data$gen_spe_2[data$gen_spe_2=="GIV"] <- "GIV.1"
data$gen_spe_1[data$gen_spe_1=="GIV"] <- "GIV.1"
table(data$gen_spe_1) # GIV and GIV.1 
table(data$gen_spe_2)

#data$gen_spe_1 = str_extract(data$specie.1,"GIV.[0-5]")# to extract the genotypes
#data$gen_spe_2 = str_extract(data$specie.2,"GIV.[0-5]")# to extract the genotypes


# data = data[!is.na(data$gen_spe_1), ] # to keep all rows of the 2nd col with no missing value
# data = data[!is.na(data$gen_spe_2), ] # to keep all rows of the 2nd col with no missing value

# data<- data[complete.cases(data), ] # to keep only lines without NAs

data$label<- paste(data$gen_spe_1,data$gen_spe_2, sep = "vs")
table(data$label)
#data <- data[!grepl("GQ261222.1 Sapovirus BD/697 BGD GI 2005", data$specie.1),] 
#data <- data[!grepl("GQ261222.1 Sapovirus BD/697 BGD GI 2005", data$specie.2),]

# outlier = data%>% filter(data$Dist>0.4)


data <- data[!grepl("AB436383.1_Sapovirus_GI.1_Nagano18-4_JPN_2004", data$specie.2),]
data <- data[!grepl("AB436383.1_Sapovirus_GI.1_Nagano18-4_JPN_2004", data$specie.1),]

# how many aa diff needed to have a new genotye
# 0,75 ----> 40 AA
# 0.13 ---? ?      # 0.13 distance is cutoff for defining new genotype
# ? = 70 Amino acid differences





class(data)
library(ggplot2)

ggplot(data, aes(Dist)) +geom_density()+ggtitle(label = "distribution aa p dist for Sav GIV ") # bi modail distribution, showing variants, genotypes and genogroup distribution

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
# 0.6661697 

cor.test(data2$aa_difference, data2$year_diff, method = "spearman", conf.level = 0.95) 
#  rho 
# 0.8178146
# p-value = 3.999e-12

#
# 
g <- ggplot(data, aes(x=Dist, y=year_diff)) 
g+ geom_bin2d() + 
  theme_bw() + labs(title = "Evolutionary pattern with time for Sav GIV(heatmap)-aa-pairwise dist") # to make a Heatmap of 2d bin counts

g <- ggplot(data, aes(x=aa_difference, y=year_diff)) 
g+ geom_bin2d() + 
  theme_bw() + labs(title = "Evolutionary pattern with time for Sav GIV(heatmap)-# aa differences") # to make a Heatmap of 2d bin counts


g <- ggplot(data, aes(x=Dist, y=year_diff)) 
g+ geom_bin2d() + 
  theme_bw() + labs(title = "Evolutionary pattern with time for Sav GIV(heatmap)-# aa differences") # to make a Heatmap of 2d bin counts


# All sapovirus GIV sequences are may have an evolving evolutionary pattern, but sample size is too low 

# Taken all together, it seems that SaV G4 has also has static evolutionnary pattern 

#  scatter plot and label by points

ggplot(data, aes(aa_difference, year_diff)) +
  geom_point(aes(colour = label)) + labs(title = "Evolutionary pattern with time for Sav GIV") 
evg4= ggplot(data, aes(Dist, year_diff)) +
  geom_point(aes(colour = label)) + 
  xlab("Distance") + 
  ylab("detection year difference") + 
  theme_bw(base_size = 15)

# if i imported the amino acid file
ggplot(data, aes(aa_difference, year_diff)) +
  geom_point(aes(colour = label)) + 
  xlab("amino acid difference") + 
  ylab("detection year difference") + 
  theme_bw(base_size = 14) 
setwd("C:/Users/Virology Tohoku/Dropbox/manuscript-plos-pathogens/supplements")

ggsave(filename = "evolution pattern genogroup IV.tiff", 
       width = 8, height = 5, dpi = 300, units = "in")


ggplot(data, aes(aa_difference, year_diff)) +
  geom_point(aes(colour = label)) + labs(title = "Evolutionary pattern with time for Sav GIV", 
                                         subtitle = "in this study p-distances were obtained using poisson method-for AA alignmenent")

# i want to label those points with Dist >0.5, to have a closer look at those sequences

ggplot(data, aes(aa_difference, year_diff)) +
  geom_point(aes(colour = label)) + labs(title = "Evolutionary pattern with time for Sav GIV")+ 
  geom_text(aes(label=label),hjust=0, vjust=0)



ggplot(data, aes(Dist, year_diff)) +
  geom_point(aes(colour = label)) + labs(title = "Evolutionary pattern with time for Sav GIV")+ 
  geom_text(aes(label=ifelse(Dist>0.5,as.character(label),'')),hjust=0.1,vjust=0)

# TAF: use alignment containing only GIV sequences
# Check all the labels again



table(data$gen_spe_1) # there are only 2 categories, GVI.1. Let me look at the variants

export<- data %>% dplyr::filter(Dist>0.5)
write.csv(export, "SaV GIV_check_recheck_classificatio_pdist.csv")














# For GV ------------------------------------------------------------------

library(readxl)

data = read_xls("pairwisebigdata_4plottingxls.xls", sheet=1, col_names = T)
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
data$Dist
data$aa_difference
class(data)
library(ggplot2)

f<- ggplot(data, aes(aa_difference))
f<- ggplot(data, aes(Dist))
f+geom_density()+ggtitle(label = "distribution of pairwise distances for Sav GV-aa") # trimodal distribution, showing variants, genotypes and genogroup distribution

# Here, everythings seems ok. OKa et al, 2012 reported that mean +- 3sd of each pairwise distance value for strain is 0-0.127, and for genotype is 0.2-0.37. 
# therefore, this distribution is in accordance with oka's observations.


# now i need to compute one column which represent the difference in isolation year 

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


table(data2$label)

hist(data$year_diff) # it is skewed to the left. most samples were isolated during the same period/ Error

g <- ggplot(data2, aes(x=Dist, y=year_diff)) 
g+ geom_bin2d() + 
  theme_bw() + labs(title = "Evolutionary pattern with time for Sav GV(heatmap)-aa") # to make a Heatmap of 2d bin counts

# it seems that SaV G2 has also has static evolutionnary pattern 

#  scatter plot  and label by points
setwd("C:/Users/Virology Tohoku/Dropbox/manuscript-plos-pathogens/supplements/")

ggplot(data2, aes(Dist, year_diff)) +
  geom_point(aes(colour = label), size=3) + 
  xlab("Distance")+ ylab("detection year difference")+theme_bw()

ggplot(data2, aes(aa_difference, year_diff)) +
  geom_point(aes(colour = label), size=3) + 
  xlab("amino acid difference")+ ylab("detection year difference")+ 
  theme_bw(base_size = 15)

evg5= ggplot(data2, aes(Dist, year_diff)) +
  geom_point(aes(colour = label), size=3) + 
  xlab("Distance")+ ylab("detection year difference")+theme_bw()
  
ggsave(filename = "evolution pattern genogroup5.tiff", 
       width = 8, height = 5, dpi = 300, units = "in")


library(patchwork)
evg1/evg2+evg4+evg5

ggsave(filename = "evolution pattern genogroup4.tiff", 
       width = 8, height = 5, dpi = 300, units = "in")



# Sapovirus Shedding measuremnts ------------------------------------------



# Sapovirus shedding - recomputing - emma ---------------------------------
# samples selecltion for NGS
# selection of samples for NGS
# we need samples from same patients with subsequent samples positive after 7 days and before 30 days. 
library(readxl)
library(tidyverse)
library(lubridate)

dataa = read_xlsx("SAPO_POS_GERARDO_TG1_2_mdfd_Emmanuel's changes 2019-04-11.xlsx", 
                 sheet=5, col_names = T)
data
colnames(data)
glimpse(data)
table(data$date_col...26==data$date_col...32) 
which(data$date_col...26!=data$date_col...32)
data %>% select(date_col...32, date_col...26) %>% slice(131:145) # PPPP ask mayuko sensei ????

data2<- data %>% select(sample_id1, sample_id2, 
                        per_id...24, date_col...26,initials, genotyping_T_1,qPCR_CT_1st_) %>% 
  group_by(per_id...24) %>% 
  filter(n() > 1) %>% # remove per-id that are unique (poeple without repeated infection)
  ungroup()

data3<- data2 %>% group_by(per_id...24) %>% 
  mutate(date_difference = base::difftime(date_col...26, first(date_col...26), units = "days")) %>% 
  ungroup()
data3 <- data3%>%  filter(genotyping_T_1=="GI.2"| genotyping_T_1=="GI.1")

data3 <- data2 %>% group_by(per_id...24)  %>%
  mutate(date_difference = difftime(date_col...26, first(date_col...26), units = "days")) %>% 
  ungroup()

data3 %>% group_by(per_id...24) %>% 
  filter(date_difference-first(date_difference) >=5)  %>% # difference between each date and the index date of the group
  # should be >5
  ungroup()
# to see subsequent samples which are >5 days after the index sample


# Sapovirus episode counting
# A day of diarrhea was defined by the presence of≥3 liquid or semiliquid stools in 24 hours.
# For infants younger than 2 months,the definition was based on the mother’s or caretaker’s assess-ment that the child had diarrhea [15]. 
# An episode ended when the child had 2 consecutive days without diarrhea.



# Re-infection matrix for GI and GII
# 2nd infection is when same persid is infected with same genotype after 30 days. 

#PPPPPP PROBLEm 2: :there are not all genotypes in the new dataset. i need them to make re-infection matrix

data2 %>% 
  data2 %>% group_by(per_id...24) %>% 
  mutate(date_difference = difftime(date_col...26, first(date_col...26), units = "days")) %>% 
  ungroup()








# AIM: Count the shedding duration spare
setwd("C:/Users/Virology Tohoku/Desktop/My desktop recovery/mes travaux de recherche/SAPOVIRUS/Sapovirus Shedding/sapovirussheddingsamplecandidates version 1")
library(ggplot2)
library(plotly)
library(tidyverse)
library(psych)
library(tibble)
library(Hmisc)
data = readxl::read_xlsx("2_SAPO2016_shedding_candidates_050419ms.xlsx", sheet=1, col_names = T)
glimpse(data)
table(data$sav)

data %>% group_by(cod_per) %>% 
  arrange(cod_per, fec_mues) %>% 
  mutate(num=1:n())


# cound shedding duration and compare GI.1 vs GI.2

data %>% filter(sav==1) %>% group_by(cod_per) %>% 
  arrange(cod_per, fec_mues) %>% 
  mutate(num=1:n()) %>% 
  mutate(sav1=)

 
 
 



data2 %>% group_by(date) %>% 
  mutate(date_difference = base::difftime(date, first(date), units = "weeks")) %>% ungroup()
summary(data3)

data3 <- data3%>%  filter(genotyping_T_1=="GI.2"| genotyping_T_1=="GI.1")

data3 <- data2 %>% group_by(per_id...24)  %>%
  mutate(date_difference = difftime(date_col...26, first(date_col...26), units = "days")) %>% 
  ungroup()

data3 %>% group_by(per_id...24) %>% 
  filter(date_difference-first(date_difference) >=5)  %>% # difference between each date and the index date of the group
  # should be >5
  ungroup()
# to see subsequent samples which are >5 days after the index sample




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
setwd("C:/Users/Virology Tohoku/Dropbox/manuscript-plos-pathogens/supplements")

ggsave(filename = "sapo.cap.pol.evol.rates.tiff", 
       width = 22, height = 15, dpi = 300, units = "cm")

