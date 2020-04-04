

setwd("~/MyRdirectory/SAPO/Sapovirus evolution over time")
library(ggplot2)
library(plotly)
library(stringr)
library(readxl)
library(psych)
library(tibble)
data = read_excel("pairwisebigdata_4plottingxls.xls", sheet=1, col_names = T)
# this dataset was obtained by using poisson distribution based pairwise amino acid comparison of the dataset containing all 
# human sapovirus genotypes in the file allSaVdata_for pairwise.fas. An excell file was generated and uploaded here. 
head(data)
glimpse(data)

# For SaV GI --------------------------------------------------------------
# Genogroup wise analyses

library(tidyverse)
library(stringi)

data <- data %>% rename(specie.1="Species 1", specie.2= "Species 2")
str_extract(data$specie.1,"[G][I-IV|1-5][.| ][0-9 | I-IV]") # to exttact g1, g2 , g4 g5
table(str_extract(data$specie.1,"[G][I-IV|1-5][.| ][0-9 | I-IV]"))
# how to use stringi https://r4ds.had.co.nz/strings.html
str_extract(data$specie.1[1],"[G][I|1][.| ][0-9 | I-IV]") # to select G1 only
table(str_extract(data$specie.1, "[G][I|1][.| ][0-9 | I-IV]"))

data$gen_spe_1 = str_extract(data$specie.1, "[G][I|1][.| ][0-9 | I-IV]")
data$gen_spe_2 = str_extract(data$specie.2, "[G][I|1][.| ][0-9 | I-IV]")# to extract the genotypes

table(data$specie.1) %>% table(data$specie.2)
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

data <- data %>% filter(specie.1 !="NA") %>% filter(specie.2 !="NA")# not working
data <- data %>% filter(!is.na(specie.1))# not working

data<- data[complete.cases(data), ] # to keep only lines without NAs

data$Dist

head(data)

# there is a group of pdist >0.5, why ? i need to checkk it out

data%>% filter(Dist>0.5) # to see all weird pairwie comparison with
# sapovirus GQ261222 in bangladesh is a recombinant sequence, thus i will remove it.

data <- data[!grepl("GQ261222.1 Sapovirus BD/697 BGD GI 2005", data$specie.1),] 
data <- data[!grepl("GQ261222.1 Sapovirus BD/697 BGD GI 2005", data$specie.2),]

# now i need to compute one column which represent the difference in isolation year 

# to extract last 4 characters
data$year1 = as.numeric(str_sub(data$specie.1,-4,-1))
data$year2= as.numeric(str_sub(data$specie.2,-4,-1))
data$label <- paste(data$gen_spe_1,data$gen_spe_2, sep = "Vs")
data$year_diff<- abs(data$year1-data$year2)# compute the difference
head(data)

g1data <- data %>% select(gen_spe_1, gen_spe_2, year_diff, Dist,label) %>% dplyr::mutate(gen_spe_1 = as.character(gen_spe_1)) %>% 
  dplyr::mutate(gen_spe_2 = as.character(gen_spe_2))
g1data <- as.data.frame(g1data)
glimpse(g1data)

# for GII -----------------------------------------------------------------

 setwd("~/MyRdirectory/SAPO/Sapovirus evolution over time")
 
 library(readxl)
 
 data = read_excel("pairwisebigdata_4plottingxls.xls", sheet=1, col_names = T)
 
 data <- data %>% dplyr::rename(specie.1="Species 1", specie.2= "Species 2")
 
 table(data$specie.1) # GII X, GII.3, GII.2, GII.1, GII.6, GII.8, GII.3, GII.5, GII.4, 
 
 table(str_extract(data$specie.1,"GII.."))# match GII and any2 following character
 table(str_extract(data$specie.1,"GII.*"))
  # data2<- data %>% filter(grepl("GII.*", specie.1))  %>%  filter(grepl("GII.*", specie.2)) # to only retain row containing GII.. pattern
 # table(data2$specie.1) # i will need this if i need to focus on some specific groups
  # str_extract(data$specie.1[1],"[G][II|2][.| ][0-9 | I-IV]") # to select G2 only
 # table(str_extract(data$specie.1, "[G][I|1][.| ][0-9 | I-IV]"))
  data$gen_spe_1 = str_extract(data$specie.1,"GII..")
 data$gen_spe_2 = str_extract(data$specie.2,"GII..")# to extract the genotypes
  table(data$gen_spe_1)
 table(data$gen_spe_2)
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
 data$label <- paste(data$gen_spe_1,data$gen_spe_2, sep = "Vs")
 class(data)
 library(ggplot2)
 #Here, everythings seems ok. OKa et al, 2012 reported that mean +- 3sd of each pairwise distance value for strain is 0-0.127, and for genotype is 0.2-0.37. 
 # therefore, this distribution is in accordance with oka's observations.
 # now i need to compute one column which represent the difference in isolation year 
 # to extract last 4 characters
 data$year1 = as.numeric(str_sub(data$specie.1,-4,-1))
 data$year2= as.numeric(str_sub(data$specie.2,-4,-1))
 
 summary(data)
 describe(data)
 table(data$label)
 
 data$year_diff<- abs(data$year1-data$year2)# compute the difference
 
 # now i need to compute one column which represent the difference in isolation year 
 
 # to extract last 4 characters
 data$year1 = as.numeric(str_sub(data$specie.1,-4,-1))
 data$year2= as.numeric(str_sub(data$specie.2,-4,-1))
 
 
 data$year_diff<- abs(data$year1-data$year2)# compute the difference
 head(data)
 
glimpse(data)
glimpse(g1data)
 
 g2data <- data %>% select(gen_spe_1, gen_spe_2, year_diff, Dist,label) %>% dplyr::mutate(gen_spe_1 = as.character(gen_spe_1)) %>% 
   dplyr::mutate(gen_spe_2 = as.character(gen_spe_2))
g2data <- as.data.frame(g2data)
 
glimpse(g2data)



# Merging ----------------------------------------------------------------

g1g2forPCA <- full_join(g1data, g2data)
 
 

# PCA demo ----------------------------------------------------------------

data(USArrests)
library(FactoMineR)
res <- PCA(USArrests)
res$eig # to see eigen values end % of variance explained
dimdesc(res, axes=1)  # show correlation of variables with 1st axis
#$`Dim.1`
#$`Dim.1`$`quanti`
#correlation  p.value
#Assault        0.918 5.76e-21
#Rape           0.856 2.40e-15
#Murder         0.844 1.39e-14
#UrbanPop       0.438 1.46e-03
res$var$coord  # show loadings associated to each axis
#Dim.1  Dim.2  Dim.3   Dim.4
#Murder   0.844 -0.416  0.204  0.2704
#Assault  0.918 -0.187  0.160 -0.3096
#UrbanPop 0.438  0.868  0.226  0.0558
#Rape     0.856  0.166 -0.488  0.0371
FactoMineR::dimdesc(res)
FactoMineR::
  
HCPC(res,nb.clust=-1,consol=FALSE,min=3,max=10,graph=TRUE, iter.max =10) # Performs an agglomerative hierarchical clustering on results from a factor analysis
plot.PCA(res, axes=c(1, 2), choix="ind", habillage="none", col.ind="black", col.ind.sup="blue", 
         col.quali="magenta", label=c("ind","ind.sup", "quali"),new.plot=TRUE, title="")
plot.PCA(res, axes=c(1, 2), choix="ind", habillage="none", col.ind="black", col.ind.sup="blue", col.quali="magenta", 
         label=c("ind", "ind.sup", "quali"),new.plot=TRUE, title="")

summary(res, nb.dec = 3, nbelements=10, nbind = 10, ncp = 3, file="")
results<- HCPC(USArrests ,nb.clust=0,consol=0,min=3,max=10,cluster.CA="rows",graph=1)

results$data.clust[,ncol(results$data.clust),drop=F]
data(decathlon)
res.pca <- PCA(decathlon, quanti.sup = 11:12, quali.sup=13)
plotellipses(res.pca,keepvar=13)


# Dummifying --------------------------------------------------------------

library(dummy)
data3 <- dummy::dummy(g1g2forPCA, int = T, verbose = T) 
head(cbind(g1g2forPCA[-3], data3))
glimpse(data3)
data3 <- cbind(g1g2forPCA[-3], data3)
write.csv(data3, "check_immunotype_clustering_g1g2.csv")


# Plotting ----------------------------------------------------------------



library(rcmdr)
data3.PCA = data3[-1]
res<- PCA(data3.PCA , scale.unit=TRUE, ncp=5, graph = FALSE)
res<-PCA(data3.PCA , scale.unit=TRUE, ncp=5, graph = T)

res.hcpc <-HCPC(res,nb.clust=-1,consol=FALSE,min=3,max=10,graph=TRUE, iter.max =10) # Performs an agglomerative hierarchical clustering on results from a factor analysis

# i need to extract all cluster1, gluster 2,???3, 4,5 and plot

# which label is corrleated to cluster1 
dimdesc(res, axes=1)  # show correlation of variables with 1st axis
dimdesc(res, axes=2) 
#$`Dim.2`
#$`Dim.2`$`quanti`
#correlation      p.value
#label_GI.1VsGI.2  0.71152078 0.000000e+00
#label_GI.2VsGI.2  0.07037004 5.565965e-07
#label_GI.2VsGI.1  0.03251524 2.085051e-02
#label_GI.3VsGI.2 -0.03006364 3.264872e-02
#year_diff        -0.05090558 2.957758e-04
#label_GI.1VsGI.1 -0.82650489 0.000000e+00


plot.PCA(res, axes=c(1, 2), choix="ind", habillage="none", col.ind="black", col.ind.sup="blue", 
         col.quali="magenta", new.plot=TRUE, title="")
summary(res, nb.dec = 3, nbelements=10, nbind = 10, ncp = 3, file="")

results<- HCPC(res ,nb.clust=0,consol=0,min=3,max=10,cluster.CA="rows",graph=1)

results$data.clust[,ncol(results$data.clust),drop=F]


# interactive plots -------------------------------------------------------

library(explor)
res.hcpc <-HCPC(res,nb.clust=-1,consol=FALSE,min=3,max=10,graph=TRUE, iter.max =10) # Performs an agglomerative hierarchical clustering on results from a factor analysis
explor::explor(res)
require(FactoMineR)

# New dataset with dummified type1type2 -----------------------------------
library(dplyr)
data2$type1 <- paste(data2$gen_spe_1,data2$gen_spe_2, sep = "+") 
head(data2)
data2 <- data2[3:5]
glimpse(data2)

library(dummy)
data3 <- as.data.frame(dummy::dummy(data2, int = T, verbose = T))
data2 <- as.data.frame(data2)
has_rownames(data2) # false

names <- data2[,3]
numbers <- as.character(1:5050)
data2$type1 <- paste(data2$type1, numbers, sep="_")
head(data2); 
data2 <- data2[-3]
data4.PCA <- cbind(data2, data3)
head(data4.PCA)
write.csv(data4.PCA, "check_immunotype_clustering2.csv")

library(rcmdr)
library(FactoMineR)
res<- PCA(data4.PCA , scale.unit=TRUE, ncp=5, graph = FALSE)
res<-PCA(data4.PCA , scale.unit=TRUE, ncp=5, graph = T)
res$eig # to see eigen values end % of variance explained
dimdesc(res, axes=1)  # show correlation of variables with 1st axis



res.hcpc <-HCPC(res,nb.clust=-1,consol=FALSE,min=3,max=10,graph=TRUE, iter.max =10) # Performs an agglomerative hierarchical clustering on results from a factor analysis

plot.PCA(res, axes=c(1, 2), choix="ind", habillage="none", col.ind="black", col.ind.sup="blue", 
         col.quali="magenta", label=c("ind","ind.sup", "quali"),new.plot=TRUE, title="") # draw the pca graphs

plot.PCA(res, axes=c(1, 2), choix="ind", habillage="none", col.ind="black", col.ind.sup="blue", col.quali="magenta", 
         label=c("ind", "ind.sup", "quali"),new.plot=TRUE, title="")
class(res)

plotellipses(res)

summary(res, nb.dec = 3, nbelements=10, nbind = 10, ncp = 3, file="")



results<-HCPC(data ,nb.clust=-1,consol=0,min=3,max=10,cluster.CA="rows",graph=1)
results$data.clust[,ncol(results$data.clust),drop=F]


# demo --------------------------------------------------------------------

## FactoMineR::MCA exploration
library(explor)
data(hobbies)
mca <- MCA(hobbies[1:1000,c(1:8,21:23)], quali.sup = 9:10, 
           quanti.sup = 11, ind.sup = 1:100, graph = FALSE)
explor(mca)

## FactoMineR::PCA exploration
data(decathlon)
d <- decathlon[,1:12]
pca <- PCA(d, quanti.sup = 11:12, graph = FALSE)
explor(pca)

## End(Not run)
## Not run: 

library(ade4)

data(bordeaux)
tab <- bordeaux
row_sup <- tab[5,-4]
col_sup <- tab[-5,4]
coa <- dudi.coa(tab[-5,-4], nf = 5, scannf = FALSE)
coa$supr <- suprow(coa, row_sup)$lisup
coa$supc <- supcol(coa, col_sup)$cosup
explor(coa)

## End(Not run)
## Not run: 

library(ade4)
data(banque)
d <- banque[-(1:100),-(19:21)]
ind_sup <- banque[1:100, -(19:21)]
var_sup <- banque[-(1:100),19:21]
acm <- dudi.acm(d, scannf = FALSE, nf = 5)
acm$supv <- supcol(acm, dudi.acm(var_sup, scannf = FALSE, nf = 5)$tab)$cosup
colw <- acm$cw*ncol(d)
X <- acm.disjonctif(ind_sup)
X <- data.frame(t(t(X)/colw) - 1)
acm$supi <- suprow(acm, X)$lisup
explor(acm)

## End(Not run)
## Not run: 

library(ade4)
data(deug)
d <- deug$tab
sup_var <- d[-(1:10), 8:9]
sup_ind <- d[1:10, -(8:9)]
pca <- dudi.pca(d[-(1:10), -(8:9)], scale = TRUE, scannf = FALSE, nf = 5)
supi <- suprow(pca, sup_ind)
pca$supi <- supi$lisup
supv <- supcol(pca, dudi.pca(sup_var, scale = TRUE, scannf = FALSE)$tab)
pca$supv <- supv$cosup
explor(pca)

## End(Not run)
