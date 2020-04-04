

# Using pairwise distance of aa matrix ------------------------------------

data = read_xlsx("SaVaminoacid_diffbigdata_4plotting.xlsx", sheet=1, col_names = T)

data <- data %>% dplyr::rename(specie.1="Species 1", specie.2= "Species 2")
table(data$specie.1) # GII X, GII.3, GII.2, GII.1, GII.6, GII.8, GII.3, GII.5, GII.4, 
# select only GI
data$gen_spe_1 = str_extract(data$specie.1, "[G][I|1][.| ][0-9 | I-IV]")
data$gen_spe_2 = str_extract(data$specie.1, "[G][I|1][.| ][0-9 | I-IV]")
# select only GII
table(str_extract(data$specie.1, "[G][I|1][.| ][0-9 | I-IV]")) # match any GI 

# Using DNA dataset -------------------------------------------------------


# Does Immunotypes A and G realy represent distinct genetic groups --------


g1g2data <-Biostrings::readDNAStringSet("sav_wholecaps_GI_GII_76.fas")

g1g2data # i want to match label and names in the fasta file. so i need to rename the fasta files names.

names(g1g2data)

names(g1g2data)[1] 

g1g2data@ranges # to see and access the ranges of each gene as well as the names. 
g1g2data@ranges

table(stri_detect_fixed(names(g1g2data), "GI.1")) # to make sure we have 39 G1.1
library(stringi)
library(dplyr)

# if G1.1 is true, then replace by immunotype A else immunotype B 

names(g1g2data) <- ifelse(stri_detect_fixed(names(g1g2data), "GI.1")==TRUE, names(g1g2data)=="A", "B")

names(g1g2data) <- factor(names(g1g2data), levels = c("FALSE","B"), labels = c("A", "B"))

# export the fasta file and run p distances in mega

Biostrings::writeXStringSet(g1g2data,format="fasta")

out23a <- tempfile()
Biostrings::writeXStringSet(g1g2data, out23a)# export the fasta file
out23a # go to the location and copy paste the file.


# after checking the p distances, i found that SaV GI.1 and GI.2 form a separate immunotype

# My tableau charts are found here 
https://public.tableau.com/profile/emmanuel3421#!/vizhome/SapovirusEvolutionresearch/Dashboard1?publish=yes


# Creating aligmenent with immunotypes A, B, C and D -----------------------------------------------
setwd("~/MyRdirectory/SAPO/Sapovirus evolution over time")

allfasta=Biostrings::readDNAStringSet("allSaVdata_for pairwise_ok.fas")
mtcars[grep("Merc", rownames(mtcars)), ] # to select all car rownames starting with "merc"

DandC<- allfasta[grep("GI.7|GI.3|GI.6|GI.5|GI.4", names(allfasta)), ] # to select only fasta sequences containing either of those genotypes
grepl("GI.7|GI.3", names(DandC))
names(DandC) %in% c("GI.3", grep("GI.7", names(DandC), value=T)) #if name contains GI.3 or GI.7, say true
table(names(DandC) %in% c("GI.3", grep("GI.7", names(DandC), value=T)))

names(DandC)<- ifelse(names(DandC) %in% c("GI.3", grep("GI.7", names(DandC), value=T))==T, "D", "C")




out23a <- tempfile()
Biostrings::writeXStringSet(DandC, out23a)# export the fasta file
out23a # go to the location and copy paste the temp file into your working directory. rename the extension of file into .fas

allfasta3<- allfasta[grep("GI.2", names(allfasta)), ] 
allfasta4<- allfasta[grep("GI.1", names(allfasta)), ] 
g1.1 <- tempfile()
g1.2<- tempfile()
writeXStringSet(allfasta4,  g1.1) # to export GI.1 sequences
writeXStringSet(allfasta3,  g1.2) # to export GI.2 sequences
g1.1 # to see fasta file containing only G1.1 sequences_very fast method, rather than doing manually
g1.2



# Creating aligmenent file with tentative immunotypes group1,2,3,4,5-------------------------------


setwd("~/MyRdirectory/SAPO/Sapovirus evolution over time")

allfasta=Biostrings::readDNAStringSet("allSaVdata_for pairwise_ok.fas")

mtcars[grep("Merc", rownames(mtcars)), ] # to select all car names starting with "merc"

group1to5 <- allfasta[grep("GI.1|GI.2|GI.7|GI.3|GI.6|GI.5|GI.4|GII.1|GII.2|GII.7|GII.3|GII.6|GII.5|GII.4|GII.8", names(allfasta)), ] # to select only fasta sequences containing either of those genotypes

grepl("GI.1", names(group1to5))

table(names(group1to5) %in% c("GI.3")

names(group1to5) %in% c("GI.3", grep("GI.7", names(group1to5), value=T)) #if name contains GI.3 or GI.7, say true

names(group1to5) <- ifelse(grepl("GI.1", names(group1to5))==T, "group1", names(group1to5)) # naming group 1
names(group1to5) <- ifelse(grepl("GI.2|GI.6|GI.5", names(group1to5))==T, "group2", names(group1to5)) # naming group 2

names(group1to5) <- ifelse(grepl("GI.7|GI.3", names(group1to5))==T, "group3", names(group1to5)) # naming group 3

names(group1to5) <- ifelse(grepl("GI.4", names(group1to5))==T, "group4", names(group1to5)) # naming group 4

names(group1to5) <- ifelse(grepl("GII.1|GII.2|GII.3|GII.4|GII.5|GII.6", names(group1to5))==T, "group5", names(group1to5)) # naming group 5
names(group1to5) <- ifelse(grepl("GII.7|GII.8", names(group1to5))==T, "group6", names(group1to5)) # naming group 6
table(names(group1to5))

out23a <- tempfile()
Biostrings::writeXStringSet(group1to5, out23a)# export the fasta file
out23a # go to the location and copy paste the temp file into your working directory. rename the extension of file into .fas
# useMEGA to compute pairwise dist and use tableau to visualize misclassified cases

BiocManager::install("Biostrings")










# clinical validation of immunotypes classification -----------------------

#open sheeet 4 in this file file https://www.dropbox.com/s/rrlkxs4y9qdgnrg/repeated%20infection%20philippines_301018.xlsx?dl=0


# Help from stacks overflow -----------------------------------------------
setwd("~/MyRdirectory/SAPO/Sapovirus evolution over time/immunotypes_g1g12")
library(readxl)
library(dplyr)
library(plyr)
data <- read_excel("mathias.xlsx", sheet = 2, 
                   col_names = T, trim_ws = T, range = "B1:C54")
data$sv_type1 <- plyr::revalue(data$sv_type1, c("GIV"="GIV.1"))
presentinboth<-intersect(data$sv_type1, data$sv_type2)
data %>% select(presentinboth)
table(data$sv_type1, data$sv_type2)



table(data$sv_type1, data$sv_type2)
#sapply(data$sv_type1, class) # see class of all items of a vector
#lapply(data$sv_type1, class) # column
table(data$sv_type1)
table(data$sv_type2)


# rename levels
data$sv_type1 <- plyr::revalue (data$sv_type1, 
                                 c("GI.1"="A" , "GI.2"="B", "GI.6"="C", "GII.1"="D",
                                   "GII.3"="E", "GII.4" = "F", 
                                   "GII.4 recombinant"="G", "GII.6"="H", "GIV.1"="I",
                                   "GV.1"="J"))
data$sv_type2 <- plyr::revalue (data$sv_type2, 
                                c("GI.1"="A" , "GI.2"="B", "GI.6"="C", 
                                  "GII.1"="D","GII.2"= "K",
                                  "GII.3"="E", "GII.4" = "F", "GIV"="I",                                   "GII.4 recombinant"="G", "GII.6"="H", "GIV.1"="I",
                                  "GV.1"="J"))
data <- data %>% rename(type1 = sv_type1, type2 = sv_type2)
table(data$sv_type1, data$sv_type2)
data2 <- data[-c(1:2),]
data2<- data2 %>% select(sv_type1, sv_type2)
data2 <- data2[1:20, ]
table(data2$type1, data2$type2)

dput(as.data.frame(data2))
# Hi to everyone. 
#I don't knnow the mathematical or reproducible method to solve the following problem. 
#I've got the following dataframe from a study we are conducting. 
# The data contain 2 columns sv_type1 and sv_type2  and items A, B, C etc. in each rows 
# The following matrix was generated from that dataframe. 
# In that matrix M1, values in the diagnonal are either 0 or 1.
# I would like to generate a new dataframe (or all possible dataframes) in which surrogate items (e.g: A+B, A+C) and corresponding matrix in which diagnonal values are either 1 or 0.
# Here is an example of expected matrix :

dataset: 
matrix:
diag(1,5)





#v1.2 - New re-infection data analysis - v1.2 -------------------------------------



setwd("~/MyRdirectory/SAPO/Sapovirus evolution over time/Immunotypes")
library(readxl)
library(dplyr)
library(plyr)
data_peru <- read_excel("reinfection_clinical_data.xlsx", sheet = 2, 
                   col_names = T, trim_ws = T)
data_phil <- read_excel("reinfection_clinical_data.xlsx", sheet = 1, 
                        col_names = T, trim_ws = T)
head(data_peru);head(data_phil)

table(data_peru$genotype1, data_peru$genotype2)
#sapply(data$sv_type1, class) # see class of all items of a vector
#lapply(data$sv_type1, class) # column
table(data_peru$genotype1)
table(data_peru$genotype1)
# rename levels
# exclude GII.N
data_peru$new_class1 <- plyr::revalue (data_peru$genotype1, 
                                c("GI.1"="group1" , "GI.2"="group2", "GI.6"="group2", "GI.7"="group3",
                                  "GII.1"="group5", "GII.2"="group5", "GII.3"="group5","GII.4"="group5", 
                                  "GII.5" = "group5", "GII.6"="group5", "GII.7"="group6", 
                                  "GII.8"="group6"))
data_peru$new_class2 <- plyr::revalue (data_peru$genotype2, 
                                  c("GI.1"="group1" , "GI.2"="group2", "GI.6"="group2", "GI.7"="group3",
                                    "GII.1"="group5", "GII.2"="group5", "GII.3"="group5", "GII.4"="group5", 
                                    "GII.5" = "group5", "GII.6"="group5", "GII.7"="group6","GII.8"="group6")) # The following `from` values were not present in `x`: GII.7

head(data_phil)
table(data_phil$sv_type1)
table(data_phil$sv_type2)


data_phil$new_class1 <- plyr::revalue (data_phil$sv_type1, 
                                       c("GI.1"="group1" , "GI.2"="group2", "GI.6"="group2", "GI.7"="group3", "GI.4"="group4",
                                         "GII.1"="group5", "GII.2"="group5", "GII.3"="group5","GII.4"="group5", "GII.4 recombinant"="group4",
                                         "GII.5" = "group5", "GII.6"="group5", "GII.7"="group6", "GIV"="GIV.1",
                                         "GII.8"="group6"))
data_phil$new_class2 <- plyr::revalue (data_phil$sv_type2, 
                                       c("GI.1"="group1" , "GI.2"="group2", "GI.6"="group2", "GI.7"="group3","GI.4"="group4",
                                         "GII.1"="group5", "GII.2"="group5", "GII.3"="group5", "GII.4"="group5", "GII.4 recombinant"="group4",
                                         "GII.5" = "group5", "GII.6"="group5", "GII.7"="group6","GII.8"="group6", "GIV"="GIV.1")) # The following `from` values were not present in `x`: GII.7


matrix_peru <- table(data_peru$new_class1, data_peru$new_class2)
matrix_phil <- table(data_phil$new_class1, data_phil$new_class2)
# v1 : i used groupA, B, C, D, E, F, G, H
# in version 2, i did what is above
# vversion 3, i what is below 
write.csv(matrix_peru, "reinfection_matrix_peru.csv")
write.csv(matrix_phil, "reinfection_matrix_phil.csv")



# v2 - New grouping ------------------------------------------------------

data_peru <- read_excel("reinfection_clinical_data.xlsx", sheet = 2, col_names = T, trim_ws = T)
data_phil <- read_excel("reinfection_clinical_data.xlsx", sheet = 1, col_names = T, trim_ws = T)

data_peru$new_class1 <- plyr::revalue (data_peru$genotype1, 
                                       c("GI.1"="group1" , "GI.2"="group2", "GI.6"="group7", "GI.7"="group3",
                                         "GII.1"="group5", "GII.2"="group5", "GII.3"="group8","GII.4"="group5", 
                                         "GII.5" = "group5", "GII.6"="group5", "GII.7"="group6", "GII.8"="group6"))
data_peru$new_class2 <- plyr::revalue (data_peru$genotype2, 
                                       c("GI.1"="group1" , "GI.2"="group2", "GI.6"="group7", "GI.7"="group3",
                                         "GII.1"="group5", "GII.2"="group5", "GII.3"="group8", "GII.4"="group5", 
                                         "GII.5" = "group5", "GII.6"="group5", "GII.7"="group6","GII.8"="group6")) # The following `from` values were not present in `x`: GII.7
data_phil$new_class1 <- plyr::revalue (data_phil$sv_type1, 
                                       c("GI.1"="group1" , "GI.2"="group2", "GI.6"="group7", "GI.7"="group3", "GI.4"="group4",
                                         "GII.1"="group5", "GII.2"="group5", "GII.3"="group8","GII.4"="group5", "GII.4 recombinant"="group5",
                                         "GII.5" = "group5", "GII.6"="group5", "GII.7"="group6", "GIV"="GIV.1","GII.8"="group6"))
data_phil$new_class2 <- plyr::revalue (data_phil$sv_type2, 
                                       c("GI.1"="group1" , "GI.2"="group2", "GI.6"="group7", "GI.7"="group3","GI.4"="group4",
                                         "GII.1"="group5", "GII.2"="group5", "GII.3"="group8", "GII.4"="group5", "GII.4 recombinant"="group5",
                                         "GII.5" = "group5", "GII.6"="group5", "GII.7"="group6","GII.8"="group6", "GIV"="GIV.1")) # The following `from` values were not present in `x`: GII.7

data_phil %>% filter(new_class1 == "group5") %>% filter(new_class2 == "group5") # tosee same genotype re-infection samples
data_phil %>% filter(new_class1 == "group2") %>% filter(new_class2 == "group2") # tosee same genotype re-infection samples
data_phil %>% filter(new_class1 == "group1") %>% filter(new_class2 == "group1") # tosee same genotype re-infection samples
data_peru %>% filter(new_class1 == "group5") %>% filter(new_class2 == "group5") # tosee same genotype re-infection samples



matrix_peru <- table(data_peru$new_class1, data_peru$new_class2)
matrix_phil <- table(data_phil$new_class1, data_phil$new_class2)
matrix_phil
matrix_peru
write.csv(matrix_peru, "2-reinfection_matrix_peru.csv")
write.csv(matrix_phil, "2-reinfection_matrix_phil.csv")






# New grouping- Version 3------------------------------------------------------------

# v3 - New grouping ------------------------------------------------------

data_peru <- read_excel("reinfection_clinical_data.xlsx", sheet = 2, col_names = T, trim_ws = T)
data_phil <- read_excel("reinfection_clinical_data.xlsx", sheet = 1, col_names = T, trim_ws = T)

data_peru$new_class1 <- plyr::revalue (data_peru$genotype1, 
                                       c("GI.1"="group1" , "GI.2"="group2", "GI.6"="group7", "GI.7"="group3",
                                         "GII.1"="group5", "GII.2"="group5", "GII.3"="group8","GII.4"="group5", 
                                         "GII.5" = "group9", "GII.6"="group5", "GII.7"="group6", "GII.8"="group6"))
data_peru$new_class2 <- plyr::revalue (data_peru$genotype2, 
                                       c("GI.1"="group1" , "GI.2"="group2", "GI.6"="group7", "GI.7"="group3",
                                         "GII.1"="group5", "GII.2"="group5", "GII.3"="group8", "GII.4"="group5", 
                                         "GII.5" = "group9", "GII.6"="group5", "GII.7"="group6","GII.8"="group6")) # The following `from` values were not present in `x`: GII.7
data_phil$new_class1 <- plyr::revalue (data_phil$sv_type1, 
                                       c("GI.1"="group1" , "GI.2"="group2", "GI.6"="group7", "GI.7"="group3", "GI.4"="group4",
                                         "GII.1"="group5", "GII.2"="group5", "GII.3"="group8","GII.4"="group5", "GII.4 recombinant"="group5",
                                         "GII.5" = "group9", "GII.6"="group5", "GII.7"="group6", "GIV"="GIV.1","GII.8"="group6"))
data_phil$new_class2 <- plyr::revalue (data_phil$sv_type2, 
                                       c("GI.1"="group1" , "GI.2"="group2", "GI.6"="group7", "GI.7"="group3","GI.4"="group4",
                                         "GII.1"="group5", "GII.2"="group5", "GII.3"="group8", "GII.4"="group5", "GII.4 recombinant"="group5",
                                         "GII.5" = "group9", "GII.6"="group5", "GII.7"="group6","GII.8"="group6", "GIV"="GIV.1")) # The following `from` values were not present in `x`: GII.7


data_phil %>% filter(new_class1 == "group1") %>% filter(new_class2 == "group1") # to see same genotype re-infection samples
data_phil %>% filter(new_class1 == "group2") %>% filter(new_class2 == "group2") # tosee same genotype re-infection samples
data_phil %>% filter(new_class1 == "group5") %>% filter(new_class2 == "group5") # tosee same genotype re-infection samples

data_peru %>% filter(new_class1 == "group1") %>% filter(new_class2 == "group1") # to see same genotype re-infection samples
data_peru %>% filter(new_class1 == "group5") %>% filter(new_class2 == "group5") # to see same genotype re-infection samples


matrix_peru_v3 <- table(data_peru$new_class1, data_peru$new_class2)
matrix_phil_v3<- table(data_phil$new_class1, data_phil$new_class2)
matrix_peru_v3
matrix_phil_v3
write.csv(matrix_peru_v3, "3-reinfection_matrix_peru.csv")
write.csv(matrix_phil_v3, "3-reinfection_matrix_phil.csv")



# v4 - New grouping ------------------------------------------------------

data_peru <- read_excel("reinfection_clinical_data.xlsx", sheet = 2, col_names = T, trim_ws = T)
data_phil <- read_excel("reinfection_clinical_data.xlsx", sheet = 1, col_names = T, trim_ws = T)

data_peru$new_class1 <- plyr::revalue (data_peru$genotype1, 
                                       c("GI.1"="group1" , "GI.2"="group2", "GI.6"="group7", "GI.7"="group3",
                                         "GII.1"="group5", "GII.2"="group10", "GII.3"="group8","GII.4"="group5", 
                                         "GII.5" = "group9", "GII.6"="group5", "GII.7"="group6", "GII.8"="group6"))
data_peru$new_class2 <- plyr::revalue (data_peru$genotype2, 
                                       c("GI.1"="group1" , "GI.2"="group2", "GI.6"="group7", "GI.7"="group3",
                                         "GII.1"="group5", "GII.2"="group10", "GII.3"="group8", "GII.4"="group5", 
                                         "GII.5" = "group9", "GII.6"="group5", "GII.7"="group6","GII.8"="group6")) # The following `from` values were not present in `x`: GII.7
data_phil$new_class1 <- plyr::revalue (data_phil$sv_type1, 
                                       c("GI.1"="group1" , "GI.2"="group2", "GI.6"="group7", "GI.7"="group3", "GI.4"="group4",
                                         "GII.1"="group5", "GII.2"="group10", "GII.3"="group8","GII.4"="group5", "GII.4 recombinant"="group5",
                                         "GII.5" = "group9", "GII.6"="group5", "GII.7"="group6", "GIV"="GIV.1","GII.8"="group6"))
data_phil$new_class2 <- plyr::revalue (data_phil$sv_type2, 
                                       c("GI.1"="group1" , "GI.2"="group2", "GI.6"="group7", "GI.7"="group3","GI.4"="group4",
                                         "GII.1"="group5", "GII.2"="group10", "GII.3"="group8", "GII.4"="group5", "GII.4 recombinant"="group5",
                                         "GII.5" = "group9", "GII.6"="group5", "GII.7"="group6","GII.8"="group6", "GIV"="GIV.1")) # The following `from` values were not present in `x`: GII.7


data_phil %>% filter(new_class1 == "group1") %>% filter(new_class2 == "group1") # to see same genotype re-infection samples
data_phil %>% filter(new_class1 == "group5") # to see same genotype re-infection samples


matrix_peru_v4 <- table(data_peru$new_class1, data_peru$new_class2)
matrix_phil_v4<- table(data_phil$new_class1, data_phil$new_class2)
matrix_peru_v4
matrix_phil_v4
write.csv(matrix_peru_v3, "4-reinfection_matrix_peru.csv")
write.csv(matrix_phil_v3, "4-reinfection_matrix_phil.csv")


# v5-Groupings ------------------------------------------------------------


data_peru <- read_excel("reinfection_clinical_data.xlsx", sheet = 2, col_names = T, trim_ws = T)
data_phil <- read_excel("reinfection_clinical_data.xlsx", sheet = 1, col_names = T, trim_ws = T)

# let us omit those what were incorrecly labelled as re-infection cases 
# remove all GI.1-GI.1 re-infections in phillipines, because after checking with mayuko sensei, we think they are all same 
data_phil <- data_phil %>% dplyr::filter(is.na(Note)) %>% filter(per_id!="LP100X") # LP100X is GII.6 and GII.6 -- tomomisan found they both mapped to a different reference
data_phil %>% filter(sv_type1== "GI.1") %>% filter(sv_type2== "GI.1") # there is none
data_peru <- data_peru %>% filter(per_id!="px029")
data_peru %>% filter(genotype1== "GI.1") %>% filter(genotype2== "GI.1")

# let us re-label the genotypes into immunotypes
data_peru$new_class1 <- plyr::revalue (data_peru$genotype1, 
                                       c("GI.1"="group1" , "GI.2"="group2", "GI.6"="group7", "GI.7"="group3",
                                         "GII.1"="group5", "GII.2"="group10", "GII.3"="group8","GII.4"="group5", 
                                         "GII.5" = "group9", "GII.6"="group5", "GII.7"="group6", "GII.8"="group6"))
data_peru$new_class2 <- plyr::revalue (data_peru$genotype2, 
                                       c("GI.1"="group1" , "GI.2"="group2", "GI.6"="group7", "GI.7"="group3",
                                         "GII.1"="group5", "GII.2"="group10", "GII.3"="group8", "GII.4"="group5", 
                                         "GII.5" = "group9", "GII.6"="group5", "GII.7"="group6","GII.8"="group6")) # The following `from` values were not present in `x`: GII.7
data_phil$new_class1 <- plyr::revalue (data_phil$sv_type1, 
                                       c("GI.1"="group1" , "GI.2"="group2", "GI.6"="group7", "GI.7"="group3", "GI.4"="group4",
                                         "GII.1"="group5", "GII.2"="group10", "GII.3"="group8","GII.4"="group5", "GII.4 recombinant"="group5",
                                         "GII.5" = "group9", "GII.6"="group5", "GII.7"="group6", "GIV"="GIV.1","GII.8"="group6"))
data_phil$new_class2 <- plyr::revalue (data_phil$sv_type2, 
                                       c("GI.1"="group1" , "GI.2"="group2", "GI.6"="group7", "GI.7"="group3","GI.4"="group4",
                                         "GII.1"="group5", "GII.2"="group10", "GII.3"="group8", "GII.4"="group5", "GII.4 recombinant"="group5",
                                         "GII.5" = "group9", "GII.6"="group5", "GII.7"="group6","GII.8"="group6", "GIV"="GIV.1")) # The following `from` values were not present in `x`: GII.7

target = c("NA","GIV.1","GV.1","GII.N") # i need to remove all cases with either NA, GIV.1, GV.1 and GII.N
data_phil <- data_phil %>% dplyr::filter(!(new_class1 %in% target)) %>% dplyr::filter(!(new_class2 %in% target)) %>% 
  select(new_class1, new_class2)
data_peru <- data_peru %>% dplyr::filter(!(new_class1 %in% target)) %>% dplyr::filter(!(new_class2 %in% target))%>% 
  select(new_class1, new_class2)

#write.csv(data_peru, "5-data_peru_v5.csv")
#write.csv(data_phil, "5-data_phil_v5.csv")


table(data_peru$new_class1) # omit group 6
table(data_peru$new_class2) # i need to omit group 8 and group3


which(data_peru$new_class2=="group3") #  4, 6
which(data_peru$new_class2=="group8") # 23
which(data_peru$new_class1=="group6") # 4 22
table(data_peru[-c(4,22), ]$new_class1) # ok
table(data_peru[-c(4,6,22), ]$new_class2) # ok

data_peru1 <- data_peru[-c(4,6,22,23), ]
colnames(data_peru1)

table(data_peru1$new_class1) 
table(data_peru1$new_class2)


table(data_phil$new_class1) 
table(data_phil$new_class2) # omit all lines containing group 10, 4

which(data_phil$new_class2=="group10") # 20
which(data_phil$new_class2=="group4") # 42,43,46
table(data_phil[-c(20,42,43,46), ]$new_class2) 
data_phil1 <- data_phil[-c(20,42,43,46), ]
colnames(data_phil1)

# TO get the contigency table for peru dataset and phil dataset

matrix_peru_v5 <- table(data_peru1$new_class1, data_peru1$new_class2); matrix_peru_v5
matrix_phil_v5<- table(data_phil1$new_class1, data_phil1$new_class2); matrix_phil_v5



write.csv(matrix_peru_v5, "5-reinfection_matrix_peru.csv")
write.csv(matrix_phil_v5, "5-reinfection_matrix_phil.csv")



## I need to MERGE Phil and peru data before making contigency table

# TO GET contingency table or re-infection matrix for the merged dataset
dim(data_phil1) #42  2
dim(data_peru1) # 24, 2 
# i expect the merged table to have 71 rows and 2 cols

data_peru_phil <- rbind(data_peru1, data_phil1)
# A tibble: 71 x 2


matrix_mix <- table(data_peru_phil$new_class1, data_peru_phil$new_class2); matrix_mix
write.csv(matrix_mix, "matrix_peru_phil_merged.csv")












# let us check the data and cut-off

# second tentative immunotypes group1,2,3,4,5-------------------------------


setwd("~/MyRdirectory/SAPO/Sapovirus evolution over time")
#BiocManager::install("Biostrings")
library(Biostrings)
getwd()
allfasta=Biostrings::readDNAStringSet("allSaV_capsid_data_for pairwise_ok.fas")

mtcars[grep("Merc", rownames(mtcars)), ] # to select all car names starting with "merc"

group1to10 <- allfasta[grep("GI.1|GI.2|GI.7|GI.3|GI.6|GI.5|GI.4|GII.1|GII.2|GII.7|GII.3|GII.6|GII.5|GII.4|GII.8", names(allfasta)), ] # to select only fasta sequences containing either of those genotypes
names(group1to10)

table(names(group1to10) %in% c("GI.3"),names(group1to10) %in% c("GI.3", grep("GI.7", names(group1to10), value=T)) #if name contains GI.3 or GI.7, say true
      
names(group1to10) <- ifelse(grepl("GI.1", names(group1to10))==T, "group1", names(group1to10)) # naming group 1
names(group1to10) <- ifelse(grepl("GI.2", names(group1to10))==T, "group2", names(group1to10)) # naming group 2
names(group1to10) <- ifelse(grepl("GI.7", names(group1to10))==T, "group3", names(group1to10)) # naming group 3
names(group1to10) <- ifelse(grepl("GI.4|GI.3|GI.5", names(group1to10))==T, "group4", names(group1to10)) # naming group 4
names(group1to10) <- ifelse(grepl("GII.1|GII.4|GII.6", names(group1to10))==T, "group5", names(group1to10)) # naming group 5
names(group1to10) <- ifelse(grepl("GII.7|GII.8", names(group1to10))==T, "group6", names(group1to10)) # naming group 6
names(group1to10) <- ifelse(grepl("GI.6", names(group1to10))==T, "group7", names(group1to10)) 
names(group1to10) <- ifelse(grepl("GII.3", names(group1to10))==T, "group8", names(group1to10)) 
names(group1to10) <- ifelse(grepl("GII.5", names(group1to10))==T, "group9", names(group1to10)) 
names(group1to10) <- ifelse(grepl("GII.2", names(group1to10))==T, "group10", names(group1to10)) 
table(names(group1to10))

      
out23a <- tempfile()
Biostrings::writeXStringSet(group1to10, out23a)# export the fasta file
out23a # go to the location and copy paste the temp file into your working directory. rename the extension of file into .fas
      # the file was saved at C:\Users\Virology Tohoku\Documents\MyRdirectory\SAPO\Sapovirus evolution over time\Immunotypes\group1to10

#based on aa p distance computed using MEGA, i found aa distance between groups is >0.16 and distance within groups is <0.14
# the aa dist file is group1to10.csv

# Confirming the p-distance cut-off ---------------------------------------
library(tibble)
library(dplyr)
library(stringr)
setwd("~/MyRdirectory/SAPO/Sapovirus evolution over time/Immunotypes/group1to10")
group1to10_aa <- as_tibble(read.csv("group1to10.csv"))
psych::headtail(group1to10_aa)
#tail(group1to10_aa) # the caption form MEGA is at the end
#group1to10_aa[14701:14709, ]
#group1to10_aa <- group1to10_aa[1:14706, ]
# df$C <- ifelse(grepl("D", df$A), "yes", "no")
#grepl("Gp 9",group1to10_aa$Species.1)
#grepl("Gp 9",group1to10_aa$Species.2)
# group1to10_aa$class <- base::ifelse(grepl("Gp 9",group1to10_aa$Species.1=="T"& grepl("Gp 9",group1to10_aa$Species.2))=="T", "intra", "inter") 
# if specie.1 and specie.2 contains "gp 9", then consider it as intergroup


head(stringr::str_sub(group1to10_aa$Species.1,1,6)) # to extract the first 6 strings
group1to10_aa$species.1<- stringr::str_sub(group1to10_aa$Species.1,1,6)
group1to10_aa$species.2<- stringr::str_sub(group1to10_aa$Species.2,1,6)
group1to10_aa <- group1to10_aa %>% 
  mutate(class= ifelse(group1to10_aa$species.1==group1to10_aa$species.2, "intra", "inter")) %>%
  mutate(class=as.factor(class)) %>% 
  mutate(name= paste(species.1, species.2, sep="-")) %>% 
  mutate(Species.1=as.factor(Species.1)) %>% 
  mutate(Species.2=as.factor(Species.2))

#group1to10_aa$class <- ifelse(group1to10_aa$Species.1==group1to10_aa$Species.2, "intra", "inter")

str(group1to10_aa)

summary(group1to10_aa)
psych :: describe(group1to10_aa)
group1to10_aa %>% group_by(class) %>% 
  summarise(avg = mean(Dist), 
            med = median(Dist),
            min = min(Dist), 
            max = max(Dist))
psych::describe.by(group1to10_aa, group=group1to10_aa$class, mat=TRUE)

# psych::multi.hist(group1to10_aa)
library(psych)

d <- outlier(group1to10_aa)
library(ggplot2)
ggplot(group1to10_aa, aes(x= class, y=Dist)) +geom_boxplot() +  
  ggtitle("Amino acid differences within and between the proposed immunotypes") + 
  xlab("category") + ylab("Distance")+
  scale_shape_manual(values = c(3,16))+
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 20),
        plot.title = element_text(color="black", size=15, vjust = 0.5 ),
        axis.title.x = element_text(color="black", size=20),
        axis.title.y = element_text(color="black", size=20),
        legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"),
        legend.justification=c(1,0),
        legend.position=c(1,0), 
        legend.text=element_text(size=20))



boxplot(Dist ~ class, data = group1to10_aa, outpch = NA) 
stripchart(Dist ~ class, data = group1to10_aa, vertical = TRUE, 
           method = "jitter", 
           pch = 21, col = "maroon", bg = "bisque", add = TRUE) 

library(ggpubr)
ggboxplot(group1to10_aa,  x = "class", y = "Dist", color = "class", palette = "jco")

ggdotplot(group1to10_aa,  x = "class", y = "Dist", 
          color = "class", binwidth=.0001), 
add_summary(bxp, fun="mean_se", error.plot = "pointrange", color="black", fill = "white")


ggline(group1to10_aa, x = "class", y = "Dist", 
       add = c("mean_se", "jitter"), title = "lineplot of aa distances")

ggbarplot(group1to10_aa, x = "class", y = "Dist", 
          add = c("mean_se", "jitter"), title = "barplot and mean+-se of aa distances")

# let us add labels
ggbarplot(group1to10_aa, x = "class", y = "Dist", 
          add = c("mean_se", "point"), title = "barplot  of aa distances") +
  ggrepel::geom_text_repel(aes(label = name))



ggpubr::ggerrorplot(group1to10_aa, x = "class", y = "Dist")
ggpubr::ggerrorplot(group1to10_aa, x = "class", y = "Dist", color = "black",
            add = "jitter", add.params = list(color = "darkgray"))




# this is my tabble of immunotypes as i made after merging peru and phil data

.Table <- matrix(c(2,2,4,7,3,2,0,0,0,1,2,0,0,0,1,0,1,2,3,2,0,6,0,5,0,3,2,3,1,0,3,2,0,2,0,0,0,2,0,3,1,0,0,1,0,0,0,0,0), 7, 7, byrow=TRUE)
chisq.test(.Table, correct=FALSE) # assumptions not met
fisher0<- fisher.test(.Table, alternative = "two.sided",simulate.p.value=TRUE,B=1e7)# us this if the chisquare assumptions are not met
fisher0 # p=0.0773

#If we remove the two same GI.1 re-infection

.Table2 <- matrix(c(0,2,4,7,3,2,0,0,0,1,2,0,0,0,1,0,1,2,3,2,0,6,0,5,0,3,2,3,1,0,3,2,0,2,0,0,0,2,0,3,1,0,0,1,0,0,0,0,0), 7, 7, byrow=TRUE)
fisher1<- fisher.test(.Table, alternative = "two.sided",simulate.p.value=TRUE,B=1e7)# us this if the chisquare assumptions are not met
fisher1 # p-value = 0.07755
