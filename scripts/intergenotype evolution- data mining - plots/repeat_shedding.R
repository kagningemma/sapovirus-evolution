


# compare CT values between SaV GI.1 and GI.2 -----------------------------


library(readxl) 
library(magrittr)
library(tidyverse)

TG1=read_xlsx(path = "20180810Sapovirus_PH (community study).xlsx",sheet = 1, col_names = T)
TG1 %>% dplyr::select(per_id,date_col, genotyping,`qPCR CT(1st)`)%>% 
  group_by(genotyping) %>% 
  summarise(avg_ct = mean(`qPCR CT(1st)`), 
            min_ct= min(`qPCR CT(1st)`),
            max_ct= max(`qPCR CT(1st)`),
            count=n())


TG2=read_xlsx(path = "20180810Sapovirus_PH (community study).xlsx",sheet = 2, col_names = T)
TG2 %>% dplyr::select(date_col, genotyping,`qPCR CT`) %>% 
  group_by(genotyping) %>% 
  summarise(avg_ct = mean(`qPCR CT`), 
            min_ct= min(`qPCR CT`),
            max_ct= max(`qPCR CT`),
            count = n())           
TG2
# considering all samples, there is no significant difference between ct values 
# or SaV GI.1 and GI.2 
# I need for each person id,  i select most recent detection dates, filter, compute mean


df <- data.frame(group=c(1,2,4,2,1,4,2,3,3),
                 ts=c("2014-02-13","2014-06-01","2014-02-14","2014-02-11","2013-02-01","2014-02-02","2014-03-21","2014-12-01","2014-02-11"),
                 letter=letters[1:9])
df$ts <- as.Date(df$ts,format='%Y-%m-%d')

dfo <- data.frame(df[order(df$ts,decreasing=F),],index=seq(1:nrow(df))) #reorder
mins <- tapply(dfo$index,dfo$group,min) # to find the min  date index for each group category
dfo[dfo$index %in% mins,] #subset the data to keep only ros with those previous minimum index 

library(dplyr)

group_by(df, group) %>% summarise(min = min(ts), letter = letter[which.min(ts)]) # shortcut
#

TG1[which.min(TG1$date_col), 3] # to see the index of ct value correspondig to the minimum detection date


# the above function groups by personid, then shows the 1st date of genotyping for each personid

group_by(TG1, per_id) %>% summarise(first_det_date = min(date_col))



group_by(TG1, per_id) %>% summarise(first_det_date = min(date_col), 
                                    Genotype = genotyping[which.min(date_col)],
                                    ct=`qPCR CT(1st)`[which.min(date_col)]) # see genotype and ct 

# i need to also group by genotype now to compare the mean ct value for each genotype

TG1_ct<- group_by(TG1, per_id, genotyping) %>% 
  summarise(first_det_date = min(date_col), ct=`qPCR CT(1st)`[which.min(date_col)])
colnames(TG1_ct)[2] <- "genotype"
table(TG1_ct$genotype)

TG1_ct$genotype <- ifelse(TG1_ct$genotype == "GIV", "GIV.1", TG1_ct$genotype)
hist(TG1_ct$ct)
anovamod1<- aov(ct~ genotype, data=TG1_ct)
summary(anovamod1) # p = 0.435
# there is no sign diff of mean ct value among genotypes in TG1 dataset


library(ggplot2)
ggplot(TG1_ct, aes(genotype, ct)) + geom_boxplot() + labs(title = "comparison of Ct value grouped by genotypes in TG1", 
                                                               subtitle = " No diff of ct vals across genotypes in TG1 dataset(p=0.4)", 
                                                               caption = " ct here are CT value for the index sample of each pers_id", 
                                                          tag = "TG1")
anovamod1<- aov(ct~ as.factor(genotype), data=TG1_ct)
summary(anovamod1) # p = 0.0323


# for TG2

TG2_ct<- group_by(TG2, per_id, `genotyping(T)`) %>% 
  summarise(first_det_date = min(date_col), ct=`qPCR CT`[which.min(date_col)])

colnames(TG2_ct)[2] <- "genotype"
table(TG1_ct$genotype)

TG2_ct$genotype <- ifelse(TG2_ct$genotype == "GIV", "GIV.1", TG2_ct$genotype)
TG2_ct$genotype = as.factor(TG2_ct$genotype)

anovamod2<- aov(ct~ genotype, data=TG2_ct)
summary(anovamod2) # p = 0.0323


hist(TG2_ct$ct)
leveneTest(ct ~ genotype, data=TG2_ct, center="median") # homogeneity of variance is ok



pairwise.t.test(write, ses, p.adj = "holm") 
TukeyHSD(anovamod2)
TukeyHSD(aov(ct ~ `genotyping(T)`, data=TG2_ct))


TG2_ct %>% group_by(`genotyping(T)`) %>%summarise(mean_ct = mean(ct), meadian_ct = median(ct),total= n())

library(ggplot2)

ggplot(TG2_ct, aes(genotype, ct)) + geom_boxplot() + labs(title = "comparison of Ct value grouped by genotypes in TG2", 
                                                          subtitle = "after posthoc test, SaV GV have a sign lower ct value compared to SaV GIV",
                                                          caption = " ct here are CT value for the index sample of each pers_id", 
                                                          tag = "TG2")






# repeated infection - reinfeection
# 1)Repeated Infection patients ---------------------------------------------------------------------------------------------------------
# reading the TG1 whole data in sheet 1
# TG2 data is in sheet 2


# how to group by dates and compute the date differences ? 
# method 1
# for each first isolation date of sample, extract



df <- tibble(~id,      ~date,      
              2380,    "10/30/12",    
              2380,   "10/31/12",    
              2380,  "11/1/12",  
              2380,    "11/2/12",  
              20100,   "10/30/12",    
              20100,   "10/31/12",   
              20100,   "11/1/12",   
              20100,   "11/2/12",   
              20103,   "10/30/12",
              20103,   "10/31/12")
df <- df %>% 
  mutate(date = mdy(date)) %>% 
  group_by(id) %>% 
  mutate(date_difference = as.numeric(date - first(date)))

# how to group by dates and compute the date differences ? 
# Method 2 write a function to extract days and compute date differences
library(stringr)

day_diff <- function(day) {
  days <- difftime(day, "2012-10-30", "days")
  str_extract(days, "\\-*\\d+\\.*\\d*")
}

# Then create a new column containing the day differences

df$date_difference2 <- unlist(lapply(df$date, day_diff))

df$date_difference == df$date_difference2


# Checking immunotypes group 1 to 5 ---------------------------------------

# for TG1

TG1=read_xlsx(path = "20180810Sapovirus_PH (community study).xlsx",sheet = 1, col_names = T)
TG1<- TG1 %>% dplyr::select(per_id,date_col, genotyping,`qPCR CT(1st)`)
TG1[table(TG1$per_id)>=2, ] # to filter only re-infection patients

TG1_reap<- TG1[which(table(TG1$per_id)>=2), ] 

# for each index in TG1$per_id)>=2, what is the difference between the second and subsequent dates ?



TG1 %>% dplyr::filter(genotyping=="GI.1") # to work only with GI.1

table(TG1$per_id)

TG1_reap %>% dplyr::filter(per_id  %in% c("LP078X", "LP039X", "LP083X", "LP097X", "LP107X", "LP125X", "LP163X", "LP177X"))# to filter person with GI.1 repeated infection

library(dplyr)
library(lubridate)
TG1_reap_1 <- TG1_reap %>%
  group_by(per_id) %>% 
  arrange(date_col) %>%
  mutate(diff = date_col - lag(date_col, default = first(date_col)))

TG1_reap_1$diff_days <- as.double(TG1_reap_1$diff, "days") # to convert the time from days 
TG1_reap_1


# are there children whihc shedding at d7, 14, 21, 

TG1_reap_1[which(TG2_reap_1$diff_days==7), ]
TG1_reap_1 %>% dplyr::filter(genotyping=="GI.1" & diff_days==6) 
TG1_reap_1 %>% dplyr::filter(genotyping=="GI.1" & diff_days==7) 
TG1_reap_1 %>% dplyr::filter(genotyping=="GI.1" & diff_days==10)  # there is shedding child with GI.1
TG1_reap_1 %>% dplyr::filter(genotyping=="GI.2" & diff_days==7) 



TG1_reap_1 %>% dplyr::filter(genotyping=="GI.2" & diff_days==6) 
TG1_reap_1 %>% dplyr::filter(genotyping=="GI.2" & diff_days==7) 
TG1_reap_1 %>% dplyr::filter(genotyping=="GI.2" & diff_days==10)  # there is shedding child with GI.1
TG1_reap_1 %>% dplyr::filter(genotyping=="GI.2" & diff_days==7) 





TG1_reinfected <- TG1_reap_1[which(TG1_reap_1$diff_days>25), ] # to see patient who have been re-infected by sapo 
# "LP060X", "LP002X", "LP013X", "LP038X", "LP066X"

TG1 %>% dplyr::filter(per_id  %in% c("LP060X", "LP002X", "LP013X", "LP038X", "LP066X"))# to filter persons with repeated infection


# For TG2


TG2=read_xlsx(path = "20180810Sapovirus_PH (community study).xlsx",sheet = 2, col_names = T) 
TG2<- TG2 %>% dplyr::select(per_id,date_col, 'genotyping(T)',`qPCR CT`)
table(TG2$per_id)
TG2[table(TG2$per_id)>=2, ] # to filter only re-infection patients from the main table

TG2_reap<- TG2[which(table(TG2$per_id)>=2), ] 

# for each index in TG1$per_id)>=2, what is the difference between the second and subsequent dates ?

TG2 %>% dplyr::filter(genotyping=="GI.1") # to work only with GI.1

library(dplyr)
library(lubridate)
TG2_reap_1 <- TG2_reap %>%
  group_by(per_id) %>% 
  arrange(date_col) %>%
  mutate(diff = date_col - lag(date_col, default = first(date_col)))

TG2_reap_1$diff_days <- as.double(TG2_reap_1$diff, "days") # to convert the time from days 
TG2_reap_1
TG2_reinfected <- TG2_reap_1[which(TG2_reap_1$diff_days>25), ] # to see patient who have been re-infected by sapo (re-infection after d25)
TG2_reinfected

TG2_reinfected_perid <- unlist(TG2_reinfected$per_id, use.names = F) # to get the list of id 

TG2_reinfected_sub <- TG2 %>% dplyr::filter(per_id  %in% TG2_reinfected_perid)%>% group_by(per_id) %>% arrange(date_col) # to filter persons with repeated infection from the main table
TG2_reinfected_sub # too many missing genotypes. therefore, focus only on TG1 table. 
TG2_reinfected_sub<- TG2_reinfected_sub %>% dplyr::select (per_id, genotyping, date_col,`qPCR CT` ) %>% group_by(per_id)

write.csv(TG2_reinfected_sub, " TG2_re_infected.csv")

# ask tomomi san or mayuko sensei. --> Those samples have no data.



# ask help from mathematicians and stackoverflow




# to identify shedding patients, that same person id, but again infected at 
# d7, day 14, day 21, day 28, day 35, day 42, day

TG2_reap_1[which(TG2_reap_1$diff_days==7), ]
TG2_reap_1 %>% dplyr::filter(genotyping=="GI.1" & diff_days==6) 
TG2_reap_1 %>% dplyr::filter(genotyping=="GI.1" & diff_days==7) 
TG2_reap_1 %>% dplyr::filter(genotyping=="GI.1" & diff_days==10)  # there is shedding child with GI.1
TG2_reap_1 %>% dplyr::filter(genotyping=="GI.2" & diff_days==7) 



TG2_reap_1 %>% dplyr::filter(genotyping=="GI.2" & diff_days==6) 
TG2_reap_1 %>% dplyr::filter(genotyping=="GI.2" & diff_days==7) 
TG2_reap_1 %>% dplyr::filter(genotyping=="GI.2" & diff_days==10)  # there is shedding child with GI.1
TG2_reap_1 %>% dplyr::filter(genotyping=="GI.2" & diff_days==7) 

# problem statement, using a table containing var1= make a matrix in which we can not diagno




### counting sapovirus episodes
library(readxl)
library(tibble)
library(tidyverse)
calendar <- read.csv("C:/Users/Virology Tohoku/Documents/MyRdirectory/PHIL_Data/TG1TG2_calendardata/calendar_lpsn_pcr_enddate.csv")
labo <- TG1=read_xlsx("C:/Users/Virology Tohoku/Documents/MyRdirectory/SAPO/Sapovirus evolution over time/SAPO_POS_GERARDO_TG1_2_mdfd_Emmanuel's changes 2019-04-11.xlsx",sheet = 1, col_names = T)
head(calendar)
calendar2<- calendar %>% filter(liqbm!=-1) %>% # remove missing values from liqbm col
  mutate(diarrhea = if_else(liqbm>=3,1,0))# create diarrhea column
colnames(calendar2)
table(calendar2$diarrhea)

nav_subset<- subset(calendar2, select = c(chldcode,date,diarrhea,mentstat,seizures,vomitday,anorexia,
                             fever,btemp,bwlmove,liqbm,bbdiahr,abdodist,abdopain,cough,brthdif,
                             exclubrstfeed,mxdfeed,fmlfeed,dmxfeed,fulcriteria,Study_Number,
                             Study_Number_2,Collection_Date,Cohort_ID,Cohort_ID_2,
                             Initials,Gender,Age_upon_collection,Type1,Type2,Type3,Types))
nav_subset2<- dplyr::sample_n(nav_subset, 30000, replace = T)
write.csv(nav_subset2, "nav_subset.csv")
nav_subset3 <- dplyr::sample_n(calendar2, 70000, replace=F)

# Create counting of episodes of diarrhea ------------------------------
# merging of genotype and clinical data will be done on basis of sampling data + childcode(personid)
# for each row, 1 diarrhea episode is defined >3 liqui/semiliqu stool /day, and no consecutive D for more than 2 days. 
# if several consecutive DDDDD and > 2 000 before another D, then count 1 episode. 


# Function to count the number of diarrhea episodes -----------------------


## make a function to count ones in a given subset which represents only one child

count_ones = function(x){
  ones = list()
  
  if (x[1,4] == 1 | x[2,4] == 1) {ones = append(ones, 1)}
  
  for (i in 3:nrow(x)){
    if (x[i,4] == 0 & x[i-1,4] == 0) {next} else ##just continue
      if (x[i,4] == 0 & x[i-1,4] == 1) {next} else ##just continue
        if (x[i,4] == 1 & x[i-1,4] == 1) {next} else ##checks if it is still in cycle
          if (x[i,4] == 1 & x[i-1,4] == 0 & x[i-2,4] == 1) {next} else
            if (x[i,4] == 1 & x[i-1,4] == 0 & x[i-2,4] == 0) {ones = append(ones, 1)} 
    
    ##check if it is new cycle
  }
  return(sum(as.numeric(ones)))
}

## now get unique names of each child code
child_codes = unique(nav_subset$chldcode)
all_cycles = list() ##create a list to add all cycles (Why not create a vector?)

## iterate through names and count ones
for (i in 1:length(child_codes)){
  sub_data = nav_subset[which(nav_subset$chldcode == child_codes[i]), ]
  if (sum(sub_data$diarrhea) == 0){
    all_cycles[[i]] = 0
  } else {
    sub_data2 = sub_data[order(sub_data$date),]
    all_cycles[[i]] = count_ones(sub_data2)
  }
  all_cyles = as.numeric(all_cycles)
}

## final result
##bind data with codes to make easier to read



all_cycles_result = cbind(child_codes, as.numeric(all_cycles))
all_cycles_result <- as.data.frame(all_cycles_result)

table(all_cycles_result$V2)













library(data.table)
dt <- as.data.table(data)
setkey(dt, id, date)
dt[, diff := value - shift(value, fill = first(value)), by = id]

data$diff <- ave(data_TG1$date_col, data_TG1$per_id, FUN=function(x) c(0, diff(x)))

data$diff <- ave(data$value, data$id, FUN=function(x) c(0, diff(x)))

# to select patients with Repeated Infectiion


# To select patient who are shedding GI.1 or GI.2 -------------------------


# loop to select SaV GI.1 with shedding on  a weekly basis 
