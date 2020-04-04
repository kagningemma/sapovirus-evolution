data <- 
  readXL("C:/Users/Virology Tohoku/Dropbox/phil_data/PRIMAL_sample_conc_ms_saka_kte.xlsx",
         rownames=FALSE, header=TRUE, na="", sheet="ivar_whole_genome", 
         stringsAsFactors=TRUE)
data <- within(data, {
  genotype <- factor(genotype, labels=c('GI.1','GI.2'))
})
data_peru <- subset(data, subset=country==2)















#plotting modelling of number of snps
# we assume that there is a linear relation
# Mathematical equation for number of snps as a function of days

# make hist of : y=snps
# x=breadth, country, day.difference
setwd("~/MyRdirectory/SAPO/Sapovirus evolution over time/iVAR")
data <- read.csv("results_ivar.csv")
head(data)
data <- within(data,{
  Country <- as.factor(Country)
})
# check linear relationship
scatter.smooth(x=data$Day.difference, 
               y=data$number.of.SNPs, 
               main="snps ~ breadth") 
library(ggplot2)
ggplot2::ggplot(data = data, aes(Day.difference,number.of.SNPs)) +
  geom_jitter()

library("SmartEDA")
ExpNumStat(data,by="A",gp=NULL,Qnt=seq(0,1,0.1),MesofShape=2,Outlier=TRUE,round=2)
ExpNumViz(data,nlim=10,Page=c(2,2),sample=8)

scatterplot(number.of.SNPs~Day.difference | Country, regLine=FALSE, 
            smooth=FALSE, boxplots=FALSE, by.groups=TRUE, data=data)

scatterplot(number.of.SNPs~Day.difference | Country, regLine=T, 
            smooth=FALSE, boxplots=FALSE, by.groups=TRUE, data=data)


scatterplot(number.of.SNPs~Day.difference | Country, regLine=F, 
            smooth=FALSE, boxplots=TRUE, by.groups=TRUE, data=data)
### this add the labels to the points using their rownames
### font = 2 is bold
text(number.of.SNPs~Day.difference, labels=data$Sample,data=data, cex=0.5, font=1)

data <- within(data,{
  Country <- as.factor(Country)
  MiseQ.run.batch. <- as.factor(MiseQ.run.batch.)
})


require("ggplot2")
.df <- data.frame(x = data$Day.difference, y = 
                    data$number.of.SNPs, z = data$Country)
.plot <- ggplot(data = .df, aes(x = x, y = y, colour = z, 
                                shape = z)) + 
  geom_point() + 
  scale_y_continuous(expand = c(0.01, 0)) + 
  xlab("Day difference") + 
  ylab("number of SNPs") + 
  labs(colour = "Country", shape = "Country") + 
  theme_bw(base_size = 14, base_family = "sans") + 
  theme(legend.position = "right")
print(.plot)

data <- within(data,{
  Country <- as.integer(Country)
  MiseQ.run.batch. <- as.integer(MiseQ.run.batch.)
})


# based on the scatterplot, we can seen that : 
# number of snps can be separated by the country.

# Using BoxPlot To Check For Outliers
# Generally, an outlier is any datapoint that lies outside the 1.5 * inter quartile range (IQR).

# IQR is calculated as the distance between the 25th percentile and 75th percentile values for that variable.

par(mfrow=c(1, 2))  # divide graph area in 2 columns
boxplot(data$number.of.SNPs, main="number of snps", 
        sub=paste("Outlier rows: ", boxplot.stats(data$number.of.SNPs)$out))  # box plot for 'speed'
boxplot(data$Day.difference, main="day difference", 
        sub=paste("Outlier rows: ", boxplot.stats(data$Day.difference)$out))  # box plot for 'distance'

# based on this plot, i can see some outlying values.
#there value will affect my  Rsquare


# Using Density Plot To Check If Response Variable Is Close To Normal
library(e1071)  # for skewness function
par(mfrow=c(1, 2))  # divide graph area in 2 columns

plot(density(na.omit(data$number.of.SNPs)), main="Density Plot: number of snps", ylab="Frequency", 
     sub=paste("Skewness:", round(e1071::skewness(data$number.of.SNPs), 2)))  # density plot for 'speed'
polygon(density(na.omit(data$number.of.SNPs)), col="red")

plot(density(data$Day.difference), main="Density Plot: Distance", ylab="Frequency", 
     sub=paste("Skewness:", round(e1071::skewness(data$Day.difference), 2)))  # density plot for 'dist'
polygon(density(na.omit(data$Day.difference)), col="red")



# multicolinearity. do not include two var that are highly correlated
cor(data[,c("Breadth.of.Coverage...","Country","CT",
            "Day.difference","diarrhea",
            "Frequency.of.mutations..snps.day.","genotype",
            "number.of.SNPs")], use="complete")

# fit single regression
slinearMod <- lm(number.of.SNPs ~ Day.difference + Country, data=data)  # build linear regression model on full data
print(slinearMod); summary(slinearMod)
# we can see that singe linear model is not suitable for this
# let us make multiple linear model
str(data)

# glm
data$genotype <- as.factor(data$genotype)
data$genotype<- relevel(data$genotype,ref="GI.1") 
mlinearMod <- glm(genotype ~ number.of.SNPs, data=data)  # build linear regression model on full data
print(mlinearMod); summary(mlinearMod)

# fit multiple regression
mlinearMod <- lm(number.of.SNPs ~ Breadth.of.Coverage... + CT+ Day.difference + Country+MiseQ.run.batch., 
                 data=data)  # build linear regression model on full data
print(mlinearMod); summary(mlinearMod)

all_model1 <- 
  lm(number.of.SNPs~Breadth.of.Coverage...+Country+Day.difference+diarrhea+MiseQ.run.batch.,
     data=data)
print(summary(all_model))

all_model2 <- 
  lm(number.of.SNPs~Breadth.of.Coverage...+Country+Day.difference+diarrhea+MiseQ.run.batch.,
     data=data)
print(summary(all_model2))

oldpar <- par(oma=c(0,0,3,0), mfrow=c(2,2))
plot(all_model)
par(oldpar)

# we can see that there is a huge effect of the country. 
# it may be better to make to subset the data and fit models separately
# but, if use phil data separately, the sample size will be too low. 
# we need >10 observations for regression
# I will omit the data from the phillipines, because it constitute outlyer.
data_phil <- subset(data, subset=Country==1)
data_peru <- subset(data, subset=Country==2)

LinearModel_peru <- lm(number.of.SNPs ~ Day.difference + 
                      Breadth.of.Coverage... + MiseQ.run.batch. + CT, data=data_peru)
summary(LinearModel_peru)

#LinearModel_phil <- lm(number.of.SNPs ~ Day.difference, data_phil)

#summary(LinearModel_phil)





# NEW ANALYSES ------------------------------------------------------------

library(readxl)
data <- 
  readXL("C:/Users/Virology Tohoku/Dropbox/phil_data/PRIMAL_sample_conc_ms_saka_kte.xlsx",
         rownames=FALSE, header=TRUE, na="", sheet="ivar_whole_genome", 
         stringsAsFactors=TRUE)
data_peru <- subset(data, subset=country==2)
data_peru <- within(data_peru, {
  genotype <- factor(genotype, labels=c('GI.1','GI.2'))
})

require("ggplot2")
.df <- data.frame(x = data_peru$Day.difference, y = 
                    data_peru$number.of.SNPs, z = data_peru$genotype)
.plot <- ggplot(data = .df, aes(x = x, y = y, colour = z, shape = z)) + 
  geom_point() + 
  stat_smooth(aes(fill = z), method = "lm", se = FALSE) + 
  xlab("Day.difference") + 
  ylab("number.of.SNPs") + 
  labs(colour = "genotype", shape = "genotype", fill = "genotype") + 
  ggthemes::theme_base(base_size = 14, base_family = "sans") + 
  theme(legend.position = "right")
print(.plot)
rm(.df, .plot)

DataExplorer::create_report(data)
DataExplorer::create_report(peru_data)
library(dplyr)
str(data_peru)
data_peru %>% group_by(genotype) %>% 
  summarise(mean.snps = mean(number.of.SNPs),
            sum of SNPs = sum(number.of.SNPs),
            mean.day.difference = mean(Day.difference), 
            sum.of day.difference = sum(Day.difference), 
            n = n()) %>% 