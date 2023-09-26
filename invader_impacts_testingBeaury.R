#Testing the binomial model of Beaury et al. with real and permuted data
#Fridley 9-26-23

##################
#create dataset

setwd("C:\\Users\\fridley\\OneDrive - Clemson University\\academic\\projects\\invasion_null_hyp\\impacts")

#datasets
spp = read.csv("SPCIS_plant_taxa.csv") #plot-spp dataset
plots = read.csv("SPCIS_plots.csv") #plot characteristics dataset

#plot-spp dataset
str(spp)
length(unique(spp$Plot)) #83,391 plots
length(unique(spp$AcceptedTaxonName)) #14,056 species
table(spp$RecordedStrata)
table(spp$NumberOfStrata)
#according to the metadata, there should always be a name in the "Original.TaxonName" column, but this is not the case:
sum(is.na(spp$Original.TaxonName))
#there are 43280 rows missing a species name. Investigate:
datNA = spp[is.na(spp$Original.TaxonName),]
head(datNA,100)
#they all seem to have cover values and are of various growth forms and durations
#need to consult the compilers; will remove for now
spp2 = spp[!is.na(spp$Original.TaxonName),]
dim(spp2)

#species-level dataset
#taxa = spp2[duplicated(spp$AcceptedTaxonName)==F,-1]
taxa = unique(spp2[,c(3,4,6,7,8)])
str(taxa)
table(taxa$NativeStatus)
#I     N    NI 
#1239 11943   359  
table(taxa$Duration)
table(taxa$GrowthHabit)

#plot-level dataset
str(plots)
length(unique(plots$Plot)) #same #plots as spp dataset
table(plots$Dataset)
table(plots$Year)
length(plots$Plot[duplicated(plots$Plot)]) #not many resamples

#merge to single dataset
dat = merge(spp2,plots,by.x=c("Plot","Year"),by.y=c("Plot","Year"))
str(dat)


#######################################################################
###A. empirical relationships (real data)

#calculate total number of native and nonnative species per plot
#limit to certain plot sizes (eg 1000 m2 and larger)
summary(dat$PlotArea.m2)
dat2 = dat[dat$PlotArea.m2>1000,]
Snat = tapply(dat2$NativeStatus,dat2$Plot,function(x)length(x[x=="N"]))
hist(log(Snat))
Sint = tapply(dat2$NativeStatus,dat2$Plot,function(x)length(x[x=="I"]))
hist(log(Sint)) #note mode is zero

#create new dataset that includes logs
df = data.frame(Snat,Sint,log.Snat=log(Snat),log.Sint=log(Sint),log.Snat1=log(Snat+1),log.Sint1=log(Sint+1))

#1. Snat and Sint
plot(df$Snat,df$Sint)
abline(lsfit(df$Snat,df$Sint))
summary(lm(Sint~Snat,df)) #strongly positive, slope is 0.19
mod = glm(cbind(Sint,Snat)~Snat,family="binomial",data=df) #binomial approach of Beaury et al.
summary(mod)
#no relationship, P=0.41

#2. log(Snat) and log(Sint), exclude -Inf (where Sint=0) where there are no aliens or no natives
df2 = df[df$Sint>0&df$Snat>0,]
plot(df2$log.Snat,df2$log.Sint)
abline(lsfit(df2$log.Snat,df2$log.Sint))
summary(lm(log.Sint~log.Snat,df2))
#strongly positive, slope is 0.20
mod = glm(cbind(log.Sint,log.Snat)~log.Snat,family="binomial",data=df2) #binomial approach of Beaury et al.
summary(mod)
#strongly negative, slope is -0.2

#3. log(Snat+1) and log(Sint+1), allowing use of all data
plot(df$log.Snat1,df$log.Sint1)
abline(lsfit(df$log.Snat1,df$log.Sint1))
summary(lm(log.Sint1~log.Snat1,df))
#strongly positive, slope is 0.22
mod = glm(cbind(log.Sint1,log.Snat1)~log.Snat1,family="binomial",data=df) #binomial approach of Beaury et al.
summary(mod)
#strongly negative, slope is -0.18
#actual model form of Beaury et al.:
mod = glm(cbind(Sint,Snat)~log(Snat+1),family="binomial",df)
summary(mod)
#even more strongly negative, slope is -0.29

#Summary: The native-alien richness relationship using plots > 1000 m is always positive using normal regression,
#with a slope of about 0.2 aliens per native species regardless of whether the data are log transformed or whether
#plots of zero aliens are ignored. Using the binomial glm that allows total richness to vary by plot, the relationship
#is either non-significant (with untransformed data) or negative (with log-transformed data)


#################################################
#B. Shuffled dataset 1 (breaking any real links between native and alien richness)

#example 1: shuffle without replacement entire NativeStatus column: ie, any individual can be N or I (and varies within sp)
#total frequency of N and I in the dataset remains constant, and plot species richness remains constant
dat2$NativeStatusShuffle = sample(dat2$NativeStatus,replace=F)
table(dat2$NativeStatus,dat2$NativeStatusShuffle)
Snat = tapply(dat2$NativeStatusShuffle,dat2$Plot,function(x)length(x[x=="N"]))
hist(log(Snat))
Sint = tapply(dat2$NativeStatusShuffle,dat2$Plot,function(x)length(x[x=="I"]))
hist(Sint,breaks=60)

#again, create new dataset that includes logs
df = data.frame(Snat,Sint,log.Snat=log(Snat),log.Sint=log(Sint),log.Snat1=log(Snat+1),log.Sint1=log(Sint+1))

#1. Snat and Sint
plot(df$Snat,df$Sint)
abline(lsfit(df$Snat,df$Sint))
summary(lm(Sint~Snat,df))
  #strongly positive, slope is 0.35
mod = glm(cbind(Sint,Snat)~Snat,family="binomial",data=df) #binomial approach of Beaury et al.
summary(mod)
  #strongly NEGATIVE, although slope is very small

#2. log(Snat) and log(Sint), exclude -Inf (where Sint=0) where there are no aliens or no natives
df2 = df[df$Sint>0&df$Snat>0,]
plot(df2$log.Snat,df2$log.Sint)
abline(lsfit(df2$log.Snat,df2$log.Sint))
summary(lm(log.Sint~log.Snat,df2))
#strongly positive, slope is 0.68
mod = glm(cbind(log.Sint,log.Snat)~log.Snat,family="binomial",data=df2) #binomial approach of Beaury et al.
summary(mod)
#strongly positive, slope is 0.17

#3. log(Snat+1) and log(Sint+1), allowing use of all data
plot(df$log.Snat1,df$log.Sint1)
abline(lsfit(df$log.Snat1,df$log.Sint1))
summary(lm(log.Sint1~log.Snat1,df))
#strongly positive, slope is 0.71
mod = glm(cbind(log.Sint1,log.Snat1)~log.Snat1,family="binomial",data=df) #binomial approach of Beaury et al.
summary(mod)
#strongly positive, slope is 0.09
#actual model form of Beaury et al.:
mod = glm(cbind(Sint,Snat)~log(Snat+1),family="binomial",df)
summary(mod)
#strongly NEGATIVE, slope is -0.08

#To summarize, there remain artifacts in native-alien richness relationships no matter which way you treat the 
#data statistically. OLS regression always finds a positive relationship in shuffled data regardless of transformation, 
#albeit one much stronger than than the real relationship. The binomial GLM result is always significant with shuffled 
#data but can be positive or negative depending on the type of transformation done; Beaury's method finds a negative
#relationship in shuffled data

#######################################################
#C. Shuffled dataset #2: restrict shuffling to across species (keep species frequency and plot membership the same, 
#but shuffle its native identity, keeping total N-I frequencies the same)
spp2$NativeStatusShuffle = sample(spp2$NativeStatus,replace=F) #shuffle first in species db
table(spp2$NativeStatusShuffle,spp2$NativeStatus) 
dat3 = merge(spp2,plots,by.x=c("Plot","Year"),by.y=c("Plot","Year")) #then add to plot db
str(dat3)
table(dat3$NativeStatusShuffle,dat3$NativeStatus) 

Snat = tapply(dat3$NativeStatusShuffle,dat3$Plot,function(x)length(x[x=="N"]))
hist(log(Snat))
Sint = tapply(dat3$NativeStatusShuffle,dat3$Plot,function(x)length(x[x=="I"]))
hist(Sint,breaks=60)

#again, create new dataset that includes logs
df = data.frame(Snat,Sint,log.Snat=log(Snat),log.Sint=log(Sint),log.Snat1=log(Snat+1),log.Sint1=log(Sint+1))

#1. Snat and Sint
plot(df$Snat,df$Sint)
abline(lsfit(df$Snat,df$Sint))
summary(lm(Sint~Snat,df))
#strongly positive, slope is 0.1
mod = glm(cbind(Sint,Snat)~Snat,family="binomial",data=df) #binomial approach of Beaury et al.
summary(mod)
#strongly NEGATIVE, although slope is very small

#2. log(Snat) and log(Sint), exclude -Inf (where Sint=0) where there are no aliens or no natives
df2 = df[df$Sint>0&df$Snat>0,]
plot(df2$log.Snat,df2$log.Sint)
abline(lsfit(df2$log.Snat,df2$log.Sint))
summary(lm(log.Sint~log.Snat,df2))
#strongly positive, slope is 0.58
mod = glm(cbind(log.Sint,log.Snat)~log.Snat,family="binomial",data=df2) #binomial approach of Beaury et al.
summary(mod)
#strongly positive, slope is 0.45

#3. log(Snat+1) and log(Sint+1), allowing use of all data
plot(df$log.Snat1,df$log.Sint1)
abline(lsfit(df$log.Snat1,df$log.Sint1))
summary(lm(log.Sint1~log.Snat1,df))
#strongly positive, slope is 0.59
mod = glm(cbind(log.Sint1,log.Snat1)~log.Snat1,family="binomial",data=df) #binomial approach of Beaury et al.
summary(mod)
#strongly positive, slope is 0.31
#actual model form of Beaury et al.:
mod = glm(cbind(Sint,Snat)~log(Snat+1),family="binomial",df)
summary(mod)
#strongly NEGATIVE, slope is -0.05

#summary: conclusion is the same regardless of whether shuffling occurs at the scale of the entire dataset or is done
#separately by species--there are built-in relationships between native and alien richness in randomized data that
#change depending on log transformation and type of regression used

#note: Beaury's mixed model version:
library(glmmTMB) 
model <- glmmTMB(cbind(Non_richness, N_richness) ~ ln_N_richness_Z*(Ecoregion+Community) + ln_Dist_landcover_m_Z + HZ_Z +  (1|Park), ziformula = ~ln_N_richness_Z + Community + (1|ln_Dist_landcover_m_Z) + (1|HZ_Z) + (1|Ecoregion) + (1|Park), data=plots, family="binomial",control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3), profile = TRUE))
summary(model)

