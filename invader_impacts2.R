
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

#example of a 'species impact score'

#limit to woody species, (and/or pick a single dataset for ease of interpretation)
#set = "CVS" #Carolina Vegetation Survey
#dat2 = dat[dat$Dataset==set,]
dat2 = dat #use all data
str(dat2)
dat2$tree = regexpr("Tree",dat2$GrowthHabit)>0
dat2$shrub = regexpr("Shrub",dat2$GrowthHabit)>0
dat2$subshrub = regexpr("Subshrub",dat2$GrowthHabit)>0
dat2$woody = (dat2$tree+dat2$shrub+dat2$subshrub)>0
dat3 = dat2[dat2$woody==1,]
str(dat3)
  #still NAs cropping up... going to delete Plot=NAs for now...
dat3 = dat3[!is.na(dat3$Plot),]
dat3$PlotYear = paste(dat3$Plot,dat3$Year,sep="-") #unique observations

#Fagus grandifolia, doesn't play well with others (particularly other woodies)
sp = "NADO"
test = dat3[dat3$SpCode==sp,]
dim(test) #1681 occurrences
#plots wth sp
splots = unique(test$Plot)
comm = dat3[is.element(dat3$Plot,splots),]
str(comm)
#same dataset but target sp removed
comm2 = comm[comm$SpCode!=sp,]

#abundance (cover) of target species per plot
target.cover = data.frame(test$PlotYear,test$PctCov)
  names(target.cover) = c("PlotYear","TargetCov")

#'impact' measure, measured for all species except target species
#ignores differences in plot sizes
#example, total cover
tot.cover = tapply(comm2$PctCov,comm2$PlotYear,function(x)sum(x,na.rm=T))
hist(tot.cover)
#example, diversity (Shannon)
shannon = tapply(comm2$PctCov,comm2$PlotYear,function(x){
  p = x/(sum(x,na.rm=T))
  return(-sum(p*log(p)))
  })
hist(shannon)
#example, total richness
richness = tapply(comm2$PctCov,comm2$PlotYear,function(x)length(x))
impact = data.frame(rownames(tot.cover),tot.cover,shannon,richness)
  names(impact)[1] = "PlotYear"; rownames(impact) = NULL

#combine data by Plot
out = merge(target.cover,impact)
head(out)

#examine, plot simple linear relationships
library(MASS)
par(mfrow=c(1,3),mar=c(5,5,1,1))
plot(out$TargetCov,out$richness,xlab="Target Cover",ylab="Plot Richness")
  lm1 = lm(richness~TargetCov,out)
  lm2 = rlm(richness~TargetCov,out)
  abline(lm1,col="blue",lwd=2)
  abline(lm2,col="green",lwd=2)
plot(out$TargetCov,out$shannon,xlab="Target Cover",ylab="Plot Shannon Diversity")
  title(main=sp)
  lm1 = lm(shannon~TargetCov,out)
  lm2 = rlm(shannon~TargetCov,out)
  abline(lm1,col="blue",lwd=2)
  abline(lm2,col="green",lwd=2)
plot(out$TargetCov,out$tot.cover,xlab="Target Cover",ylab="Total Plot Cover")
  lm1 = lm(tot.cover~TargetCov,out)
  lm2 = rlm(tot.cover~TargetCov,out)
  abline(lm1,col="blue",lwd=2)
  abline(lm2,col="green",lwd=2)
  
#Frequency of all woodies
freq = table(comm$SpCode)
rev(sort(freq))

#compare to Cornus florida, a common understory species
sp = "COFL2"
test = dat3[dat3$SpCode==sp,]
dim(test) #1681 occurrences
#plots wth sp
splots = unique(test$Plot)
comm = dat3[is.element(dat3$Plot,splots),]
str(comm)
#same dataset but target sp removed
comm2 = comm[comm$SpCode!=sp,]

#abundance (cover) of target species per plot
target.cover = data.frame(test$PlotYear,test$PctCov)
names(target.cover) = c("PlotYear","TargetCov")

#'impact' measure, measured for all species except target species
#ignores differences in plot sizes
#example, total cover
tot.cover = tapply(comm2$PctCov,comm2$PlotYear,function(x)sum(x,na.rm=T))
hist(tot.cover)
#example, diversity (Shannon)
shannon = tapply(comm2$PctCov,comm2$PlotYear,function(x){
  p = x/(sum(x,na.rm=T))
  return(-sum(p*log(p)))
})
hist(shannon)
#example, total richness
richness = tapply(comm2$PctCov,comm2$PlotYear,function(x)length(x))
impact = data.frame(rownames(tot.cover),tot.cover,shannon,richness)
names(impact)[1] = "PlotYear"; rownames(impact) = NULL

#combine data by Plot
out = merge(target.cover,impact)
head(out)

#examine, plot simple linear relationships
par(mfrow=c(1,3),mar=c(5,5,1,1))
plot(out$TargetCov,out$richness,xlab="Target Cover",ylab="Plot Richness")
lm1 = lm(richness~TargetCov,out)
abline(lm1,col="blue",lwd=2)
plot(out$TargetCov,out$shannon,xlab="Target Cover",ylab="Plot Shannon Diversity")
title(main=sp)
lm1 = lm(shannon~TargetCov,out)
abline(lm1,col="blue",lwd=2)
plot(out$TargetCov,out$tot.cover,xlab="Target Cover",ylab="Total Plot Cover")
lm1 = lm(tot.cover~TargetCov,out)
abline(lm1,col="blue",lwd=2)


#calculate 'impact slopes' for all species found in at least 20 plots
sp.select = rownames(freq)[freq>19] #sample of 200
  #note that many species will never be that abundant, so slope can be poor estimate of impact in that case

spout = NULL
spout$SpCode = sp.select
spout$richness.slope = NULL
spout$diversity.slope = NULL
spout$cover.slope = NULL
spout$median.cover = NULL
spout$freq = NULL
spout$cover.t = NULL
spout$diversity.t = NULL
spout$richness.t = NULL

#species loop
for(i in 1:length(sp.select)) {
  print(i)

sp = sp.select[i]
test = dat3[dat3$SpCode==sp,]
splots = unique(test$Plot)
comm = dat3[is.element(dat3$Plot,splots),]
#same dataset but target sp removed
comm2 = comm[comm$SpCode!=sp,]
#abundance (cover) of target species per plot
target.cover = data.frame(test$PlotYear,test$PctCov)
names(target.cover) = c("PlotYear","TargetCov")
#'impact' measure, measured for all species except target species
#ignores differences in plot sizes
#example, total cover
tot.cover = tapply(comm2$PctCov,comm2$PlotYear,function(x)sum(x,na.rm=T))
#example, diversity (Shannon)
shannon = tapply(comm2$PctCov,comm2$PlotYear,function(x){
  p = x/(sum(x,na.rm=T))
  return(-sum(p*log(p)))
})
#example, total richness
richness = tapply(comm2$PctCov,comm2$PlotYear,function(x)length(x))
impact = data.frame(rownames(tot.cover),tot.cover,shannon,richness)
names(impact)[1] = "PlotYear"; rownames(impact) = NULL
#combine data by Plot
out = merge(target.cover,impact)

spout$richness.slope[i] = lm(richness~TargetCov,out)$coef[2]
spout$diversity.slope[i] = lm(shannon~TargetCov,out)$coef[2]
spout$cover.slope[i] = lm(tot.cover~TargetCov,out)$coef[2]
spout$median.cover[i] = median(out$TargetCov,na.rm=T)
spout$freq[i] = dim(out)[1]
spout$cover.t[i] = summary(lm(tot.cover~TargetCov,out))$coef[2,3]
spout$richness.t[i] = summary(lm(richness~TargetCov,out))$coef[2,3]
spout$diversity.t[i] = summary(lm(shannon~TargetCov,out))$coef[2,3]

}

#species information
species = data.frame("SpCode"=sp.select)
sp.info = unique(merge(species,taxa))
spout = merge(spout,sp.info,by="SpCode")

#transformations: inverse logit 
spout$cover.slope.inv = exp(spout$cover.slope)/(1+exp(spout$cover.slope))
spout$diversity.slope.inv = exp(spout$diversity.slope)/(1+exp(spout$diversity.slope))
spout$richness.slope.inv = exp(spout$richness.slope)/(1+exp(spout$richness.slope))

#relationships
plot(spout$median.cover,spout$cover.slope.inv,log="x",xlab="Focal Species Median % Cover",ylab="Associated Cover Score")
points(spout$median.cover[spout$NativeStatus=="I"],spout$cover.slope.inv[spout$NativeStatus=="I"],col="red",pch=19)
text(spout$median.cover[spout$NativeStatus=="I"],spout$cover.slope.inv[spout$NativeStatus=="I"]+.03,spout$SpCode[spout$NativeStatus=="I"],col="red",cex=.7)
abline(h=.5,col="gray",lty=2)
text(.1,.45,"Negative association")
text(.1,.55,"Positive association")
text(12,.1,"Native",cex=.9)
text(12,.05,"Introduced",cex=.9,col="red")

#linear slope not a good way to assess diversity impacts
plot(spout$median.cover,spout$diversity.slope)
plot(spout$median.cover,spout$richness.slope)


##3-1-23 addition, JDF
#can use regression t statistic (slope/SE) rather than raw slope to control for median cover issue in previous graph
  #t stat also much better than slope for diversity association statistics
  #t stats saved in above loop

par(mfrow=c(1,3))
plot(spout$median.cover,spout$cover.t,log="x",xlab="Focal Species Median % Cover",ylab="Associated Cover Score")
points(spout$median.cover[spout$NativeStatus=="I"],spout$cover.t[spout$NativeStatus=="I"],col="red",pch=19)
text(spout$median.cover,spout$cover.t+.7,spout$SpCode,col="red",cex=.4)
abline(h=0,col="gray",lty=2)
abline(h=c(-2,2),col="darkgray",lty=3)

plot(spout$median.cover,spout$diversity.t,log="x",xlab="Focal Species Median % Cover",ylab="Associated Shannon Score")
points(spout$median.cover[spout$NativeStatus=="I"],spout$diversity.t[spout$NativeStatus=="I"],col="red",pch=19)
text(spout$median.cover,spout$diversity.t+.7,spout$SpCode,col="red",cex=.4)
abline(h=0,col="gray",lty=2)
abline(h=c(-2,2),col="darkgray",lty=3)

plot(spout$median.cover,spout$richness.t,log="x",xlab="Focal Species Median % Cover",ylab="Associated Richness Score")
points(spout$median.cover[spout$NativeStatus=="I"],spout$richness.t[spout$NativeStatus=="I"],col="red",pch=19)
text(spout$median.cover,spout$richness.t+.7,spout$SpCode,col="red",cex=.4)
abline(h=0,col="gray",lty=2)
abline(h=c(-2,2),col="darkgray",lty=3)

#tot cover vs. shannon and richness scores: built-in correlation
par(mfrow=c(1,2))
plot(spout$cover.t,spout$diversity.t,log="",xlab="Associated Cover Score",ylab="Associated Shannon Score")
points(spout$cover.t[spout$NativeStatus=="I"],spout$diversity.t[spout$NativeStatus=="I"],col="red",pch=19)
abline(h=0,lty=2,col="gray")
abline(v=0,lty=2,col="gray")
abline(lsfit(spout$cover.t,spout$diversity.t),lwd=2,lty=2)
plot(spout$cover.t,spout$richness.t,log="",xlab="Associated Cover Score",ylab="Associated Richness Score")
points(spout$cover.t[spout$NativeStatus=="I"],spout$richness.t[spout$NativeStatus=="I"],col="red",pch=19)
abline(h=0,lty=2,col="gray")
abline(v=0,lty=2,col="gray")
abline(lsfit(spout$cover.t,spout$richness.t),lwd=2,lty=2)

#general graph
par(mfrow=c(1,1))
plot(1,1,col="white",xlim=c(-20,20),ylim=c(-20,20),xlab="Associated Cover Score",ylab="Residual Diversity Score")
abline(h=0,lty=2,col="gray")
abline(v=0,lty=2,col="gray")
arrows(0,0,0,10,length=.1,angle=45,lwd=2)
text(0,12,"More Diverse",cex=1.3)
arrows(0,0,0,-10,length=.1,angle=45,lwd=2)
text(0,-12,"Less Diverse",cex=1.3)
arrows(0,0,6,0,length=.1,angle=45,lwd=2)
text(14,0,"More Cover",cex=1.3)
arrows(0,0,-6,0,length=.1,angle=45,lwd=2)
text(-14,0,"Less Cover",cex=1.3)


#take diversity score as residual from relationship with cover
par(mfrow=c(1,2))
div.res = residuals(lm(spout$diversity.t~spout$cover.t))
plot(spout$cover.t,div.res,log="",xlab="Associated Cover Score",ylab="Residual Diversity Score",ylim=c(-15,10),xlim=c(-20,20))
abline(h=0,lty=2,col="gray")
abline(v=0,lty=2,col="gray")
points(spout$cover.t[spout$NativeStatus=="I"],div.res[spout$NativeStatus=="I"],col="red",pch=19)
#text(spout$cover.t,div.res+.2,spout$SpCode,col="gray",cex=.4)
text(spout$cover.t[spout$NativeStatus=="N"&abs(spout$cover.t)>5]+.3,div.res[spout$NativeStatus=="N"&abs(spout$cover.t)>5]+.5,spout$SpCode[spout$NativeStatus=="N"&abs(spout$cover.t)>5],col="gray",cex=.7)
text(spout$cover.t[spout$NativeStatus=="N"&abs(div.res)>5]+.3,div.res[spout$NativeStatus=="N"&abs(div.res)>5]+.5,spout$SpCode[spout$NativeStatus=="N"&abs(div.res)>5],col="gray",cex=.7)
text(spout$cover.t[spout$NativeStatus=="I"]+.3,div.res[spout$NativeStatus=="I"]+.5,spout$SpCode[spout$NativeStatus=="I"],col="red",cex=.7)
#arrows(0,-15,8,-15,angle=45,lwd=3,length=.1)
#arrows(0,-15,-8,-15,angle=45,lwd=3,length=.1)
#text(20,0,"Diversity Neutral")

#take richness score as residual from relationship with cover
div.res = residuals(lm(spout$richness.t~spout$cover.t))
plot(spout$cover.t,div.res,log="",xlab="Associated Cover Score",ylab="Residual Richness Score",ylim=c(-15,10),xlim=c(-20,20))
abline(h=0,lty=2,col="gray")
abline(v=0,lty=2,col="gray")
points(spout$cover.t[spout$NativeStatus=="I"],div.res[spout$NativeStatus=="I"],col="red",pch=19)
#text(spout$cover.t,div.res+.2,spout$SpCode,col="gray",cex=.4)
text(spout$cover.t[spout$NativeStatus=="N"&abs(spout$cover.t)>5]+.3,div.res[spout$NativeStatus=="N"&abs(spout$cover.t)>5]+.5,spout$SpCode[spout$NativeStatus=="N"&abs(spout$cover.t)>5],col="gray",cex=.7)
text(spout$cover.t[spout$NativeStatus=="N"&abs(div.res)>5]+.3,div.res[spout$NativeStatus=="N"&abs(div.res)>5]+.5,spout$SpCode[spout$NativeStatus=="N"&abs(div.res)>5],col="gray",cex=.7)
text(spout$cover.t[spout$NativeStatus=="I"]+.3,div.res[spout$NativeStatus=="I"]+.5,spout$SpCode[spout$NativeStatus=="I"],col="red",cex=.7)

