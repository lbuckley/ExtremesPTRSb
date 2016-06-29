#ANALYSIS FOR AUSTRALIA

library(msm) # for rtnorm
library(zoo) # for interpolating NAs

#------------------------------------------------------------------------
#Set up AUST PLOT

#Read Huey sites
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\Extremes\\Data\\BiologicalData\\Drosophila\\")
locs= read.csv("OzLocalitiesMel_17Sep2014.csv")

#FIND NEARBY ACORN SITES
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\Extremes\\Data\\ClimateData\\AustACORN\\")
acorn= read.csv("AustACORNstations.csv")
acorn= subset(acorn, acorn[,"Use"]=="y")
stats= acorn[,"Number"]

#==========================================================
### ANALYSIS
#=================================================================

#temp data
temp=  array(NA, dim=c(length(stats),1,15))
temp.baseline=  array(NA, dim=c(length(stats),1,15))

for(stat.k in 1:length(stats)){
#for(stat.k in 1:length(stats)){
print(stat.k)

filename= paste("acorn.sat.maxT.0", stats[stat.k], ".daily.txt", sep="")
dat.all= read.table(filename, na.strings="99999.9", skip=1)
names(dat.all)=c("YYYYMMDD","Tmax")

#RESTRICT TO YEARS
dat.all$Year= as.numeric(substr(dat.all$YYYYMMDD,1,4))
dat.all$Month= substr(dat.all$YYYYMMDD,5,6)
dat.all$YrMo= paste(dat.all$Year, dat.all$Month, sep="_")

dat= dat.all[dat.all$Year %in% 1962:2013,]

temp.dist= dat[,"Tmax"]
recent.ind.yr= which.max(dat[,"Year"]>1990)

#STORE TEMPERATURE DATA
temp.dist.cut= dat[dat$Year %in% 1991:2013,"Tmax"]
temp[stat.k,1,1:8]= c(mean(temp.dist.cut, na.rm=TRUE), sd(temp.dist.cut, na.rm=TRUE), median(temp.dist.cut, na.rm=TRUE), quantile(temp.dist.cut, probs=c(0.05, 0.95), na.rm=TRUE), range(dat[,"Year"]), length(unique(dat$Year)) )

#Generalized extreme value distribution
dat1= na.omit(temp.dist.cut)
try(mod.gev<- gev.fit(dat1, show=FALSE) ) #stationary
#mod.gev<- gev.fit.wrap(dat1, ydat=as.matrix(dat0$year), sigl=1, show=FALSE) # nonstationary 
#gev.diag(mod.gev)
if(class(mod.gev)!="try-error") temp[stat.k,1, 9]<-mod.gev$nllh
if(class(mod.gev)!="try-error") temp[stat.k,1, 10:12]<-mod.gev$mle #add another for non-stat
if(class(mod.gev)!="try-error") temp[stat.k,1, 14]<-mod.gev$conv #add another for non-stat

#Generalized pareto distribution, for number of times exceeds threshold
thresh= 35
mod.gpd <-gpd.fit(dat1, thresh, npy=365) #stationary
temp[stat.k,1, 15]<-mod.gpd$rate

#---------------------------------------
temp.dist.cut= dat[dat$Year %in% 1962:1990,"Tmax"]
temp.baseline[stat.k,1,1:8]= c(mean(temp.dist.cut, na.rm=TRUE), sd(temp.dist.cut, na.rm=TRUE), median(temp.dist.cut, na.rm=TRUE), quantile(temp.dist.cut, probs=c(0.05, 0.95), na.rm=TRUE), range(dat[,"Year"]), length(unique(dat$Year)) )

#Generalized extreme value distribution
dat1= na.omit(temp.dist.cut)
try(mod.gev<- gev.fit(dat1, show=FALSE) ) #stationary
#mod.gev<- gev.fit.wrap(dat1, ydat=as.matrix(dat0$year), sigl=1, show=FALSE) # nonstationary 
#gev.diag(mod.gev)
if(class(mod.gev)!="try-error") temp.baseline[stat.k,1, 9]<-mod.gev$nllh
if(class(mod.gev)!="try-error") temp.baseline[stat.k,1, 10:12]<-mod.gev$mle #add another for non-stat
if(class(mod.gev)!="try-error") temp.baseline[stat.k,1, 14]<-mod.gev$conv #add another for non-stat

#Generalized pareto distribution, for number of times exceeds threshold
thresh= 40
try(mod.gpd <-gpd.fit(dat1, thresh, npy=365) ) #stationary
if(class(mod.gpd)!="try-error") temp.baseline[stat.k,1, 15]<-mod.gpd$rate

} #end loop station

#=========================
## PLOT TOGETHER
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ExtremesPhilTrans\\figures\\")

file<-paste("AustGEV.pdf" ,sep="", collapse=NULL)
pdf(file,height = 8, width = 11)

sites=cbind(acorn,temp[,1,])
colnames(sites)[18:24]= c("gev.nllh", "gev.loc", "gev.scale", "gev.shape", "gev.mle4", "conv", "rate")

## PLOT SITE CLIMATE INFO
par(mfrow=c(2,2), cex=1.4, mar=c(2, 3, 1, 1), mgp=c(2, 1, 0), oma=c(2,0,0,0), lwd=2)

sites$st.lat= sites$Latitude
sites.coast= sites[sites$cline=="Coast",]
sites.cont= sites[sites$cline=="Cont",]

sites.temps= sites[order(sites$Latitude),]
temps.coast= sites.temps[sites.temps$cline=="Coast",]
temps.cont= sites.temps[sites.temps$cline=="Cont",]

#CLIMATE
plot(temps.coast$Latitude, temps.coast$gev.loc, col="black", type="b", xlim= range(sites.temps$Latitude), ylim= range(sites$gev.loc), ylab="GEV location", xlab="Latitude (°S)" )
points(temps.cont$Latitude, temps.cont$gev.loc, lty="dashed", col="grey", type="b")

legend("bottomright", legend=c("Coastal","Continental"), lty="solid", col=c("black","grey"), bty='n', cex=1)

plot(temps.coast$Latitude, temps.coast$gev.scale, col="black", type="b", xlim= range(sites.temps$Latitude), ylim= range(sites$gev.scale), ylab="GEV scale")
points(temps.cont$Latitude, temps.cont$gev.scale, lty="dashed", col="grey", type="b")

plot(temps.coast$Latitude, temps.coast$gev.shape, col="black", type="b", xlim= range(sites.temps$Latitude), ylim= range(sites$gev.shape), ylab="GEV shape")
points(temps.cont$Latitude, temps.cont$gev.shape, lty="dashed", col="grey", type="b")

plot(temps.coast$Latitude, temps.coast$rate, col="black", type="b", xlim= range(sites.temps$Latitude), ylim= range(sites$rate), ylab="Annual rate of exceeding 40°C")
points(temps.cont$Latitude, temps.cont$rate, lty="dashed", col="grey", type="b")

mtext("Latitude (°)", side=1, line = 0, cex=1.3, outer=TRUE)

dev.off()
