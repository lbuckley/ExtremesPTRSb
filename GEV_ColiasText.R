#CALCULATE EXTREMES ACROSS GEOGRAPHY
count= function(x) length(na.omit(x))

fdir= "C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ExtremesPhilTrans\\"

#library(sp)
library(ismev) #for gev
library(RAtmosphere)
library(insol)
library(date) #for julian functions
library(reshape)
library(maptools) #for mapping
library(evd) #for extremes value distributions
library(extRemes)
library(fExtremes) # generate gev

library(rnoaa) #http://recology.info/2015/07/weather-data-with-rnoaa/
library(dplyr)

library(RColorBrewer)
library(colorRamps)

library(zoo)
library(TeachingDemos) #for inset

#Wrapper for gev function
gev.fit.wrap= function(xdat, ydat, sigl=1, show=FALSE) try(gev.fit(xdat, ydat, sigl=1, show=FALSE), silent=TRUE)

#==================================================
#STATIONS
#STATION MAP: http://www.ncdc.noaa.gov/cdo-web/datatools/findstation

#Montrose: 38.4783° N, 107.8762° W
#GHCND:USC00055722
dat.co=ghcnd_search("USC00055722", var = "TMAX")

#Sacramento: 38.5816° N, 121.4944° W
#GHCND:USW00023271
#GHCND:USW00023232
dat.ca=ghcnd_search("USW00023271", var = "TMAX")

dat.co=dat.co$tmax
dat.ca=dat.ca$tmax

    #split date 
    date= as.Date(dat.co$date, "%Y-%m-$d")
    dat.co$year=as.numeric(format(date, "%Y"))
    dat.co$month=as.numeric(format(date, "%m"))
    dat.co$day=as.numeric(format(date, "%d"))
    
    date= as.Date(dat.ca$date, "%Y-%m-$d")
    dat.ca$year=as.numeric(format(date, "%Y"))
    dat.ca$month=as.numeric(format(date, "%m"))
    dat.ca$day=as.numeric(format(date, "%d"))
    
    ##Subset to summer
   dat.co= dat.co[dat.co$month %in% c(6,7,8,9),]
   dat.ca= dat.ca[dat.ca$month %in% c(6,7,8,9),]
    
    #Format data
    dat.co$tmax[which(dat.co$tmax==-9999)]= NA
    dat.co$tmax= dat.co$tmax/10 #correct for tenths of degrees or mm
    
    dat.ca$tmax[which(dat.ca$tmax==-9999)]= NA
    dat.ca$tmax= dat.ca$tmax/10 #correct for tenths of degrees or mm
    
#Matrix of parameters for extreme value distribution
ns.ext.stat= matrix(NA, 4, 9)
colnames(ns.ext.stat)= c("gev.nllh", "gev.loc", "gev.scale", "gev.shape", "gev.mle4", "se.loc", "se.scale", "se.shape", "rates")
rownames(ns.ext.stat)= c("CO.1961-1971", "CO.2001-2011", "CA.1961-1971", "CA.2001-2011")

#----------------------------
par(mfrow=c(1,2), cex=1.2, lwd=2, mar=c(1.5, 3, 1.5, 0.2), mgp=c(1.5, 0.5, 0), oma=c(2,0,0,0), bty="l")
labs=c("CO","CO","CA","CA")

for(stat.k in 1:4){

  if(stat.k==1) dat1= na.omit(dat.co[dat.co$year %in% 1961:1971,]$tmax)
  if(stat.k==2) dat1= na.omit(dat.co[dat.co$year %in% 2001:2011,]$tmax)
  if(stat.k==3) dat1= na.omit(dat.ca[dat.ca$year %in% 1961:1971,]$tmax)
  if(stat.k==4) dat1= na.omit(dat.ca[dat.ca$year %in% 2001:2011,]$tmax)
  
#Generalized extreme value distribution
try(mod.gev<- gev.fit(dat1, show=FALSE) ) #stationary
#mod.gev<- gev.fit.wrap(dat1, ydat=as.matrix(dat0$year), sigl=1, show=FALSE) # nonstationary 
#gev.diag(mod.gev)
if(class(mod.gev)!="try-error") ns.ext.stat[stat.k, 1]<-mod.gev$nllh
if(class(mod.gev)!="try-error") ns.ext.stat[stat.k, 2:4]<-mod.gev$mle #add another for non-stat
#if(class(mod.gev)!="try-error") ns.ext.stat[stat.k, 6]<-mod.gev$conv #add another for non-stat
if(class(mod.gev)!="try-error") ns.ext.stat[stat.k,6:8]<-mod.gev$se

#Generalized pareto distribution, for number of times exceeds threshold
if(stat.k %in% c(1,3)) thresh= quantile(dat1, 0.95)

mod.gpd <-gpd.fit(dat1, thresh, npy=122) #stationary
ns.ext.stat[stat.k, 9]<-mod.gpd$rate

#-----------------
#PLOT FIT
if(stat.k %in% c(1,3)) plot(density(dat1), main=labs[stat.k])
if(stat.k %in% c(2,4)) points(density(dat1), col="red", type="l")

if(class(mod.gev)!="try-error") if(mod.gev$mle[2]>0){
r = rgev(n = 1000, xi =mod.gev$mle[3], mu=mod.gev$mle[1], beta=mod.gev$mle[2])
#plot(r, type = "l", col = "steelblue") #xi = shape, mu = location, beta = scale,
points(density(r), type="l", col="grey", lty="dashed")
}


} #end loop stat

ns.ext.stat= as.data.frame(ns.ext.stat)

#compute CIs
ns.ext.stat$loc.5CI= ns.ext.stat$gev.loc -1.96*ns.ext.stat$se.loc
ns.ext.stat$loc.95CI= ns.ext.stat$gev.loc +1.96*ns.ext.stat$se.loc

ns.ext.stat$shape.5CI= ns.ext.stat$gev.shape -1.96*ns.ext.stat$se.shape
ns.ext.stat$shape.95CI= ns.ext.stat$gev.shape +1.96*ns.ext.stat$se.shape

ns.ext.stat$scale.5CI= ns.ext.stat$gev.scale -1.96*ns.ext.stat$se.scale
ns.ext.stat$scale.95CI= ns.ext.stat$gev.scale +1.96*ns.ext.stat$se.scale

