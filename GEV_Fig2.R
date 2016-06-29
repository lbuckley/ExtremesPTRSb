#CALCULATE EXTREMES ACROSS GEOGRAPHY
count= function(x) length(na.omit(x))

fdir= "C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ExtremesPhilTrans\\"

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

#Wrapper for gev function
gev.fit.wrap= function(xdat, ydat, sigl=1, show=FALSE) try(gev.fit(xdat, ydat, sigl=1, show=FALSE), silent=TRUE)

#==================================================
#Choose Variables 
var="TMAX"

#===================================
#DOWNLOAD DATA

#1: Set sites
#Examine extreme value for sites along longitude line
#NA 100W and SA 65W
noa= cbind( seq(20,65,1), rep(-100,length(seq(20,65,1))) )
asia= cbind( seq(10,75,1), rep(77.5,length(seq(10,75,1))) )

sites= rbind(noa, asia)
colnames(sites)= c("Lat", "Lon")
sites= as.data.frame(sites)
sites= sites[order(sites$Lon,sites$Lat),]

#===================================
#Read GGCN database inventory http://www.ncdc.noaa.gov/oa/climate/ghcn-daily/
setwd(paste(fdir,"data\\",sep=""))
stations=read.table("ghcnd-inventory.txt")
names(stations)=c("Id","Lat","Lon","Var","FirstYear","LastYear")
stations.un= unique(stations$Id)
stats.elev= read.csv("ghcnd-stations.csv")

#Restrict stations to those with recent data
stations= stations[which(stations$LastYear>2010),] #stations$FirstYear<1996 & 

#choose var
stations= stations[which(stations$Var==var),]

#Restrict elevations
##get elevations
stations$elev= stats.elev[match(stations$Id, stats.elev$Id),"Elev"]
stations= stations[which(stations$elev<500),]

#------------
stat.coords= cbind(stations$Lon, stations$Lat)

#Matrix of parameters for extreme value distribution
ns.ext.stat= matrix(NA, nrow(sites), 10)
colnames(ns.ext.stat)= c("gev.nllh", "gev.loc", "gev.scale", "gev.shape", "gev.mle4", "gev.conv2", "gev.loc2", "gev.scale2", "gev.shape2", "gev.conv2")

#Placeholders for extremes and medians
sites$Ta.99=NA

#---------------------------
sites$Id=NA; sites$st.lat=NA; sites$st.lon=NA

#add continent
sites$cont=NA
sites$cont[which(sites$Lon==-100)]="noa"
sites$cont[which(sites$Lon==77.5)]="asia"

#============================================
#start analysis
#============================================
#Calculate distribution

#set up colors by continent
conts= c("asia", "noa")  

sites$col=NA
for(con in 1:length(conts) ){
  csites= which(sites[,"cont"]==conts[con])
  sites[csites,"col"]= 1:length(csites)
}

cols.lat=blue2green(max(sites[,"col"]))

#set up plots
setwd(paste(fdir,"figures\\",sep=""))
file<-paste(var, "GEV_fullyear.pdf" ,sep="", collapse=NULL)
pdf(file,height = 11, width = 8)

#set up density plots
par(mfrow=c(4,4), cex=1.2, lwd=2, mar=c(1, 1, 1.5, 0.2), mgp=c(1.5, 0.5, 0), oma=c(2,2,0,0), bty="l")

for(stat.k in 1:nrow(sites) ){  #1:nrow(sites)
  print(stat.k)
  
  min.dist<- order(spDistsN1(stat.coords, as.numeric(sites[stat.k,c("Lon","Lat")]), longlat = TRUE))[1:100]
  min.site= stations$Id[min.dist]
  drop.site= which(min.site %in% sites$Id)
  if(length(drop.site)>0) min.site= min.site[-drop.site]
  
  ind=0
  years=NA
  while(length(years)<10){
    ind= ind + 1
    dat=ghcnd_search(min.site[ind], var = var)
    while(is.null(dat)){ ind=ind+1
    dat=ghcnd_search(min.site[ind], var = var)
    }
    
    if(var=="TMAX") dat= dat$tmax
    if(var=="PRCP") dat= dat$prcp
    names(dat)[2]="var" 
    
    #split date 
    date= as.Date(dat$date, "%Y-%m-$d")
    dat$year=as.numeric(format(date, "%Y"))
    dat$month=as.numeric(format(date, "%m"))
    dat$day=as.numeric(format(date, "%d"))
    
    ##Subset to summer
   nlat= sites$Lat[stat.k]>0
    #in northern lat
  if(nlat) dat= dat[dat$month %in% c(6,7,8),]
    #if southern lat
  if(!nlat) dat= dat[dat$month %in% c(1,2,12),]
    
    #Format data
    dat$var[which(dat$var==-9999)]= NA
    dat$var= dat$var/10 #correct for tenths of degrees or mm
    
    ## FIND YEARS WITH NEARLY COMPLETE DATA
    dat.agg= aggregate(dat, list(dat$year),FUN=count)
    years= dat.agg$Group.1[which(dat.agg$var>80)]
    dat= dat[which(dat$year %in% years),]
  } #end check data
  
  #RECORD SITE DATA
  sites$Id[stat.k]= as.character(min.site[ind])
  sites$st.lat[stat.k]= stations$Lat[min.dist[ind]]
  sites$st.lon[stat.k]= stations$Lon[min.dist[ind]]
  
#-----------------------
#ANALYZE EXTREMES 

dat$var[dat$var< (-50)]=NA
dat0=dat[!is.na(dat$var),]
dat1=dat0$var
dat= dat$var

if(nrow(dat0)>0){ #check data

#calculate quantiles
sites$Ta.99[stat.k]= quantile(dat, na.rm = TRUE, probs=0.99)

#Generalized extreme value distribution
try(mod.gev<- gev.fit(dat1, show=FALSE) ) #stationary
#mod.gev<- gev.fit.wrap(dat1, ydat=as.matrix(dat0$year), sigl=1, show=FALSE) # nonstationary 
if(class(mod.gev)!="try-error") ns.ext.stat[stat.k, 1]<-mod.gev$nllh
if(class(mod.gev)!="try-error") ns.ext.stat[stat.k, 2:4]<-mod.gev$mle #add another for non-stat
if(class(mod.gev)!="try-error") ns.ext.stat[stat.k, 6]<-mod.gev$conv #add another for non-stat
#if(class(mod.gev)!="try-error") ns.ext.stat[stat.k,6:9]<-mod.gev$se

#-----------------
#PLOT FIT
plot(density(dat1), main=sites$Id[stat.k] )

if(class(mod.gev)!="try-error") if(mod.gev$mle[2]>0){
r = rgev(n = 1000, xi =mod.gev$mle[3], mu=mod.gev$mle[1], beta=mod.gev$mle[2])
points(density(r), type="l", col="red")
}

} #check dat
} #end loop stat

dev.off()
#---------------------------
#Add info and write out
ns.ext.stat= as.data.frame(ns.ext.stat)
sites= cbind(sites, ns.ext.stat)

#---------------------------
#set up plots
setwd(paste(fdir,"figures\\",sep=""))
file<-file<-paste("Fig2_",var,"ByLat.pdf" ,sep="", collapse=NULL)
pdf(file,height = 8, width = 11)

#5. SUMMARY PLOTS
par(mfrow=c(2,3), cex=1.2, lwd=2, mar=c(1.5, 3, 1.5, 0.2), mgp=c(1.5, 0.5, 0), oma=c(2,0,0,0), bty="l")

#---------------------------------
#PLOT EXAMPLE GEVs

rcol= gray.colors(3, end=0.9)

r1= rgev(n = 50000, xi = 0, mu=30, beta=3)
r2= rgev(n = 50000, xi = 0, mu=50, beta=3)
r3= rgev(n = 50000, xi = 0, mu=30, beta=4.5)
r4= rgev(n = 50000, xi = 0, mu=50, beta=4.5)
r5= rgev(n = 50000, xi = 0, mu=30, beta=6)
r6= rgev(n = 50000, xi = 0, mu=50, beta=6)
plot(density(r1), col=rcol[1], type="l", main="",xlab="Maximum daily temperature (°C)", xlim=range(18,95))
points(density(r2), col=rcol[1], type="l", lty="dashed")
points(density(r3), col=rcol[2], type="l")
points(density(r4), col=rcol[2], type="l", lty="dashed")
points(density(r5), col=rcol[3], type="l")
points(density(r6), col=rcol[3], type="l", lty="dashed")

legend(55,0.12, lty="solid", col=rcol, cex=0.8,legend=c("scale=3", "scale=4.5", "scale=6"), bty="n")
legend(55,0.09, lty=c("solid","dashed"), col="black", cex=0.8,legend=c("location=30", "location=50"), bty="n")

#---------------------
#PLOT LAT PATTERNS
#divide by continents
conts= unique(sites$cont)
#cols= rainbow(length(conts))
cols= c("black","gray")

for(cn in 1:length(conts)){
sites1= subset(sites, sites$cont==conts[cn])
sites1= subset(sites1, exp(sites1$gev.scale)>0 & sites1$gev.loc>0)
if(cn==1) plot(sites1$Lat, sites1$Ta.99, ylab="99% quantile", xlab="", main="", col=cols[cn], type="p", xlim=range(sites$st.lat), ylim=range(18,47.2))  #ylim=range(sites$Ta.99, na.rm=TRUE)
if(cn>1) points(sites1$Lat, sites1$Ta.99, col=cols[cn], type="p")

if(cn==1) legend("bottomleft", col=cols, cex=1,pch=1,legend=c("North America","Asia"), bty="n")
}

for(cn in 1:length(conts)){
sites1= subset(sites, sites$cont==conts[cn])
sites1= subset(sites1, sites1$gev.scale>0 & sites1$gev.loc>0)
if(cn==1) plot(sites1$Lat, sites1$gev.loc, ylab="GEV location", xlab="", col=cols[cn], type="p", xlim=range(sites$st.lat), ylim=range(sites$gev.loc, na.rm=TRUE))
if(cn>1) points(sites1$Lat, sites1$gev.loc, col=cols[cn], type="p")
}

#---------------------------
#EXAMPLE GEV
rcol= gray.colors(4, end=0.9)

r1= rgev(n = 50000, xi = -0.3, mu=40, beta=4.5)
r2= rgev(n = 50000, xi = -0.1, mu=40, beta=4.5)
r3= rgev(n = 50000, xi = 0.0, mu=40, beta=4.5)
r4= rgev(n = 50000, xi = 0.1, mu=40, beta=4.5)
plot(density(r1), col=rcol[1], type="l", main="",xlab="", xlim=range(25,65))
points(density(r2), col=rcol[2], type="l")
points(density(r3), col=rcol[3], type="l")
points(density(r4), col=rcol[4], type="l")

legend("topright", lty="solid", col=rcol, cex=0.8,legend=c("shape= -0.3","shape= -0.1","shape= 0.0","shape= 0.1"), bty="n")

#---------------------------

for(cn in 1:length(conts)){
  sites1= subset(sites, sites$cont==conts[cn])
  sites1= subset(sites1, sites1$gev.scale>0 & sites1$gev.loc>0)
  if(cn==1) plot(sites1$Lat, sites1$gev.scale, ylab="GEV scale", xlab="", col=cols[cn], type="p", xlim=range(sites$st.lat), ylim=range(sites$gev.scale, na.rm=TRUE))
  if(cn>1) points(sites1$Lat, sites1$gev.scale, col=cols[cn], type="p")
}

for(cn in 1:length(conts)){
sites1= subset(sites, sites$cont==conts[cn])
sites1= subset(sites1, sites1$gev.scale>0 & sites1$gev.loc>0)
if(cn==1) plot(sites1$Lat, sites1$gev.shape, ylab="GEV shape", xlab="", col=cols[cn], type="p", xlim=range(sites$st.lat), ylim=range(sites$gev.shape, na.rm=TRUE)) 
if(cn>1) points(sites1$Lat, sites1$gev.shape, col=cols[cn], type="p")
}

mtext(c("Maximum daily temperature (°C)","Latitude (°)"), at=c(0.2,0.7), side=1, line=0.5, outer=TRUE, cex=1.2)

dev.off() #end output to pdf
