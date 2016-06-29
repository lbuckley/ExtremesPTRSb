#CALCULATE EXTREMES ACROSS TIME SCALES
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
library(zoo) #for rollmean
library(RColorBrewer)
library(colorRamps)
library(climdex.pcic)

#Wrapper for gev function
gev.fit.wrap= function(xdat, ydat, sigl=1, show=FALSE) try(gev.fit(xdat, ydat, sigl=1, show=FALSE), silent=TRUE)
julian.wrap= function(date) julian(date, origin=as.Date(paste(format(date, "%Y"),"-01-01",sep="")))


#-----------------------------
#Read in site data
setwd(paste(fdir,"out\\",sep=""))
temp.dat= read.csv("TMAXExtOut.csv")

#select sites
temp.dat= temp.dat[c(1,6),1:14]

temp.dat$var="TMAX"

sites= rbind(temp.dat)

#Set up data storage
rl= array(NA, dim=c(nrow(sites),8,3))
ns.ext.stat= matrix(NA, nrow(sites), 10)
colnames(ns.ext.stat)= c("gev.nllh", "gev.loc", "gev.scale", "gev.shape", "gev.mle4", "gev.conv2", "gev.loc2", "gev.scale2", "gev.shape2", "gev.conv2")
#--------------------------------

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

#============================================
#Calculate distribution

#set up plots
setwd(paste(fdir,"figures\\",sep=""))
file<-paste("Fig3_TimeScales_sim.pdf" ,sep="", collapse=NULL)
pdf(file,height = 11, width = 8)

#set up density plots
par(mfrow=c(3,2), cex=1.2, lwd=2, mar=c(2.5, 1, 1.5, 0.2), mgp=c(1.5, 0.5, 0), oma=c(0,2,0,0), bty="l")
cont.last="mars"

for(stat.k in 1:2 ){  
  var= sites$var[stat.k]
  min.dist<- order(spDistsN1(stat.coords, as.numeric(sites[stat.k,c("Lon","Lat")]), longlat = TRUE))[1:100]
  min.site= stations$Id[min.dist]
  
  dat=ghcnd_search(sites$Id[stat.k], var = var)
    
    if(var=="TMAX") dat= dat$tmax
    if(var=="PRCP") dat= dat$prcp
    names(dat)[2]="var" 
    
    #split date 
    date= as.Date(dat$date, "%Y-%m-$d")
    dat$year=as.numeric(format(date, "%Y"))
    dat$month=as.numeric(format(date, "%m"))
    dat$day=as.numeric(format(date, "%d"))
    dat$J= sapply(date, FUN=julian.wrap)
        
    ##Subset to summer
    nlat= sites$Lat[stat.k]>0
    #in northern lat
    if(nlat) dat= dat[dat$month %in% c(6,7,8),]
    #if southern lat
    if(!nlat) dat= dat[dat$month %in% c(1,2,12),]
    
    #Format data
    dat$var[which(dat$var==-9999)]= NA
    if(var=="TMAX") dat$var= dat$var/10 #correct for tenths of degrees
    #sort by date
    dat= dat[order(dat$date),]
    
    ## FIND YEARS WITH NEARLY COMPLETE DATA
    dat.agg= aggregate(dat, list(dat$year),FUN=count)
    years= dat.agg$Group.1[which(dat.agg$var>80)]
    dat= dat[which(dat$year %in% years),]
  
  #-----------------------
    
    dat.week= rollmean(na.approx(dat$var), 7) #!Change to better handle NAs
    dat.month= aggregate(dat, list(dat$month,dat$year),FUN=mean, na.rm=TRUE)
    dat.year= aggregate(dat, list(dat$year),FUN=mean, na.rm=TRUE) #!change tomake sure enough years
    
    #---------------------------------
    #Return times
    
    dat0=dat[!is.na(dat$var),]
    dat1=dat0$var
    #CLEAN DATA
    #dat1= dat1[which(dat1>-10)]
    
    #Generalized extreme value distribution
    try(mod.gev<- gev.fit(dat1, show=FALSE) ) #stationary
    #mod.gev<- gev.fit.wrap(dat1, ydat=as.matrix(dat0$year), sigl=1, show=FALSE) # nonstationary 
    #gev.diag(mod.gev)
    if(class(mod.gev)!="try-error") ns.ext.stat[stat.k, 1]<-mod.gev$nllh
    if(class(mod.gev)!="try-error") ns.ext.stat[stat.k, 2:4]<-mod.gev$mle #add another for non-stat
    if(class(mod.gev)!="try-error") ns.ext.stat[stat.k, 6]<-mod.gev$conv #add another for non-stat
    
    #Generalized pareto distribution, for number of times exceeds threshold
    #mod.gpd <-gpd.fit(Ta$value, 40, npy=92) #stationary
    ## nonstationary 
    try(mod.gpd<-gpd.fit(dat1, 40, npy=92, ydat=as.matrix(dat0$year), sigl=1),silent = FALSE) 
  
    #RETURN LEVELS:  MLE Fitting of GPD - package extRemes
    thresh= quantile(dat.month$var, 0.9, na.rm=TRUE)
    mpers= c(2,5,10,20,50,100)
    for(m in 1:length(mpers)){
    pot.day= fpot(dat0$var, threshold=30, npp=365.25, mper=mpers[m] )
    pot.week= fpot(dat.week, threshold=30, npp=365.25, mper=mpers[m] )
    if(m>4)pot.month= fpot(dat.month$var, threshold=30, npp=3, mper=mpers[m] )
    #plot(plot.day)
    
    rl[stat.k, m, 1]=pot.day$estimate[1]
    rl[stat.k, m, 2]=pot.week$estimate[1]
    if(m>4)rl[stat.k, m, 3]=pot.month$estimate[1]
    }
    rl[stat.k, 7, 1]=pot.day$pat
    rl[stat.k, 7, 2]=pot.week$pat
    rl[stat.k, 7, 3]=pot.month$pat
    
    #-------------------------------------
    #simulate intervals
    r = rgev(n = 100000, xi =ns.ext.stat[stat.k,4], mu=ns.ext.stat[stat.k,2], beta=ns.ext.stat[stat.k,3])
    
    #single heat waves
    thresh= quantile(r, 0.90)
    dat.thresh= rep(0, length(r))
    dat.thresh[which(r>thresh)]=1
    dat.rl= rle(dat.thresh)
    dat.rl= dat.rl$lengths[dat.rl$values==0]
    #plot(density(dat.rl))
    hw.list[[stat.k+2]]=dat.rl
    
    thresh= quantile(r, 0.95)
    dat.thresh= rep(0, length(r))
    dat.thresh[which(r>thresh)]=1
    dat.rl= rle(dat.thresh)
    dat.rl= dat.rl$lengths[dat.rl$values==0]
    #points(density(dat.rl), type="l")
    hw.list[[stat.k+4]]=dat.rl
    
    thresh= quantile(r, 0.98)
    dat.thresh= rep(0, length(r))
    dat.thresh[which(r>thresh)]=1
    dat.rl= rle(dat.thresh)
    dat.rl= dat.rl$lengths[dat.rl$values==0]
    #points(density(dat.rl), type="l")
    hw.list[[stat.k+6]]=dat.rl
    #-----------------------
  #DENSITY PLOTS
    d.day= density( na.omit(dat$var)) 
    d.week=  density( na.omit(dat.week)) 
    d.month=  density( na.omit(dat.month$var)) 
    d.year=  density( na.omit(dat.year$var)) 
    
    cols=blue2green(4)
    #cols=brewer.pal(n = 5, name = "YlGnBu")
    
    m.labs=c("North America, 24°N","North America, 45°N")
    
    plot(d.day, xlim=range(13, 45), ylim=range(0,0.25), main=m.labs[stat.k], xlab="Maximum temperature (°C)", ylab="", cex.main=1, col=cols[1]) 
    points(d.week, col=cols[2], type="l")
    points(d.month, col=cols[3], type="l")
    points(d.year, col=cols[4], type="l")
    
    leg.lab=c("day","week","month","year")
    if(stat.k==1) legend("topleft", lty="solid",col=cols, legend=leg.lab, bty="n")
     
    #-----------------------------
    
} #end loop stat

#-----------------------

#PLOT RETURN LEVELS
rp= c(2,5,10,20,50,100)

for(k in 1:2){
 labs= "Temperature (°C)"
 
plot(rp, rl[k, 1:6, 1], type="l", ylim=range(na.omit(rl[k,1:6,])), col=cols[1], ylab=labs, xlab="Return Period (years)", main="", xlim=c(0,80))
points(rp, rl[k, 1:6, 2], type="l", col=cols[2] )
points(rp[4:6], rl[k, 4:6, 3], type="l", col=cols[3] )
}

#----------------------------
#PLOT SIMULATED

plot(density(hw.list[[3]]), col=cols[1], type="l", main="", xlim=range(0,200), ylim= range(0,0.082), xlab= "Days between heat events", lty="dotted") # xlim=range(0,60), ylim= range(0,0.1), 

points(density(hw.list[[5]]), col=cols[1], type="l", lty="dashed") 
points(density(hw.list[[7]]), col=cols[2], type="l") 

legend("topright", lty=c("dotted","dashed","solid"),col=cols[1], legend=c("90 %ile", "95 %ile","98 %ile"), bty="n")

#------
plot(density(hw.list[[4]]), col=cols[1], type="l", main="", xlim=range(0,200), ylim= range(0,0.082), xlab= "Days between heat events", lty="dotted") # xlim=range(0,60), ylim= range(0,0.1), 

points(density(hw.list[[6]]), col=cols[1], type="l", lty="dashed") 
points(density(hw.list[[8]]), col=cols[2], type="l") 

mtext(c("Density","Maximum temperature (°C)", "Density"), side = 2, line = 1, outer = TRUE, cex=1.5, at=c(0.85,0.5,0.15))

dev.off()
