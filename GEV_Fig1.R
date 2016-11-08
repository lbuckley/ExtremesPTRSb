## Code from Kingsolver JG and Buckley LB. 2017. Quantifying thermal extremes and biological variation to predict evolutionary responses to changing climate. Phil. Trans. R. Soc. B. 
## Code for Figure 1
##
###########################################

#CALCULATE EXTREMES ACROSS GEOGRAPHY
count= function(x) length(na.omit(x))

#set base directory
fdir= "C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ExtremesPhilTrans\\"

#library(sp)
library(ismev) #for gev
library(RAtmosphere)
#see also package evir
library(insol)
#library(GhcnDaily)
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
library(plotrix)

#Wrapper for gev function
gev.fit.wrap= function(xdat, ydat, sigl=1, show=FALSE) try(gev.fit(xdat, ydat, sigl=1, show=FALSE), silent=TRUE)

#Make wrapper for Julian function
julian.wrap=function(date1){ 
  #date1= as.Date(date1, "%Y-%m-$d")
  origin1 = as.Date(paste("1-1-",as.numeric(format(date1, "%Y")), sep=""), format="%m-%d-%Y")
  julian(date1, origin=origin1)
}

#==================================================
#Choose Variables
var="TMAX"

#------------------------
#DOWNLOAD DATA

#1: Set sites
#Examine extreme value for sites along longitude line
#NA 100W and SA 65W
noa= cbind( seq(20,65,5), rep(-100,length(seq(20,65,5))) )
asia= cbind( seq(10,75,5), rep(77.5,length(seq(10,75,5))) )

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
stations= stations[which(stations$LastYear>2010),] 

#choose var
stations= stations[which(stations$Var==var),]

#Restrict elevations
##get elevations
stations$elev= stats.elev[match(stations$Id, stats.elev$Id),"Elev"]
stations= stations[which(stations$elev<500),]

#------------
stat.coords= cbind(stations$Lon, stations$Lat)

#Matrix of parameters for extreme value distribution
ns.ext.stat= matrix(NA, nrow(sites), 17)
colnames(ns.ext.stat)= c("gev.nllh", "gev.loc", "gev.scale", "gev.shape", "gev.mle4", "gev.loc", "gev.shape", "gev.scale", "gev.se4", "gpd.nllh", "gpd.rate", "gpd.mle1", "gpd.mle2", "gpd.mle3", "gev.se1", "gev.se2", "gev.se3")

#Placeholders for extremes and medians
sites$Ta.99=NA
sites$Ta.95=NA
sites$Ta.50=NA
sites$Ta99.50=NA
sites$Trange=NA

#---------------------------
sites$Id=NA; sites$st.lat=NA; sites$st.lon=NA

#Plot sites
data(wrld_simpl)
plot(wrld_simpl, xlim=c(-180,180), ylim=c(-60,60))
points(sites$st.lon, sites$st.lat, col="red", cex=1)

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
file<-paste(var,"Hist_1980_rb.pdf" ,sep="", collapse=NULL)
pdf(file,height = 8, width = 11)

#set up density plots
par(mfrow=c(1,2), cex=1.2, lwd=2, mar=c(1, 1, 1.5, 0.2), mgp=c(3, 1, 0), oma=c(2,2,0,0), bty="l")
cont.last="mars"

for(stat.k in 1:nrow(sites) ){  #1:nrow(sites)
  min.dist<- order(spDistsN1(stat.coords, as.numeric(sites[stat.k,c("Lon","Lat")]), longlat = TRUE))[1:100]
  min.site= stations$Id[min.dist]
  
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
    dat$j= unlist(lapply(date, FUN="julian.wrap"))
    
    ##Subset to summer
    nlat= sites$Lat[stat.k]>0
    #in northern lat
    if(nlat) dat= dat[dat$month %in% c(6,7,8),]
    #if southern lat
    if(!nlat) dat= dat[dat$month %in% c(1,2,12),]
    
    #subset to last 30 years
    dat= dat[dat$year>1979,]
    
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
  sites$st.elev[stat.k]= stations$elev[min.dist[ind]]
  
  #-----------------------
  #ANALYZE EXTREMES 
  
  dat0=dat[!is.na(dat$var),]
  dat1=dat0$var
  dat= dat$var
  
  #calculate quantiles
  sites$Ta.99[stat.k]= quantile(dat, na.rm = TRUE, probs=0.99)
  sites$Ta.95[stat.k]= quantile(dat, na.rm = TRUE, probs=0.95)
  sites$Ta.50[stat.k]= quantile(dat, na.rm = TRUE, probs=0.50)
  
  #PLOT DENSITIES BY CONTINENT
  d.Ta= density(dat1) 
  if(var=="TMAX") xlabs= "Maximum daily temperature (?C)"
  if(var=="PRCP") xlabs= "Precipitation (mm)"
  
  if(sites$cont[stat.k]=="asia")  mainlab= "Asia"
  if(sites$cont[stat.k]=="noa")  mainlab= "North America"
  
  if(sites[stat.k,"cont"]!=cont.last) {
    plot(d.Ta, xlim=range(-5, 50), ylim=range(0,0.27), main=mainlab, xlab="", ylab="", col=cols.lat[sites[stat.k,"col"]], cex.main=1) 
    
    site.g= as.numeric(sites[ which(sites[,"cont"]==sites[stat.k,"cont"]),"st.lat"])
    
    if(sites[stat.k,"cont"]=="noa") lat.lab1= c("24?N",rep("", length(site.g)-2) ,"65?N")
    if(sites[stat.k,"cont"]=="asia") lat.lab1= c("11?N",rep("", length(site.g)-2) ,"74?N")
    legend("topleft", fill=rev(cols.lat),border=rev(cols.lat), cex=1,legend=rev(lat.lab1), bty="n", y.intersp=0.5)
    
  }
  if(sites[stat.k,"cont"]==cont.last) points(d.Ta, col=cols.lat[sites[stat.k,"col"]], type="l")
  
  #Plot 95% quantiles
  abline(v=sites$Ta.99[stat.k], lty="dashed", col=cols.lat[sites[stat.k,"col"]], lwd=0.5)
  
  #update continent
  cont.last= sites$cont[stat.k]
  
} #end loop stat

mtext(xlabs, side=1, line=1, outer=TRUE, cex=1.5)
mtext("Density", side=2, line=1, outer=TRUE, cex=1.5)

dev.off()
