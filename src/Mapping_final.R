#####################################################################################################################################################
##
## File     : Mapping _final.R
## History  : 2018/10/10  Created by Áki Jarl Láruson (AJL)
##          : 2018/10/31  Modified and updated by Lindsay Veazey (LV), AJL added additional data subsetting commands for visualization, file uploaded to GitHub
##
#####################################################################################################################################################
##
## This script subsets and visualizes historic, modern, and future predicted oceonographic data, arranged and provided by Liqing Jiang (LQJ)
## 
#####################################################################################################################################################

# IDW: http://pro.arcgis.com/en/pro-app/help/analysis/geostatistical-analyst/how-inverse-distance-weighted-interpolation-works.htm
# tmap vignette: https://cran.r-project.org/web/packages/tmap/tmap.pdf
# setwd(../data/)
dat<-read.csv("Data_Lotterhos.csv")

#for(i in 1:nrow(dat)){
#  if(dat$Lon[i]>360){
#    dat$Lon[i]<-dat$Lon[i]-360
#  }
#}

#Get 1800's data (1800, 1810, 1820, 1830)
EJ<-dat[dat$Year<1970,]

EJ<-EJ[EJ$Month==1,]
EJA<-EJ[,c(2,3,7)] #6 for SST, 7 for Arag, 8 for Calc
rownames(EJA)<-1:nrow(EJA)

#Get 1900's data (1970, 1980, 1990, 2000)
NJ<-dat[which(dat$Year>1830 & dat$Year<2070),]
NJ<-NJ[NJ$Month==1,]
NJA<-NJ[,c(2,3,7)]
rownames(NJA)<-1:nrow(NJA)

#Get 2000's data (2070, 2080, 2090, 2100)
TJ<-dat[dat$Year>2000,]
TJ<-TJ[TJ$Month==1,]
TJA<-TJ[,c(2,3,7)]
rownames(TJA)<-1:nrow(TJA)

for(i in 1:nrow(EJA)){
  if(EJA$Lon[i]>360){
    EJA$Lon[i]<-EJA$Lon[i]-360
  }
}

#EJA <- SpatialPoints(EJA , proj4string=CRS("+proj=longlat +datum=WGS84")) # this is your spatial points data frame


library(raster)
library(rgdal)
library(tmap)
library(gstat) 
library(sp)   

head(EJA)
EJA$Lon <- EJA$Lon - 180 # convert to WGS 1984 bounds
EJA <- SpatialPoints(EJA) # this is your spatial points df

# Project sp object to WGS 84
proj4string(EJA) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

# Create an empty grid; n = number of cells
# Increase n to increase resolution
gr <- as.data.frame(spsample(EJA, 'regular', n  = 50000))
names(gr) <- c('X', 'Y')
coordinates(gr) <- c('X', 'Y')
gridded(gr) <- TRUE  
fullgrid(gr) <- TRUE  # Create SpatialGrid object
proj4string(gr) <- proj4string(EJA) # If this line throws an error, run this
# chunk again. It is a random grid.

# Interpolate the grid cells using power value = 2
# Arag ~ 1 = simple kriging
EJA.idw <- idw(Arag ~ 1, EJA, newdata = gr, idp = 2)

# Convert to raster
r <- raster(EJA.idw)

# Load continents
data('World', package = 'tmap') 
# World is a {sf} df of countries

# Plot
# Per the tmap vignette, I think I can layer shapes..
    tm_shape(r) + 
    tm_raster(n = 12, palette = 'Blues', # n = 10 may be better
            title = "Aragonite concentration \n in Jan. 1800 (units here)") + 
  tm_legend(legend.outside = TRUE) +
  tm_shape(World) + 
      tm_fill()
#
