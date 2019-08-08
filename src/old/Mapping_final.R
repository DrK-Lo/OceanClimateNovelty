#####################################################################################################################################################
##
## File     : Mapping _final.R
## History  : 2018/10/10  Created by Áki Jarl Láruson (AJL)
##          : 2018/10/31  Modified and updated with code from Lindsay Veazey (LV) and AJL
##
#####################################################################################################################################################
##
## This script subsets and visualizes historic, modern, and future predicted oceonographic data, arranged and provided by Liqing Jiang (LQJ)
## Have SST, Arag, Calc, & pH
#####################################################################################################################################################

setwd("~/Desktop/PostDoc/OceanClimateNovelty/")
source("src/Novelty_Oceans_Functions.R")

# IDW: http://pro.arcgis.com/en/pro-app/help/analysis/geostatistical-analyst/how-inverse-distance-weighted-interpolation-works.htm
# tmap vignette: https://cran.r-project.org/web/packages/tmap/tmap.pdf
# setwd(../data/)

dat<-fread("data/large_files/Katie_T_Ar_Ca_pH_RCP85.txt", sep = ",")
#dat<-fread("data/large_files/Katie_T_Ar_Ca_pH_RCP45.txt", sep = ",")

#for(i in 1:nrow(dat)){
#  if(dat$Lon[i]>360){
#    dat$Lon[i]<-dat$Lon[i]-360
#  }
#}

#Get 1800's data (1800, 1810, 1820, 1830)
EJ<-dat[dat$Year<1970,]

#Modify to column names
EJ<-EJ[EJ$Month==1,]
EJA<-EJ[,c("Lon","Lat","SST")] #Options are SST, Arag, Calc, and pH
rownames(EJA)<-1:nrow(EJA)

#Get 1900's data (1970, 1980, 1990, 2000)
NJ<-dat[which(dat$Year>1830 & dat$Year<2070),]
NJ<-NJ[NJ$Month==1,]
NJA<-NJ[,c("Lon","Lat","SST")] #Options are SST, Arag, Calc, and pH
rownames(NJA)<-1:nrow(NJA)

#Get 2000's data (2070, 2080, 2090, 2100)
TJ<-dat[dat$Year>2000,]
TJ<-TJ[TJ$Month==1,]
TJA<-TJ[,c("Lon","Lat","SST")] #Options are SST, Arag, Calc, and pH
rownames(TJA)<-1:nrow(TJA)

EJA[Lon>360,Lon:=Lon-360]
NJA[Lon>360,Lon:=Lon-360]
TJA[Lon>360,Lon:=Lon-360]

#EJA <- SpatialPoints(EJA , proj4string=CRS("+proj=longlat +datum=WGS84")) # this is your spatial points data frame

#head(EJA)
EJA$Lon <- EJA$Lon - 180 # convert to WGS 1984 bounds
EJA <- SpatialPoints(EJA) # this is your spatial points df

NJA$Lon <- NJA$Lon - 180 # convert to WGS 1984 bounds
NJA <- SpatialPoints(NJA) # this is your spatial points df

TJA$Lon <- TJA$Lon - 180 # convert to WGS 1984 bounds
TJA <- SpatialPoints(TJA) # this is your spatial points df

# Project sp object to WGS 84
proj4string(EJA) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
proj4string(NJA) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
proj4string(TJA) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

# Create an empty grid; n = number of cells
# Increase n to increase resolution

#gr <- as.data.frame(spsample(EJA, 'regular', n  = 50000))
#gr <- as.data.frame(spsample(NJA, 'regular', n  = 50000))
#gr <- as.data.frame(spsample(TJA, 'regular', n  = 50000))

r<-NULL

MakeRast<-function(data){
  while(class(r)=="NULL"){
  gr <- as.data.frame(spsample(data, 'regular', n  = 50000))
  names(gr) <- c('X', 'Y')
  coordinates(gr) <- c('X', 'Y')
  gridded(gr) <- TRUE  
  fullgrid(gr) <- TRUE 
  proj4string(gr) <- proj4string(data)# If this line throws an error, run this
  # chunk again. It is a random grid.
  data.idw <- idw(pH ~ 1, data, newdata = gr, idp = 2)
  r <- raster(data.idw)
  }
}
  
MakeRast(EJA)

names(gr) <- c('X', 'Y')
coordinates(gr) <- c('X', 'Y')
gridded(gr) <- TRUE  
fullgrid(gr) <- TRUE  # Create SpatialGrid object
#proj4string(gr) <- proj4string(EJA) # If this line throws an error, run this
# chunk again. It is a random grid.
#proj4string(gr) <- proj4string(NJA)
proj4string(gr) <- proj4string(TJA)

# Interpolate the grid cells using power value = 2
EJA.idw <- idw(pH ~ 1, EJA, newdata = gr, idp = 2)

# Convert to raster
r <- raster(EJA.idw)

#require(doParallel)
#cores <- 3
#cl <- makeCluster(cores)
#registerDoParallel(cl)

# Interpolate the grid cells using power value = 2
# Arag ~ 1 = simple kriging

EJA.idw <- idw(pH ~ 1, EJA, newdata = gr, idp = 2)

NJA.idw <- idw(pH ~ 1, NJA, newdata = gr, idp = 2)

TJA.idw <- idw(pH ~ 1, TJA, newdata = gr, idp = 2)

#stopCluster()

# Convert to raster
r <- raster(EJA.idw)
r <- raster(NJA.idw)
r <- raster(TJA.idw)

# Load continents
data('World', package = 'tmap') 
# World is a {sf} df of countries

# Plot
tm_shape(r) + 
    tm_raster(n = 12, palette = 'Blues', # n = 10 may be better
            #title = "pH \n in Jan. 1800") + 
            #title = "pH \n in Jan. 1970") + 
            title = "pH \n in Jan. 2070") + 
  tm_legend(legend.outside = TRUE) +
  tm_shape(World) + 
      tm_fill()
#
