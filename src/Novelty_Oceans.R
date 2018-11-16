### KE Lotterhos
### Oct 2018
### Northeastern University
### Mod. Áki - Nov. 2018

## Based on basic code for calculating and mapping climatic novelty using the sigma dissimilarity metric. 
## written by Colin Mahony
## Supp Material of "A closer look at novel climates: new methods and insights at continental to landscape scales"

# Plans
# Compare 1800 analog (A) to today surface (B), with today's surface variability (C)
# Compare today surface analog (A) to 2100 (B), with today's surface variability (C)
# Compare today all depths analog (A) to 2100 (B), with today's surface variability (C)

## Specify location of data
#setwd("/Users/katie/Desktop/OceanClimateNovelty/") 
setwd("~/Desktop/PostDoc/NovelOceanClim/OceanClimateNovelty/")

#Create function that removes previous user installed packages to avoid masking
clean_pkgs<-function(){
  lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE)
}
#clean_pkgs() #Remove all non-essential previously called packages

source("src/Novelty_Oceans_Functions.R")

##-----------------------------
#### Read in the input data ####
# -----------------------------
##these files were created in a separate script "Novelty_NA_LocalStPCA_InputData_Feb2016.R", except for the cru surrogates, which were created in the Oct2015 version. 

<<<<<<< HEAD
dat <- fread("data/large_files/Katie_T_Ar_Ca_pH_RCP85.txt", sep = ",")
=======
#dat <- fread("data/large_files/Data_OceNov.txt", sep = ",")
dat <- fread("data/large_files/ESM2M_2000_RCP8.5.txt", sep = ",")
>>>>>>> 2f7bb9b4ddcbe7004e37d5725fc5a7a2fb58dc1b
head(dat)
unique(dat$Year)
cond <- dat$Lat > 40 & dat$Lat < 50
d2 <- dat[cond,]

#ggplot(d2, aes( Month, SST))+
#  geom_hex(binwidth = c(0.5, 0.5))

#ggplot(dat, aes( Lon, Lat))+
#  geom_hex(binwidth = c(5, 5))

### Load depth data
  #dat_depth <- fread("data/Data_subsurface.txt", sep=",")
  #head(dat_depth,100)
  #hist(dat_depth$Depth)
  #sort(unique(dat_depth$Year))
  #dat_depth$MonthDay <- dat_depth$Month + dat_depth$Day/30.5
  #sort(unique(dat_depth$Depth))

#cond <- dat_depth$No==5000 & dat_depth$Depth < 30
#plot(dat_depth$MonthDay[cond], dat_depth$Temp[cond])

#ggplot(dat_depth, aes( Lon, Lat))+
#  geom_hex(binwidth = c(5, 5))


dat_1800 <- dat %>% filter(Year<1850)
dim(dat_1800)
dat_2000 <-   dat %>% filter(Year>1960 & Year<2010)
dim(dat_2000)  
dat_2100 <-   dat %>% filter(Year>2050)
dim(dat_2100) 

#--------------------------------
#### 40-year climate normals ####
#--------------------------------
  # summer  
    # (if Lat > 0, months 6,7,8)
    # (if Lat < 0, months 12,1,2)
  # winter
    # (if Lat < 0, months 6,7,8)
    # (if Lat > 0, months 12,1,2)
  # means: SST_summer, SST_winter, Arag_summer, Arag_winter, 
          # Calc_summer, Calc_winter

norm_1800 <- calculate_normals(dat_1800)
norm_2000 <- calculate_normals(dat_2000)
norm_2100 <- calculate_normals(dat_2100)
head(norm_1800)
head(norm_2000)
head(norm_2100)

#--------------------------------  
### data frame to link station number to Lat Long ####
#--------------------------------
head(dat)
min(which(dat$No==2))
stations <- unique(dat$No)
stationInfo = data.frame(stations=stations, lat=NA, long=NA)
unik <- which(!duplicated(dat$No))
head(unik)
stationInfo$lat <- dat$Lat[unik]
stationInfo$long <- dat$Lon[unik]

head(stationInfo)
which(!complete.cases(stationInfo))




#--------------------------------  
### 1800 analog to today ####
#--------------------------------
length(norm_1800$No)
length(norm_2000$No)
  # subset the 2000 data to the same stations as the 1800 data

A <- norm_1800
  # 1800-1830 climate normals
  head(A)
  dim(A)
  
B <- norm_2000[norm_2000$No %in% norm_1800$No,]
  # 1970-2000 climate normals 
  head(B)
  dim(B)
  # same columns as A
  
  # sanity check to make sure stations in right order
  identical(A$No, B$No) # should be true
  
  head(dat_2000)  
C <- data.frame(dat_2000[,c(1,6,7,8)], dat_2000[,c(6,7,8)])
head(C)

C.id <- C$No
proxy <- B$No

###-----------------------
### Calculation of sigma dissimilarity of today from 1800 ####
###-----------------------

# Principal component truncation rule
trunc.SDs <- 0.1 #truncation 
  
#initiate the data frame to store the projected sigma dissimilarity of best analogs for each grid cell. 
NN.sigma <- rep(NA,length(proxy))

for(j in sort(unique(proxy))){       
  # run the novelty calculation once for each ICV proxy. 
  # Takes about 1.5 sec/iteration on a typical laptop. 
  
  ## Select data relevant to ICV proxy j
  Bj <- B[which(proxy==j),]   # future climates
    # select locations for which ICV proxy j is the closest ICV proxy. 
  Cj <- C[which(C.id==j),]    # reference period ICV at ICV proxy j
  
  ## Step 1: express climate data as standardized anomalies of reference period 
  #  ICV at ICV proxy j. 
  Cj.sd <- apply(Cj,2,sd, na.rm=T)  #standard deviation of interannual variability in each climate variable, ignoring missing years
    #standard deviation of variability in each climate 
    # variable, ignoring missing years
  A.prime <- sweep(A[,2:7],MARGIN=2,Cj.sd[2:7],`/`) #standardize the reference ICV
    # a <- matrix(c(1,2,3,4,5,6), nrow=2)
    # sweep(a, MARGIN =2, STATS=c(2,3,4)) # subtracts STATs from each column
    # sweep(a, MARGIN =2, STATS=c(2,3,4), FUN=`/`) # divides each column by STATS
  Bj.prime <- sweep(Bj[,2:7],MARGIN=2,Cj.sd[2:7],`/`) #standardize the analog pool    
  Cj.prime <- sweep(Cj[,2:7],MARGIN=2,Cj.sd[2:7],`/`) #standardize the projected future conditions of grid cells represented by ICV proxy j
  
  colnames(Cj.prime) <- colnames(A.prime)
  ## Step 2: Extract the principal components (PCs) of the reference period ICV 
    # and project all data onto these PCs
  PCA <- prcomp(Cj.prime[!is.na(apply(Cj.prime,1,mean)),])   
    # Principal components analysis. The !is.na(apply(...)) term is there 
    # simply to select all years with complete observations in all variables. 
  PCA$rotation
  
  #plot(PCA$rotation[,1], PCA$rotation[,2], xlim=c(-0.43, -0.39))
  #text(PCA$rotation[,1], PCA$rotation[,2], rownames(PCA$rotation))
  # SST summer and winter in lower right of PC space, 
  # Arag and Calc in upper left of PC space
  
  #plot(PCA$rotation[,1], PCA$rotation[,3], xlim=c(-0.43, -0.39))
  #text(PCA$rotation[,1], PCA$rotation[,3], rownames(PCA$rotation))
  # separates the three variables
  
  #plot(PCA$sdev)
  #round(PCA$sdev, 2)
  
  PCs <- max(which(unlist(summary(PCA)[1])>trunc.SDs))    
    # find the number of PCs to retain using the PC truncation 
    # rule of eigenvector stdev > the truncation threshold
  X <- as.data.frame(predict(PCA,A.prime))   
    # project the analog pool onto the PCs
    head(X)
    
  Yj <- as.data.frame(predict(PCA,Bj.prime)) 
    # project the projected future conditions onto the PCs
    
  Zj <- as.data.frame(predict(PCA,Cj.prime)) 
    # project the reference ICV onto the PCs
    
  
    #plot(X[,1], X[,2]) # analog
    #points(Zj[,1], Zj[,2], pch=19, col=rgb(1,0,0, 0.5)) # reference ICV
    #points(Yj[,1], Yj[,2], pch=19, col=rgb(0,1,0)) # future conditions
    
        
  ## Step 3a: express PC scores as standardized anomalies of reference interannual variability 
  Zj.sd <- apply(Zj,2,sd, na.rm=T)     
    #standard deviation of 1951-1990 interannual variability in each principal component, ignoring missing years
    #Zj.sd
  X.prime <- sweep(X,MARGIN=2,Zj.sd,`/`) 
    #standardize the analog pool  
    #head(X.prime)
  Yj.prime <- sweep(Yj,MARGIN=2,Zj.sd,`/`) 
    #standardize the projected conditions   
    #Yj.prime
  ## Step 3b: find the sigma dissimilarity of each projected condition with 
  # its best analog (Euclidean nearest neighbour) in the observed analog pool.
    X.prime <- X.prime[complete.cases(X.prime),]
  NN.dist <- as.vector(get.knnx(data=X.prime[,1:PCs],
                                query=Yj.prime[,1:PCs],
                                k=1,algorithm="brute")[[2]]) 
    # Euclidean nearest neighbour distance in the z-standardized PCs of 
    # interannual climatic variability, i.e. the Mahalanobian nearest neighbour. 
  NN.chi <- pchi(NN.dist,PCs) # percentile of the nearest neighbour 
    # distance on the chi distribution with degrees of freedom 
    # equaling the dimensionality of the distance measurement (PCs)
  NN.sigma[which(proxy==j)] <- qchi(NN.chi,1) 
    # values of the chi percentiles on a standard half-normal distribution (chi distribution with one degree of freedom)
  
  if(j%%10==0){print(j)}
}

dim(A)
dim(B)
head(NN.sigma)

NN.sigma[which(is.infinite(NN.sigma))] <- NA
#which(is.na(NN.sigma))
tail(sort(NN.sigma))
length(NN.sigma)

B2 <- data.frame(No=B$No, NN.sigma)
B2 <- merge(B2, stationInfo, by.x="No", by.y="stations", all.x=TRUE)

# Visualize
world <- map_data("world2")
dim(stationInfo)

Plot_nonInt(B2$lat, B2$long, 
            B2$NN.sigma, world, "sigma dis.")
#write.csv(NN.sigma,"NN.sigma.RCP45.GlobalMean.2085.csv", row.names=FALSE)
#write.csv(B2,"Sigma.RCP85.today_1800.csv", row.names=FALSE)
B2<-fread("./data/Sigma.RCP85.today_1800.csv")

##Interpolation for visualization
B2a<-B2[!is.na(B2$NN.sigma),]

B2a<-B2a[,c(4,3,2)]

for(i in 1:nrow(B2a)){
  if(B2a$long[i]>360){
    B2a$long[i]<-B2a$long[i]-360
  }
}

EB2 <- SpatialPoints(B2a) # this is your spatial points df

# Project sp object to WGS 84
proj4string(EB2) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

# Create an empty grid; n = number of cells
# Increase n to increase resolution
gr <- as.data.frame(spsample(EB2, 'regular', n  = 50000))
names(gr) <- c('X', 'Y')
coordinates(gr) <- c('X', 'Y')
gridded(gr) <- TRUE  
fullgrid(gr) <- TRUE  # Create SpatialGrid object
proj4string(gr) <- proj4string(EB2) 
#If this line throws an error, run this
# chunk again. It is a random grid.

# Interpolate the grid cells using inverse distance weighing power = 8
# NN.sigma ~ 1 = simple kriging
EB2.idw <- idw(NN.sigma ~ 1, EB2, newdata = gr, idp = 8)

# Convert to raster
r <- raster(EB2.idw)

world<-map("world2", fill=T,plot=F)
y<-map2SpatialPolygons(world, IDs = sapply(strsplit(world$names, ":"), function(x) x[1]), proj4string=CRS("+proj=longlat +datum=WGS84"))
z<-st_as_sf(y)
wr <- raster(z, res = 0.01)
wrld_r <- fasterize(z, wr)
gplot_wrld_r <- gplot_data(wrld_r)

gplot_r <- gplot_data(r)

#Change scale_fill_gradient value to name of variable
ggplot() +
  geom_tile(data = gplot_r, 
            aes(x = x, y = y, fill = value)) +
  geom_tile(data = dplyr::filter(gplot_wrld_r, !is.na(value)), 
            aes(x = x, y = y), fill = "grey20") +
  xlim(0,360) +
  ylim(-77,90) + 
  xlab("Long") +
  ylab("Lat") +
  ggtitle("Sigma dissimilarity: Today from 1800") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_gradient2(expression(paste(sigma," dis.")),
                       low = 'blue', mid = "yellow", high = 'red',
                       midpoint = 4,
                       limits=c(0,8.1),
                       na.value = NA) +
  coord_quickmap()

#--------------------------------  
### Today analog to 2100 ####
#--------------------------------
identical(norm_2000$No, norm_2100$No)

A <- norm_2000
# 1970-2000 climate normals
head(A)
dim(A)

B <- norm_2100
# 2070-2100
head(B)
dim(B)
# same columns as A

plot(A$SST_sum, A$SST_win)
points(B$SST_sum, B$SST_win, pch=19, col=adjustcolor("blue", 0.1))

plot(A$Arag_sum, A$Arag_win)
points(B$Arag_sum, B$Arag_win, pch=19, col=adjustcolor("blue", 0.5))

plot(A$pH_sum, A$pH_win)
points(B$pH_sum, B$pH_win, pch=19, col=adjustcolor("blue", 0.5))

plot(A$SST_sum, A$Arag_sum)
points(B$SST_sum, B$Arag_sum, pch=19, col=adjustcolor("blue", 0.5))

plot(A$SST_sum, A$pH_sum)
points(B$SST_sum, B$pH_sum, pch=19, col=adjustcolor("blue", 0.5))

plot(A$Arag_sum, A$pH_sum)
points(B$Arag_sum, B$pH_sum, pch=19, col=adjustcolor("blue", 0.5))


# sanity check to make sure stations in right order
identical(A$No, B$No) # should be true

head(dat_2000)  
C <- data.frame(dat_2000[,c(1,6:9)], dat_2000[,c(6:9)])
head(C)

C.id <- C$No
proxy <- B$No



###-----------------------
###### Calculation of sigma dissimilarity of today from 2100 ####
###-----------------------

# Principal component truncation rule
trunc.SDs <- 0.1 #truncation 

#initiate the data frame to store the projected sigma dissimilarity of best analogs for each grid cell. 
NN.sigma <- rep(NA,length(proxy))


for(j in sort(unique(proxy))){       
  # run the novelty calculation once for each ICV proxy. 
  # Takes about 1.5 sec/iteration on a typical laptop. 
  
  ## Select data relevant to ICV proxy j
  Bj <- B[which(proxy==j),]   # future climates
  # select locations for which ICV proxy j is the closest ICV proxy. 
  Cj <- C[which(C.id==j),]    # reference period ICV at ICV proxy j
  
  ## Step 1: express climate data as standardized anomalies of reference period 
  #  ICV at ICV proxy j. 
  Cj.sd <- apply(Cj,2,sd, na.rm=T)  #standard deviation of interannual variability in each climate variable, ignoring missing years
  #standard deviation of variability in each climate 
  # variable, ignoring missing years
  A.prime <- sweep(A[,2:7],MARGIN=2,Cj.sd[2:7],`/`) #standardize the reference ICV
  # a <- matrix(c(1,2,3,4,5,6), nrow=2)
  # sweep(a, MARGIN =2, STATS=c(2,3,4)) # subtracts STATs from each column
  # sweep(a, MARGIN =2, STATS=c(2,3,4), FUN=`/`) # divides each column by STATS
  Bj.prime <- sweep(Bj[,2:7],MARGIN=2,Cj.sd[2:7],`/`) #standardize the analog pool    
  Cj.prime <- sweep(Cj[,2:7],MARGIN=2,Cj.sd[2:7],`/`) #standardize the projected future conditions of grid cells represented by ICV proxy j
  
  colnames(Cj.prime) <- colnames(A.prime)
  ## Step 2: Extract the principal components (PCs) of the reference period ICV 
  # and project all data onto these PCs
  PCA <- prcomp(Cj.prime[!is.na(apply(Cj.prime,1,mean)),])   
  # Principal components analysis. The !is.na(apply(...)) term is there 
  # simply to select all years with complete observations in all variables. 
  PCA$rotation
  
  #plot(PCA$rotation[,1], PCA$rotation[,2], xlim=c(-0.43, -0.39))
  #text(PCA$rotation[,1], PCA$rotation[,2], rownames(PCA$rotation))
  # SST summer and winter in lower right of PC space, 
  # Arag and Calc in upper left of PC space
  
  #plot(PCA$rotation[,1], PCA$rotation[,3], xlim=c(-0.43, -0.39))
  #text(PCA$rotation[,1], PCA$rotation[,3], rownames(PCA$rotation))
  # separates the three variables
  
  #plot(PCA$sdev)
  #round(PCA$sdev, 2)
  
  PCs <- max(which(unlist(summary(PCA)[1])>trunc.SDs))    
  # find the number of PCs to retain using the PC truncation 
  # rule of eigenvector stdev > the truncation threshold
  X <- as.data.frame(predict(PCA,A.prime))   
  # project the analog pool onto the PCs
  head(X)
  
  Yj <- as.data.frame(predict(PCA,Bj.prime)) 
  # project the projected future conditions onto the PCs
  
  Zj <- as.data.frame(predict(PCA,Cj.prime)) 
  # project the reference ICV onto the PCs
  
  
  #plot(X[,1], X[,2]) # analog
  #points(Zj[,1], Zj[,2], pch=19, col=rgb(1,0,0, 0.5)) # reference ICV
  #points(Yj[,1], Yj[,2], pch=19, col=rgb(0,1,0)) # future conditions
  
  
  ## Step 3a: express PC scores as standardized anomalies of reference interannual variability 
  Zj.sd <- apply(Zj,2,sd, na.rm=T)     
  #standard deviation of 1951-1990 interannual variability in each principal component, ignoring missing years
  #Zj.sd
  X.prime <- sweep(X,MARGIN=2,Zj.sd,`/`) 
  #standardize the analog pool  
  #head(X.prime)
  Yj.prime <- sweep(Yj,MARGIN=2,Zj.sd,`/`) 
  #standardize the projected conditions   
  #Yj.prime
  ## Step 3b: find the sigma dissimilarity of each projected condition with 
  # its best analog (Euclidean nearest neighbour) in the observed analog pool.
  X.prime <- X.prime[complete.cases(X.prime),]
  NN.dist <- as.vector(get.knnx(data=X.prime[,1:PCs],
                                query=Yj.prime[,1:PCs],
                                k=1,algorithm="brute")[[2]]) 
  # Euclidean nearest neighbour distance in the z-standardized PCs of 
  # interannual climatic variability, i.e. the Mahalanobian nearest neighbour. 
  NN.chi <- pchi(NN.dist,PCs) # percentile of the nearest neighbour 
  # distance on the chi distribution with degrees of freedom 
  # equaling the dimensionality of the distance measurement (PCs)
  NN.sigma[which(proxy==j)] <- qchi(NN.chi,1) 
  # values of the chi percentiles on a standard half-normal distribution (chi distribution with one degree of freedom)
  
  if(j%%10==0){print(j)}
}

dim(A)
dim(B)
head(NN.sigma)
which(is.infinite(NN.sigma))
NN.sigma[which(is.infinite(NN.sigma))] <- NA
which(is.na(NN.sigma))
tail(sort(NN.sigma))
hist(NN.sigma)
length(NN.sigma)

#Visualize

world <- map_data("world2")
dim(stationInfo)
B2 <- data.frame(No=B$No, NN.sigma)
B2 <- merge(B2, stationInfo, by.x="No", by.y="stations", all.x=TRUE)
Plot_nonInt(B2$lat, B2$long, 
            B2$NN.sigma, world, "sigma dis.")
#write.csv(NN.sigma,"NN.sigma.RCP45.GlobalMean.2085.csv", row.names=FALSE)
<<<<<<< HEAD


### Subset
head(dat_2000)
cond_indop <- function(x){x$long<200 & x$long>150 & 
  abs(x$lat) < 5}
cond_antar <- function(x){x$long < 250 & x$long > 150 & x$lat < -60}

head(norm_2000)
head(stationInfo)

A2 <- merge(norm_2000, stationInfo, by.x="No", by.y="stations", all.x=TRUE)
B2 <- merge(norm_2100, stationInfo, by.x="No", by.y="stations", all.x=TRUE)


## Plot indopacific
whichindop_2000 <- which(cond_indop(A2))
whichindop_2100 <-  which(cond_indop(B2))
plot(A$SST_sum, A$Arag_sum, xlim=c(-5, 35), ylim=c(-0.5,6))
  points(A$SST_sum[whichindop_2000], A$Arag_sum[whichindop_2000], col=adjustcolor("green",0.5))
  points(B$SST_sum[whichindop_2100], B$Arag_sum[whichindop_2100], pch=19, col=adjustcolor("darkgreen", 0.5))
  ## Plot antartica
  whichantar_2000 <- which(cond_antar(A2))
  whichantar_2100 <-  which(cond_antar(B2))
  points(A$SST_sum[whichantar_2000], A$Arag_sum[whichantar_2000], xlim=c(-5, 10), ylim=c(-0.5, 3),  col=adjustcolor("lightblue", 0.5))
  points(B$SST_sum[whichantar_2100], B$Arag_sum[whichantar_2100], pch=19, col=adjustcolor("blue", 0.5))

# where is future antartic going to look like?
fa <- A2[which(A$SST_sum<10 &A$Arag_sum<1),] 
Plot_nonInt(fa$lat, fa$long, 
            1, world, "sigma dis.")
=======
#write.csv(B2,"Sigma.RCP85.today_2100.csv", row.names=FALSE)
B2<-fread("./data/Sigma.RCP85.today_2100.csv")

#Visualize with interpolation
B2a<-B2[!is.na(B2$NN.sigma),]

B2a<-B2a[,c(4,3,2)]

for(i in 1:nrow(B2a)){
  if(B2a$long[i]>360){
    B2a$long[i]<-B2a$long[i]-360
  }
}

EB2 <- SpatialPoints(B2a) # this is your spatial points df

# Project sp object to WGS 84
proj4string(EB2) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

# Function to create a random grid empty grid
gr<-makeGrid(EB2)

# Interpolate the grid cells using power value = 2
# NN.sigma ~ 1 = simple kriging

EB2.idw <- idw(NN.sigma ~ 1, EB2, newdata = gr, idp = 8)

# Convert to raster
r <- raster(EB2.idw)

world<-map("world2", fill=T,plot=F)
y<-map2SpatialPolygons(world, IDs = sapply(strsplit(world$names, ":"), function(x) x[1]), proj4string=CRS("+proj=longlat +datum=WGS84"))
z<-st_as_sf(y)
wr <- raster(z, res = 0.01)
wrld_r <- fasterize(z, wr)
gplot_wrld_r <- gplot_data(wrld_r)

gplot_r <- gplot_data(r)

#Change scale_fill_gradient value to name of variable
ggplot() +
  geom_tile(data = gplot_r, 
            aes(x = x, y = y, fill = value)) +
  geom_tile(data = dplyr::filter(gplot_wrld_r, !is.na(value)), 
            aes(x = x, y = y), fill = "grey20") +
  #xlim(2,358) +
  ylim(-78,90) + 
  xlab("Long") +
  ylab("Lat") +
  ggtitle("Sigma dissimilarity: Today from 2100") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_gradient2(expression(paste(sigma," dis.")),
                      low = 'blue', mid = "yellow", high = 'red',
                      midpoint = 4,
                      limits=c(0,8.1),
                      na.value = NA) +
  coord_quickmap()

>>>>>>> 2f7bb9b4ddcbe7004e37d5725fc5a7a2fb58dc1b
