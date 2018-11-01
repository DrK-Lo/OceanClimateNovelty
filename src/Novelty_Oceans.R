### KE Lotterhos
### Oct 2018
### Northeastern University

## Based on basic code for calculating and mapping climatic novelty using the sigma dissimilarity metric. 
## written by Colin Mahony
## Supp Material of "A closer look at novel climates: new methods and insights at continental to landscape scales"

library(raster)
library(FNN)
library(RColorBrewer)
library(colorRamps)
library(adehabitatLT)   #provides the chi distribution
library(data.table)
library(tidyverse)
library(ggplot2)
library(fields)
if("ggplot2"%in%installed.packages()){
  require(ggplot2)}
##-----------------------------
#### Read in the input data ####
# -----------------------------
##these files were created in a separate script "Novelty_NA_LocalStPCA_InputData_Feb2016.R", except for the cru surrogates, which were created in the Oct2015 version. 

setwd("/Users/katie/Desktop/OceanClimateNovelty/data/") 
  ## specify location of data

dat <- fread("large_files/Data_OceNov.txt", sep = ",")
head(dat)
unique(dat$Year)
cond <- dat$Lat > 40 & dat$Lat < 50
d2 <- dat[cond,]

ggplot(d2, aes( Month, SST))+
  geom_hex(binwidth = c(0.5, 0.5))

ggplot(dat, aes( Lon, Lat))+
  geom_hex(binwidth = c(5, 5))

dat_depth <- fread("Data_subsurface.txt", sep=",")
head(dat_depth,100)
sort(unique(dat_depth$Year))
dat_depth$MonthDay <- dat_depth$Month + dat_depth$Day/30.5

sort(unique(dat_depth$Depth))

#cond <- dat_depth$No==5000 & dat_depth$Depth < 30
#plot(dat_depth$MonthDay[cond], dat_depth$Temp[cond])

ggplot(dat_depth, aes( Lon, Lat))+
  geom_hex(binwidth = c(5, 5))

# Plans
# Compare 1800 analog (A) to today surface (B), with today's surface variability (C)
# Compare today surface analog (A) to 2100 (B), with today's surface variability (C)
# Compare today all depths analog (A) to 2100 (B), with today's surface variability (C)


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

calculate_normals <- function(dat1){
  # input the data frame for the span of years you want the normals for
  month1 <- c(1,2,3, 4)
  month2 <- c(7, 8,9,10)
  
  ## Summer calculations
  x_sum <- dat1 %>% filter((Lat <0 & Month %in% month1)|
                        (Lat >0 & Month %in% month2))
  SST_sum <- tapply(x_sum$SST,INDEX = x_sum$No,mean, rm.na=TRUE)
    length(SST_sum)
  Arag_sum <- tapply(x_sum$Arag,INDEX = x_sum$No,mean, rm.na=TRUE)
    length(Arag_sum)
  Calc_sum <- tapply(x_sum$Calc,INDEX = x_sum$No,mean, rm.na=TRUE)
    length(Calc_sum)
  identical(names(SST_sum), names(Arag_sum))
  identical(names(Arag_sum), names(Calc_sum))
  smr <- data.frame(No=as.integer(names(SST_sum)), SST_sum,Arag_sum,Calc_sum  )
  head(smr)
  
  ## Winter calculations    
  x_win <- dat1 %>% filter((Lat >0 & Month %in% month1)|
                            (Lat <0 & Month %in% month2))
  SST_win <- tapply(x_win$SST,INDEX = x_win$No,mean, rm.na=TRUE)
  Arag_win <- tapply(x_win$Arag,INDEX = x_win$No,mean, rm.na=TRUE)
  Calc_win <- tapply(x_win$Calc,INDEX = x_win$No,mean, rm.na=TRUE)
  identical(names(SST_win), names(Arag_win))
  identical(names(Arag_win), names(Calc_win))
  wnt <- data.frame(No=as.integer(names(SST_win)), SST_win,Arag_win,Calc_win  )
  head(wnt)
  
  # merge summer and winter data frames
  identical(names(SST_win), names(SST_sum))
  normals <- full_join(smr, wnt, by="No")
  head(normals)
  nrow(normals) == nrow(wnt)
  nrow(normals) == nrow(smr)
  return(normals[order(normals$No),])
}

norm_1800 <- calculate_normals(dat_1800)
norm_2000 <- calculate_normals(dat_2000)
norm_2100 <- calculate_normals(dat_2100)
head(norm_1800)
head(norm_2000)

#--------------------------------  
### Detrend for ICV and get sd ####
### Note I don't think this is a conservative way to do this
#--------------------------------

get_detrend <- function(dat1){
  
  dat1$Month2 <- dat1$Month^2
  dat1$Month3 <- dat1$Month^3
  
  get_sd_resids <- function(i){
    # dat1 <- dat_2000
    # plot(dat1$SST[dat1$No==i]~dat1$Month[dat1$No==i])
    SST_sd <- sd(residuals(lm(SST~Month + Month2 + Month3,
              data=dat1[dat1$No==i,])))
    Arag_sd <- sd(residuals(lm(Arag~Month + Month2 + Month3,
                            data=dat1[dat1$No==i,])))
    Calc_sd <- sd(residuals(lm(Calc~Month + Month2 + Month3,
                            data=dat1[dat1$No==i,])))
    #points(dat1$Month[dat1$No==i],predict(mod), pch=19)
    return(c(SST_sd, Arag_sd, Calc_sd))
  }
  
  stations <- sort(unique(dat1$No))
    sds<- sapply(stations, get_sd_resids)
      # this takes ~10 minutes to run
  out <- data.frame(stations, t(sds), t(sds))
  colnames(out) <- c("No", "SST_sum", "Arag_sum", "Calc_sum",
                     "SST_win", "Arag_win", "Calc_win")
}


#--------------------------------  
### Get annual ICV
#--------------------------------

get_ICV <- function(dat1){
  get_sd_resids <- function(i){
    # dat1 <- dat_2000
    # plot(dat1$SST[dat1$No==i]~dat1$Month[dat1$No==i])
    SST_sd <- sd(dat1$SST[dat1$No==i])
    Arag_sd <- sd(dat1$Arag[dat1$No==i])
    Calc_sd <- sd(dat1$Calc[dat1$No==i])
    #points(dat1$Month[dat1$No==i],predict(mod), pch=19)
    return(c(SST_sd, Arag_sd, Calc_sd))
  }
  stations <- sort(unique(dat1$No))
  sds<- sapply(stations, get_sd_resids)
  # this takes ~10 minutes to run
  out <- data.frame(stations, t(sds), t(sds))
  colnames(out) <- c("No", "SST_sum", "Arag_sum", "Calc_sum",
                     "SST_win", "Arag_win", "Calc_win")
  return(out)
}

sd_stations <- get_ICV(dat_2000)
  # this takes about 10 min to run
head(sd_stations) 

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
### function for plotting ####
#--------------------------------
Plot_nonInt<-function(lat, long, var, refMap, legend_name){
  sampData <- data.frame(lat, long, var)

    ggplot()+
      geom_polygon(data = refMap, aes(x=long, y = lat, group=group))+
      stat_summary_2d(data=sampData, aes(x=long, y = lat, z= var), bins=80, alpha = 0.8)+
      #theme(legend.position=c(50,100))+
      #scale_fill_brewer(palette="Dark2") +
      scale_fill_gradientn(name=legend_name, 
                           colours=two.colors(40,start = "blue", 
                                              end="red", middle="orange")) +
      coord_fixed() 
}

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

B

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
which(is.na(NN.sigma))
tail(sort(NN.sigma))
length(NN.sigma)
world <- map_data("world2")
dim(stationInfo)
B2 <- data.frame(No=B$No, NN.sigma)
B2 <- merge(B2, stationInfo, by.x="No", by.y="stations", all.x=TRUE)
Plot_nonInt(B2$lat, B2$long, 
            B2$NN.sigma, world, "sigma dis.")
#write.csv(NN.sigma,"NN.sigma.RCP45.GlobalMean.2085.csv", row.names=FALSE)

#-----------------------------
### 