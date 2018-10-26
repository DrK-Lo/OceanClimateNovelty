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
##################
#### Read in the input data 
##################
##these files were created in a separate script "Novelty_NA_LocalStPCA_InputData_Feb2016.R", except for the cru surrogates, which were created in the Oct2015 version. 

setwd("/Users/katie/Desktop/OceanClimateNovelty/data/") 
  ## specify location of data

dat <- fread("large_files/Data_OceNov.txt", sep = ",")
head(dat)
unique(dat$Year)

ggplot(dat, aes( Lon, Lat))+
  geom_hex(binwidth = c(5, 5))

dat_depth <- fread("Data_subsurface.txt", sep=",")
head(dat_depth)
sort(unique(dat_depth$Year))

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

# 40-year climate normals
  # summer  
    # (if Lat > 0, months 6,7,8)
    # (if Lat < 0, months 12,1,2)
  # winter
    # (if Lat < 0, months 6,7,8)
    # (if Lat > 0, months 12,1,2)
  # means: SST_summer, SST_winter, Arag_summer, Arag_winter, 
          # Calc_summer, Calc_winter
  
  ## Summer calculations
  x_sum <- dat %>% filter((Lat <0 & Month %in% c(12,1,2))|
                        (Lat >0 & Month %in% c(6,7,8)))
  SST_sum <- tapply(x_sum$SST,INDEX = x_sum$No,mean, rm.na=TRUE)
    length(SST_sum)
  Arag_sum <- tapply(x_sum$Arag,INDEX = x_sum$No,mean, rm.na=TRUE)
    length(Arag_sum)
  Calc_sum <- tapply(x_sum$Calc,INDEX = x_sum$No,mean, rm.na=TRUE)
    length(Calc_sum)
  identical(names(SST_sum), names(Arag_sum))
  identical(names(Arag_sum), names(Calc_sum))
  smr <- data.frame(No=names(SST_sum), SST_sum,Arag_sum,Calc_sum  )
  head(smr)
  
  ## Winter calculations    
  x_win <- dat %>% filter((Lat >0 & Month %in% c(12,1,2))|
                            (Lat <0 & Month %in% c(6,7,8)))
  SST_win <- tapply(x_win$SST,INDEX = x_win$No,mean, rm.na=TRUE)
  Arag_win <- tapply(x_win$Arag,INDEX = x_win$No,mean, rm.na=TRUE)
  Calc_win <- tapply(x_win$Calc,INDEX = x_win$No,mean, rm.na=TRUE)
  identical(names(SST_win), names(Arag_win))
  identical(names(Arag_win), names(Calc_win))
  wnt <- data.frame(No=names(SST_win), SST_win,Arag_win,Calc_win  )
  head(wnt)
  
## Climate data. note that all precipitation variables have been log-transformed
A <- read.csv("X.NAnaec8.ref.csv") 
  # 1971-2000 climate normals for all land cells of the DEM 
  head(A)
  dim(A)
  # max and min temperatures for each season; precip for each season

B <- read.csv("X.NAnaec8.proj_GlobalMean_RCP45_2085.csv") 
  # 2071-2100 climate normals (RCP4.5 ensemble mean projection) for all land cells of the DEM. 
  head(B)
  dim(B)
  # same columns as A
  
C <- read.csv("X.stn_detrended.csv") 
  # ICV proxy data. Linearly detrended 1951-1990 annual time series at selected 
  # CRU TS3.23 climate stations. these time series are used as proxies for 
  # local interannual climate variability ("ICV proxies")
  # same columns as A and B
  head(C, 50)
  dim(C)


C.id <- read.csv("A.stn_detrended.csv")[,1] 
  # the id number for each ICV proxy. 
  length(C.id) # same rows as C
  nlevels(as.factor(C.id))

proxy <- read.csv("grid4nn_NAnaec8.csv")[,1]  
  # the ICV proxy used for each grid cell. 
  # need to read paper again for this one too - 
  head(proxy)
  length(proxy) # same rows as A and B
  nlevels(as.factor(proxy)) # it appears some C.id's not used?

## subsample the analog pool to reduce processing time
subsample <- read.csv("subsample.NAnaec8.csv")[,1]  
  # need to read paper again for this one, assume kept one of multiple sites with same analog
A <- A[subsample,] 
dim(A)

########################
### Calculation of sigma dissimilarity
########################

# Principal component truncation rule
trunc.SDs <- 0.1 #truncation 

#initiate the data frame to store the projected sigma dissimilarity of best analogs for each grid cell. 
NN.sigma <- rep(NA,length(proxy))

for(j in sort(unique(proxy))){       
  # run the novelty calculation once for each ICV proxy. 
  # Takes about 1.5 sec/iteration (1 hour total) on a typical laptop. 
  
  ## Select data relevant to ICV proxy j
  Bj <- B[which(proxy==j),]   # future climates
    # select locations for which ICV proxy j is the closest ICV proxy. 
  Cj <- C[which(C.id==j),]    # reference period (1951-1990) ICV at ICV proxy j
  
  ## Step 1: express climate data as standardized anomalies of reference period 
  # (1951-1990) ICV at ICV proxy j. 
  Cj.sd <- apply(Cj,2,sd, na.rm=T)  
    #standard deviation of 1951-1990 interannual variability in each climate 
    # variable, ignoring missing years
  A.prime <- sweep(A,MARGIN=2,Cj.sd,`/`) #standardize the reference ICV
    # a <- matrix(c(1,2,3,4,5,6), nrow=2)
    # sweep(a, MARGIN =2, STATS=c(2,3,4)) # subtracts from each column
    # sweep(a, MARGIN =2, STATS=c(2,3,4), FUN=`/`) # divides each column by STATS
  Bj.prime <- sweep(Bj,MARGIN=2,Cj.sd,`/`) #standardize the analog pool    
  Cj.prime <- sweep(Cj,MARGIN=2,Cj.sd,`/`) #standardize the projected future conditions of grid cells represented by ICV proxy j
  
  ## Step 2: Extract the principal components (PCs) of the reference period ICV 
    # and project all data onto these PCs
  PCA <- prcomp(Cj.prime[!is.na(apply(Cj.prime,1,mean)),])   
    # Principal components analysis. The !is.na(apply(...)) term is there 
    # simply to select all years with complete observations in all variables. 
  PCs <- max(which(unlist(summary(PCA)[1])>trunc.SDs))    
    # find the number of PCs to retain using the PC truncation 
    # rule of eigenvector stdev > the truncation threshold
  X <- as.data.frame(predict(PCA,A.prime))   
    # project the analog pool onto the PCs
  Yj <- as.data.frame(predict(PCA,Bj.prime)) 
    # project the projected future conditions onto the PCs
  Zj <- as.data.frame(predict(PCA,Cj.prime)) 
    # project the reference ICV onto the PCs
  
  ## Step 3a: express PC scores as standardized anomalies of reference interannual variability 
  Zj.sd <- apply(Zj,2,sd, na.rm=T)     
    #standard deviation of 1951-1990 interannual variability in each principal component, ignoring missing years
  X.prime <- sweep(X,MARGIN=2,Zj.sd,`/`) 
    #standardize the analog pool    
  Yj.prime <- sweep(Yj,MARGIN=2,Zj.sd,`/`) 
    #standardize the projected conditions   
  
  ## Step 3b: find the sigma dissimilarity of each projected condition with 
  # its best analog (Euclidean nearest neighbour) in the observed analog pool. 
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
  
  print(j) 
}

write.csv(NN.sigma,"NN.sigma.RCP45.GlobalMean.2085.csv", row.names=FALSE)

########################
### Map of sigma dissimilarity
########################

## import data
dem <- raster("dem8.tif") # the digital elevation model (DEM) used to generate the input data. This is used as a template raster for mapping the results. 
land <- read.csv("land.NAnaec8.csv")[,1] # DEM cells that make up the B matrix. this list includes a coastal buffer that is excluded from the analog pool (A matrix).  
# NN.sigma <- read.csv("NN.sigma.RCP45.GlobalMean.2085.csv")[,1]  # read in if not already in memory

## map (exported to working directory via the png() and dev.off() calls)
png(filename=paste("NoveltyMap.png",sep="_"), type="cairo", units="in", width=9, height=8, pointsize=16, res=800)
par(mar=c(0,0,0,0))
xl <- 2550000; yb <- -3000000; xr <- 2850000; yt <- -200000
breakseq <- c(0,2,4)
breakpoints <- c(seq(breakseq[1], breakseq[3], 0.01),9); length(breakpoints)
ColScheme <- colorRampPalette(c("light gray", "black"))(length(breakpoints)-1)  #color scheme used in the manuscript
# ColScheme <- colorRampPalette(c("#0000CD", blue2red(length(breakpoints)-1), "#CD0000"))(length(breakpoints)-1) #alternate color scheme 
X <- dem
values(X) <- NA
values(X)[land] <- NN.sigma
plot(X, xaxt="n", yaxt="n", col=ColScheme, breaks=breakpoints, legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 
rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme, border=NA)
rect(xl,  yb,  xr,  yt,  col=NA)
text(rep(xr,3),c(yb,mean(c(yt,yb)),yt-100000),sapply(c(bquote(.(breakseq[1])*sigma), bquote(.(breakseq[2])*sigma),bquote(.(breakseq[3])*sigma)),as.expression),pos=4,cex=1,font=1)  
text(xl-180000,mean(c(yb,yt)), "Dissimilarity of best analog", font=2, cex=1, srt=90) 
rect(xl,  yt+100000,  xr,  yt+300000,  col=ColScheme[length(ColScheme)])
text(xr,  yt+200000,  bquote(">"*.(breakseq[3])*sigma),pos=4,cex=1,font=1)  
box(col="white", lwd=1.5)
dev.off()


#####################
####### Code for preparing ClimateNA files for use as an alternate B matrix for other time periods, other model projections, other climate variables etc.
#####################

# ClimateNA can be obtained here: http://cfcg.forestry.ubc.ca/projects/climate-data/climatebcwna/#ClimateNA
# use the NAnaec8.csv file provided in the working directory to download projected normals for North America ("NA") in a North American equidistant conic ("naec") projection at 8-km resolution. 

# Specify the ClimateNA file attributes
models <-  c("ACCESS1-0","CanESM2","CCSM4","CESM1-CAM5","CNRM-CM5","CSIRO-Mk3-6-0", "GFDL-CM3","GISS-E2R", "GlobalMean","HadGEM2-ES", "INM-CM4", "IPSL-CM5A-MR", "MIROC-ESM", "MIROC5", "MPI-ESM-LR","MRI-CGCM3")    
grid <- "NAnaec8"
model <- models[9]
RCP <- "RCP45"
proj.year <- 2085 #central year of the normal period, e.g., 2085 indicates the 2071-2100 normal period. 
VarCode <- "S" # type of variables selected from ClimateNA. Y is annual, S is seasonal, M is monthly. 

#read in ClimateNA output file
file <- paste(grid,"_",model,"_",RCP,"_",proj.year,VarCode,".csv",sep="")
grid.proj <- read.csv(file, strip.white = TRUE, na.strings = c("NA","",-9999) )
nonCNA <- which(is.na(grid.ref[,6]))  # dem cells outside climateWNA extent

#select predictor variables
predictors <- names(grid.proj)[grep("Tmin|Tmax|PPT",names(grid.proj))]
X.grid.proj <- grid.proj[-nonCNA,which(names(grid.proj)%in%predictors)]  #removes NA cells and selects analysis variables. 

##log-transform precipitation
for(i in grep("PPT",names(X.grid.proj))){X.grid.proj[which(X.grid.proj[,i]==0),i] <- 1}  #set zero values to one, to facilitate log-transformation
X.grid.proj[,grep("PPT",names(X.grid.proj))] <- log(X.grid.proj[,grep("PPT",names(X.grid.proj))]) #log-transform the precipitation variables 

write.csv(X.grid.proj,paste("X.",grid,".proj_",model,"_",RCP,"_",proj.year,".csv",sep=""), row.names=FALSE)  #This file can be used as the B matrix to calculate novelty for other time periods and/or other model projections. 



