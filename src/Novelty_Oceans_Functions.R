### Novelty_Oceans Functions
### KE Lotterhos
### Oct 2018 - mod. June 2019
### Northeastern University
### Mod. Áki Nov. 2018


#Create function that removes previous user installed packages to avoid masking
#clean_pkgs<-function(){
#  lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE)
#}
#clean_pkgs() #Remove all non-essential previously called packages

packages_needed <- c("raster", "FNN", "RColorBrewer", "colorRamps", "adehabitatLT",
                     "data.table", "tidyverse", "fields", "ggplot2", "hexbin",
                     "rgdal", "tmap", "gstat", "sp", "maptools", "sf", "fasterize",
                     "fansi", "raster", "tmap", "gstat"
                     )

for (i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
}

#  require(ggplot2)}

for (i in 1:length(packages_needed)){
  library( packages_needed[i], character.only = TRUE)
}

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
# Calc_summer, Calc_winter, pH_summer, pH_winter

calculate_normals <- function(dat1){
  # input the data frame for the span of years you want the normals for
  
  month1 <- c(6,7,8)
  month2 <- c(12,1,2)
  
  ## Summer calculations
  x_sum <- dat1 %>% filter((Lat < 0 & Month %in% month1)|
                             (Lat > 0 & Month %in% month2))
  SST_sum <- tapply(x_sum$SST,INDEX = x_sum$No,mean, rm.na=TRUE)
  #length(SST_sum)
  Arag_sum <- tapply(x_sum$Arag,INDEX = x_sum$No,mean, rm.na=TRUE)
  #length(Arag_sum)
  Calc_sum <- tapply(x_sum$Calc,INDEX = x_sum$No,mean, rm.na=TRUE)
  #length(Calc_sum)
  pH_sum <- tapply(x_sum$pH,INDEX = x_sum$No,mean, rm.na=TRUE)
  
  
  # Do we need to delete this?
    #Long <- aggregate(Lon~No, x_sum, paste, simplify = F) #Get Long
    #Lon <- as.numeric(lapply(Long$Lon, `[[`, 1))
    #Lati <- aggregate(Lat~No, x_sum, paste, simplify = F) #Get Lat
    #Lat <- as.numeric(lapply(Lati$Lat, `[[`, 1))
  
  # check that the column names are identical:
  if(!(identical(names(SST_sum), names(Arag_sum)) |
       identical(names(Arag_sum), names(Calc_sum)) |
       identical(names(pH_sum), names(SST_sum)))
    ){break}
  
  # create summer data frame for each station
  smr <- data.frame(No=as.integer(names(SST_sum)), SST_sum,Arag_sum,Calc_sum, pH_sum  )
  head(smr)
  
  ## Winter calculations    
  x_win <- dat1 %>% filter((Lat >0 & Month %in% month1)|
                             (Lat <0 & Month %in% month2))
  SST_win <- tapply(x_win$SST,INDEX = x_win$No,mean, rm.na=TRUE)
  Arag_win <- tapply(x_win$Arag,INDEX = x_win$No,mean, rm.na=TRUE)
  Calc_win <- tapply(x_win$Calc,INDEX = x_win$No,mean, rm.na=TRUE)
  pH_win <- tapply(x_win$pH,INDEX = x_win$No,mean, rm.na=TRUE)

  wnt <- data.frame(No=as.integer(names(SST_win)), SST_win,Arag_win,Calc_win,pH_win  )
  head(wnt)
  
  # Do we need to delete this?
  #Long <- aggregate(Lon~No, x_win, paste, simplify=F) # Get Lon
  #Lon <- as.numeric(lapply(Long$Lon, `[[`, 1))
  #Lati <- aggregate(Lat~No, x_win, paste, simplify=F) # Get Lat
  #Lat <- as.numeric(lapply(Lati$Lat, `[[`, 1))
  
  # merge summer and winter data frames
  normals <- full_join(smr, wnt, by="No")
  #head(normals)
  cond <- which(!complete.cases(normals))
  normals[cond,]
  return(normals[order(normals$No),])
}


#--------------------------------  
### Calculate Sigma Dissimilarity ####
#--------------------------------
#initiate the data frame to store the projected sigma dissimilarity of best analogs for each grid cell. 

loop_sigma_D <- function(A, B, C, append=""){
  NN.sigma <- data.frame(No=A$No,NN.sigma=NA, NN.station=NA, NN.Mdist=NA)
  
  for(j in 1:nrow(NN.sigma)){ 
    NN.sigma[j,2:4] <- calc_sigma_D(A, B, C, NN.sigma$No[j])
    if(j%%10==0){print(c(j, "of", nrow(NN.sigma)))}
  }
  names(NN.sigma)[2:4] <- paste0(names(NN.sigma_today_2100_4.5)[2:4],append)
  return(NN.sigma)
}

  
calc_sigma_D <- function(A, B, C, whichStation){
  # A is past climate
  # B is future climate
  # C is the data frame used to calculate the ICV
  
  if(!identical(dim(A), dim(B))){print("Error A and B different dimensions")}
  if(ncol(A) != ncol(C)){print("Error A and C different number of columns")}
  
  C.id <- C$No
  proxy <- B$No
  length(proxy)
  proxy2 <- sort(unique(proxy))
  if(!identical(proxy, proxy2)){break}
  
  # Principal component truncation rule
  trunc.SDs <- 0.1 #truncation 

  j <- whichStation
      
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
    A.prime <- sweep(A[,-1],MARGIN=2,Cj.sd[-1],`/`) #standardize the reference ICV
    # a <- matrix(c(1,2,3,4,5,6), nrow=2)
    # sweep(a, MARGIN =2, STATS=c(2,3,4)) # subtracts STATs from each column
    # sweep(a, MARGIN =2, STATS=c(2,3,4), FUN=`/`) # divides each column by STATS
    Bj.prime <- sweep(Bj[,-1],MARGIN=2,Cj.sd[-1],`/`) #standardize the analog pool    
    Cj.prime <- sweep(Cj[,-1],MARGIN=2,Cj.sd[-1],`/`) #standardize the projected future conditions of grid cells represented by ICV proxy j
    
    colnames(Cj.prime) <- colnames(A.prime)
    ## Step 2: Extract the principal components (PCs) of the reference period ICV 
    # and project all data onto these PCs
    PCA <- prcomp(Cj.prime[!is.na(apply(Cj.prime,1,mean)),])   
    # Principal components analysis. The !is.na(apply(...)) term is there 
    # simply to select all years with complete observations in all variables. 
    PCA$rotation
    
    #plot(PCA$rotation[,1], PCA$rotation[,2], xlim=c(-0.6, 0.1))
    #text(PCA$rotation[1:4,1], PCA$rotation[1:4,2], rownames(PCA$rotation)[1:4])
    # SST right of PC space, 
    # Arag and Calc and pH in upper left of PC space
    
    #plot(PCA$rotation[,1], PCA$rotation[,3], xlim=c(-0.6, 0.1))
    #text(PCA$rotation[1:4,1], PCA$rotation[1:4,3], rownames(PCA$rotation)[1:4])
    # separates pH from Calc and Arag
    
    #plot(PCA$rotation[,1], PCA$rotation[,4], xlim=c(-0.6, 0.1))
    #text(PCA$rotation[1:4,1], PCA$rotation[1:4,4], rownames(PCA$rotation)[1:4])
    # separates Calc from Arag
    
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
    nnd <- get.knnx(data=X.prime[,1:PCs],
                    query=Yj.prime[,1:PCs],
                    k=1,algorithm="brute")
    NN.dist <- as.vector(nnd[[2]]) 
    # Euclidean nearest neighbour distance in the z-standardized PCs of 
    # interannual climatic variability, i.e. the Mahalanobian nearest neighbour. 
    NN.chi <- pchi(NN.dist,PCs) # percentile of the nearest neighbour 
    # distance on the chi distribution with degrees of freedom 
    # equaling the dimensionality of the distance measurement (PCs)
    if(NN.chi>=1){NN.chi=1-1e-16}
      # sometimes with rounding error the NN.chi is greater than 1
      # if NN.chi equals 1 than NN.sigma is infinite
      # this slight transformation gives the largest possible value of NN.sigma
    NN.sigma <- qchi(NN.chi,1) 
    # values of the chi percentiles on a standard half-normal distribution (chi distribution with one degree of freedom)
  
    NN.station <- A$No[as.vector(nnd[[1]])]
    
    return(data.frame(NN.sigma, NN.station, NN.Mdist=NN.dist))
}


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


## Function to create a random grid empty grid 
# n = number of cells, increase n to increase resolution
makeGrid<-function(EB2){
  repeat{
    gr <- as.data.frame(spsample(EB2, 'regular', n  = 50000))
    names(gr) <- c('X', 'Y')
    coordinates(gr) <- c('X', 'Y')
    gridded(gr) <- TRUE  
    fullgrid(gr) <- TRUE  # Create SpatialGrid object
    try(proj4string(gr) <- proj4string(EB2))
    if(is.na(proj4string(gr))==FALSE) return(gr)
  }
}

#Function to convert raster objects as tibbles, written by Sébastien Rochette
gplot_data <- function(x, maxpixels = 100000)  {
  x <- raster::sampleRegular(x, maxpixels, asRaster = TRUE)
  coords <- raster::xyFromCell(x, seq_len(raster::ncell(x)))
  ## Extract values
  dat <- utils::stack(as.data.frame(raster::getValues(x))) 
  names(dat) <- c('value', 'variable')
  
  dat <- dplyr::as.tbl(data.frame(coords, dat))
  
  if (!is.null(levels(x))) {
    dat <- dplyr::left_join(dat, levels(x)[[1]], 
                            by = c("value" = "ID"))
  }
  dat
}
