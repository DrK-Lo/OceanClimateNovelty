### Novelty_Oceans Functions
### KE Lotterhos
### Oct 2018 - mod. June 2019
### Northeastern University
### Mod. Áki Nov. 2018


#Create function that removes previous user installed packages to avoid masking
clean_pkgs<-function(){
  lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE)
}
clean_pkgs() #Remove all non-essential previously called packages

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
