### Novelty_Oceans Functions
### KE Lotterhos
### Oct 2018
### Northeastern University

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
  pH_sum <- tapply(x_sum$pH,INDEX = x_sum$No,mean, rm.na=TRUE)
  length(Calc_sum)
  identical(names(SST_sum), names(Arag_sum))
  identical(names(Arag_sum), names(Calc_sum))
  smr <- data.frame(No=as.integer(names(SST_sum)), SST_sum,Arag_sum,Calc_sum, pH_sum  )
  head(smr)
  
  ## Winter calculations    
  x_win <- dat1 %>% filter((Lat >0 & Month %in% month1)|
                             (Lat <0 & Month %in% month2))
  SST_win <- tapply(x_win$SST,INDEX = x_win$No,mean, rm.na=TRUE)
  Arag_win <- tapply(x_win$Arag,INDEX = x_win$No,mean, rm.na=TRUE)
  Calc_win <- tapply(x_win$Calc,INDEX = x_win$No,mean, rm.na=TRUE)
  pH_win <- tapply(x_win$pH,INDEX = x_win$No,mean, rm.na=TRUE)
  identical(names(SST_win), names(Arag_win))
  identical(names(Arag_win), names(Calc_win))
  wnt <- data.frame(No=as.integer(names(SST_win)), SST_win,Arag_win,Calc_win,pH_win  )
  head(wnt)
  
  # merge summer and winter data frames
  identical(names(SST_win), names(SST_sum))
  normals <- full_join(smr, wnt, by="No")
  head(normals)
  nrow(normals) == nrow(wnt)
  nrow(normals) == nrow(smr)
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