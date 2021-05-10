# plot sd

-#--------------------------------  
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
identical(stationInfo$stations, sd_stations$No)
stationInfo <- data.frame(stationInfo, sd_stations)

#--------------------------------  
### plot sd ####
#--------------------------------
world <- map_data("world2") #Read in reference map from ggplot2
Plot_nonInt(stationInfo$lat, stationInfo$long, 
            stationInfo$SST_sum, world, "sd SST")

Plot_nonInt(stationInfo$lat, stationInfo$long, 
            stationInfo$Arag_sum, world, "sd Arag")


Plot_nonInt(stationInfo$lat, stationInfo$long, 
            stationInfo$Calc_sum, world, "sd Calc")

identical(stationInfo$No, 1:nrow(stationInfo))
# the code below works because the code above is true
Plot_nonInt(stationInfo$lat[norm_1800$No], 
            stationInfo$long[norm_1800$No], 
            norm_1800$SST_sum, world, 'mean sum 1800-SST')

Plot_nonInt(stationInfo$lat[norm_1800$No], 
            stationInfo$long[norm_1800$No], 
            norm_1800$SST_win, world, 'mean win 1800-SST')
