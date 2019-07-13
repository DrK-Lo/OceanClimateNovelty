### KE Lotterhos
### Oct 2018 - July 2019
### Northeastern University
### Mod. Áki - Nov. 2018 - June 2019

## Based on basic code for calculating and mapping climatic novelty using the sigma dissimilarity metric. 
## written by Colin Mahony
## Supp Material of "A closer look at novel climates: new methods and insights at continental to landscape scales"

# What this code does
# Compare 1800 analog (A) to today surface (B), with today's surface variability (C)
# Compare today surface analog (A) to 2100-RCP8.5 (B), with today's surface variability (C)
# Compare today surface analog (A) to 2100-RCP4.5 (B), with today's surface variability (C)

# This code also compares the converse (e.g. the current climate as B to the future climate as A)
# to assess disappearing climates

##################################
#### System setup ####
##################################
## Specify location of data ####
  setwd("/Users/katie/Desktop/Repos/OceanClimateNovelty/") 
  #setwd("~/Desktop/PostDoc/OceanClimateNovelty/")

  source("src/0-Novelty_Oceans_Functions.R")

##################################
#### Read in the input data ####
##################################
dat1 <- fread("/Users/katie/Google Drive/katie_research/2018-OceanClimateNovelty/Data/Lotterhos/Katie_Temp_Arag_1800_2000.txt", sep = ",")
dat4.5 <- fread("/Users/katie/Google Drive/katie_research/2018-OceanClimateNovelty/Data/Lotterhos/Katie_Temp_Arag_2070_2100_RCP45.txt", sep = ",")
dat8.5 <- fread("/Users/katie/Google Drive/katie_research/2018-OceanClimateNovelty/Data/Lotterhos/Katie_Temp_Arag_2070_2100_RCP85.txt", sep = ",")

head(dat1)
head(dat8.5)
head(dat4.5)
unique(dat1$Year)
unique(dat8.5$Year)
unique(dat4.5$Year)



#boxplot(dat8.5$pH ~ dat8.5$Year)
#boxplot(dat4.5$pH ~ dat4.5$Year)
# note that these two datasets are exactly the same for years < 2000


#--------------------------------
#### 40-year climate normals ####
#--------------------------------
  # for example, 1930 represents from 1/1/1925 to 12/31/1934.  

  dat_1800 <- dat1 %>% filter(Year<1850)
  dim(dat_1800)
  dat_2000 <-   dat1 %>% filter(Year>1960 & Year<2010)
  dim(dat_2000)  
  dat_2100_8.5 <-   dat8.5 
  dim(dat_2100_8.5) 
  dat_2100_4.5 <-   dat4.5 
  dim(dat_2100_4.5) 

  sum(!complete.cases(dat_1800))
  sum(!complete.cases(dat_2000))
  sum(!complete.cases(dat_2100_4.5))
  sum(!complete.cases(dat_2100_8.5))
  hist(dat_1800$Month)
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
  norm_2100_8.5 <- calculate_normals(dat_2100_8.5)
  norm_2100_4.5 <- calculate_normals(dat_2100_4.5)

  head(norm_1800)
  dim(norm_1800)
  dim(norm_2000)
  dim(norm_2100_8.5)
  dim(norm_2100_4.5)

#--------------------------------  
### data frame to link station number to Lat Long ####
#--------------------------------
  head(dat8.5)
  stations <- unique(dat8.5$No)
  stationInfo = data.frame(stations=stations, lat=NA, long=NA)
  unik <- which(!duplicated(dat8.5$No))
  head(unik)
  stationInfo$lat <- dat8.5$Lat[unik]
  stationInfo$long <- dat8.5$Lon[unik]
  names(stationInfo)[1] <- "No"
  head(stationInfo)
  which(!complete.cases(stationInfo))


#--------------------------------  
### 1800 analog to today ####
# What are today's novel climates compared to 1800?
# Novel climates are identified by comparing the 
# future climate realization (B) for each gridpoint 
# to the past (A) climate realizations (Williams)
#--------------------------------

  A <- norm_1800
  # 1800 climate normals
  head(A)
  dim(A)
  
  # subset the 2000 data to the same stations as the 1800 data
  B <- norm_2000[norm_2000$No %in% norm_1800$No,]
  # 2000 climate normals 
  head(B)
  dim(B)
  # same columns as A
  
  # sanity check to make sure stations in right order
  identical(A$No, B$No) # should be true
  
  head(dat_2000)  
  (whichcols <- which(names(dat_2000) %in% c("SST", "Arag", "pH")))
  C <- data.frame(dat_2000[,c(1, whichcols)], dat_2000[,c( whichcols)])
  head(C)
  
  NN.sigma_1800_today <- loop_sigma_D(A, B, C, "_1800_today")
  head(NN.sigma_1800_today)
  
  final_dat <- full_join(NN.sigma_1800_today, stationInfo, by="No")
  head(final_dat)
  
  final_dat[which(is.infinite(final_dat$NN.sigma_1800_today)),]
  
  final_dat[which(!complete.cases(final_dat)),]
  
  #NN.sigma_1800_today[which(is.infinite(NN.sigma_1800_today))] <- NA
 
  identical(nrow(final_dat), nrow(stationInfo)) 
    # should be true, same number as stations as stationInfo
  
  # Visualize
  world <- map_data("world2")
  Plot_nonInt(final_dat$lat, final_dat$long, 
              final_dat$NN.sigma_1800_today, world, "sigma dis.")
  
  # The NN.sigma_1800_today represents the novelty of today's
    # climate compared to 1800
  # NN.dist_1800_today represents the Mahalanobis distance from the climate
    # in today to it's nearest neighbor in all 1800 data
  # No--> NN.station_1800_today represents the station that was queried 2000 (B) 
      # and it's nearest data in the global reference NN.station_1800_today (1800, A)
  
#--------------------------------  
### Today analog to 1800 ####
# What are 1800 disappearing climates?
# disappearing climates are identified by comparing
# each past gridpoint (A) to all future climate realizations (B)
#--------------------------------
  calc_sigma_D(A, B, C, 1)
  calc_sigma_D(B, A, C, 1)

  B <- norm_1800[complete.cases(norm_1800),]
  A <- norm_2000[which(norm_2000$No %in% B$No),]
  identical(A$No, B$No) # should be true
  identical(sort(A$No), A$No) # should be true
  NN.sigma_today_1800 <- loop_sigma_D(A, B, C, "_today_1800")
  
  head(NN.sigma_today_1800)
  
  final_dat2 <- full_join(final_dat, NN.sigma_today_1800)
  head(final_dat2)
  dim(final_dat2)
  Plot_nonInt(final_dat2$lat, final_dat2$long, 
              final_dat2$NN.sigma_today_1800, world, "sigma dis.")
  
  ggplot(final_dat2) + geom_point(aes(x=NN.sigma_1800_today, 
                                      y=NN.sigma_today_1800,
                                      color=lat), alpha=0.5) +
    scale_color_gradient2(low="red", high="blue", mid="grey") +
    theme_classic() + xlab("Novel climates today from 1800") +
    ylab("Disappearing 1800 climates")
  
  # The NN.sigma_1800_disappear represents the novelty of 1800
    # compared to today (e.g disappearing climates)
  # NN.dist_1800_disappear represents the Mahalanobis distance from the climate
  # in 1800 to it's nearest neighbor in all of today's data
  # No--> NN.station_1800_disappear represents the the station that was queried 1800 (B) 
  # and it's nearest data in the global reference NN.station_1800_disappear (2000, A)
  
  # Add the info for lat/long of nearest neighbors
  stationInfo2 <- stationInfo
  names(stationInfo2) <- paste0(names(stationInfo2),"_1800_today")
  head(stationInfo2)
  dim(stationInfo)
  names(stationInfo2)[1] <-"NN.station_1800_today"
  final_dat3 <- left_join(final_dat2, stationInfo2)
  dim(final_dat3)
  head(final_dat3)
  tail(final_dat3)
  
  
  stationInfo3 <- stationInfo
  names(stationInfo3) <- paste0(names(stationInfo3),"_today_1800")
  names(stationInfo3)[1] <- "NN.station_today_1800"
  final_dat4 <- left_join(final_dat3, stationInfo3)
  dim(final_dat4)
  head(final_dat4)
  
  # cond <- which(final_dat4$NN.sigma_1800_today>5 |  final_dat4$NN.sigma_1800_disappear>5)
  # length(cond)
  # final_dat4[cond,]
  
  ggplot(final_dat4) + geom_point(aes(y=lat, 
                                      x=lat_1800_today,
                                      color=NN.sigma_1800_today), alpha=0.2) +
    scale_color_gradient2(low="red", high="blue", mid="grey") +
    theme_classic() + ylab("Latitude in 2000") +
    xlab("Latitude of nearest neighbor in 1800") + geom_abline(intercept=0,slope=1) +
    geom_abline(intercept=0,slope=-1)
  # This plot makes more sense in the way we typically think about it.
  # If we consider a location in 2000, where was the climate it came from in 1800?
  
  
  
  # ggplot(final_dat4) + geom_point(aes(x=lat_1800_disappear,
  #                                     y=lat,
  #                                     color=NN.sigma_1800_disappear), alpha=0.2) +
  #   scale_color_gradient2(low="red", high="blue", mid="grey") +
  #   theme_classic() + ylab("Latitude in 1800") +
  #   xlab("Latitude of nearest neighbor in 2000") + geom_abline(intercept=0,slope=1) +
  #   geom_abline(intercept=0,slope=-1)
  # this is weird
  
  #n <- final_dat4 %>% filter(lat_1800_disappear > 0 & lat_1800_disappear < 10)
  #Plot_nonInt(n$lat_1800_disappear, n$long_1800_disappear, 
  #            n$NN.sigma_1800_disappear, world, "sigma dis.")
  hist(final_dat4$lat_1800_disappear, breaks=seq(-90,90, by=1))
  hist(final_dat4$lat_1800_today, breaks=seq(-90,90, by=1))
  hist(final_dat4$lat, breaks=seq(-90,90, by=1))
  
#--------------------------------  
### Today analog to 2100 RCP 8.5 ####
#--------------------------------
  identical(norm_2000$No, norm_2100_8.5$No)

  A1 <- norm_2000
  # 1970-2000 climate normals
  head(A1)
  dim(A1)

  B1 <- norm_2100_8.5
  # 2070-2100
  head(B1)
  dim(B1)
  # same columns as A

  # sanity check to make sure stations in right order
  identical(A1$No, B1$No) # should be true
  identical(A1$No, sort(A1$No)) # should also be true

  # C is the same defined above
  NN.sigma_today_2100_8.5 <- loop_sigma_D(A1, B1, C, "_today_2100_8.5")
  
  which(is.infinite(NN.sigma_today_2100_8.5$NN.sigma_today_2100_8.5))
  which(is.na(NN.sigma_today_2100_8.5$NN.sigma_today_2100_8.5))
  
  final_dat5 <- full_join(NN.sigma_today_2100_8.5, final_dat4, by="No")
  head(final_dat5)
  
  dim(final_dat5)
  stationInfo5 <- stationInfo
  names(stationInfo5) <- paste0(names(stationInfo),"_today_2100_8.5")
  names(stationInfo5)[1] <- "NN.station_today_2100_8.5"
  head(stationInfo5)
  final_dat6 <- left_join(final_dat5, stationInfo5)
  dim(final_dat6)
  head(final_dat6)
  
  # Visualize
  Plot_nonInt(final_dat6$lat, final_dat6$long, 
              final_dat6$NN.sigma_today_2100_8.5, world, "sigma dis.")
  
#--------------------------------  
### 2100 RCP 8.5 to today, what are today's climates that ####
### will disappear in 2100 RCP 8.5?
#--------------------------------
  NN.sigma_2100_8.5_today <- loop_sigma_D( B1, A1, C, "_2100_8.5_today")
  final_dat7 <- full_join(NN.sigma_2100_8.5_today, final_dat6, by="No")
  head(final_dat7)
  Plot_nonInt(final_dat7$lat, final_dat7$long, 
              final_dat7$NN.sigma_2100_8.5_today, world, "sigma dis.")
  
  stationInfo7 <- stationInfo
  names(stationInfo7) <- paste0(names(stationInfo),"_2100_8.5_today")
  names(stationInfo7)[1] <- "NN.station_2100_8.5_today"
  head(stationInfo7)
  final_dat8 <- left_join(final_dat7, stationInfo7)
  dim(final_dat8)
  head(final_dat8)
  
#--------------------------------  
### 2100 RCP 4.5 against today (novelty) ####
#--------------------------------
  identical(norm_2000$No, norm_2100_4.5$No)
  
  # A1 the same as above
  
  B2 <- norm_2100_4.5
  # 2070-2100
  head(B2)
  dim(B2)
  # same columns as A
  
  # sanity check to make sure stations in right order
  identical(A1$No, B2$No) # should be true
  identical(A1$No, sort(A1$No)) # should also be true
  
  # C is the same defined above
  NN.sigma_today_2100_4.5 <- loop_sigma_D(A1, B2, C, "_today_2100_4.5")

  which(is.infinite(NN.sigma_today_2100_4.5$NN.sigma_today_2100_4.5))
  which(is.na(NN.sigma_today_2100_4.5$NN.sigma_today_2100_4.5))
  
  
  final_dat9 <- full_join(NN.sigma_today_2100_4.5, final_dat8, by="No")
  head(final_dat9)
  

  # Visualize
  Plot_nonInt(final_dat9$lat, final_dat9$long, 
              final_dat9$NN.sigma_today_2100_4.5, world, "sigma dis.")
  
  
  stationInfo9 <- stationInfo
  names(stationInfo9) <- paste0(names(stationInfo),"_today_2100_4.5")
  names(stationInfo9)[1] <- "NN.station_today_2100_4.5"
  head(stationInfo9)
  final_dat10 <- left_join(final_dat9, stationInfo9)
  dim(final_dat10)
  head(final_dat10)
  
#--------------------------------  
### today against 2100 RCP 4.5 (disappearing climate) ####
#--------------------------------
  NN.sigma_2100_4.5_today <- loop_sigma_D(B2, A1, C, "_2100_4.5_today")
  which(is.infinite(NN.sigma_2100_4.5_today$NN.sigma_2100_4.5_today))
  which(is.na(NN.sigma_2100_4.5_today$NN.sigma_2100_4.5_today))
  
  
  final_dat11 <- full_join(NN.sigma_2100_4.5_today, final_dat10, by="No")
  head(final_dat11)
  
  # Visualize
  Plot_nonInt(final_dat11$lat, final_dat11$long, 
              final_dat11$NN.sigma_2100_4.5_today, world, "sigma dis.")
  
  stationInfo11 <- stationInfo
  names(stationInfo11) <- paste0(names(stationInfo),"_2100_4.5_today")
  names(stationInfo11)[1] <- "NN.station_2100_4.5_today"
  head(stationInfo11)
  final_dat12 <- left_join(final_dat11, stationInfo11)
  dim(final_dat12)
  head(final_dat12)
    
  sum(!complete.cases(final_dat12))

  
  
### Write to file ####
write.csv(final_dat12, "results/SigmaD.csv", row.names = FALSE)
  
dat_Nov <- final_dat12
for_liqing <- data.frame(No=dat_Nov$No, lat=dat_Nov$lat,long=dat_Nov$long,A=dat_Nov$NN.sigma_today_1800, B=dat_Nov$NN.sigma_1800_today, C=dat_Nov$NN.sigma_2100_4.5_today, D=dat_Nov$NN.sigma_today_2100_4.5, E=dat_Nov$NN.sigma_2100_8.5_today, F=dat_Nov$NN.sigma_today_2100_8.5)
write.csv(for_liqing,  "../results/SigmaD_plot1.csv", row.names=FALSE)

for_liqing_Md <- data.frame(No=dat_Nov$No, lat=dat_Nov$lat,long=dat_Nov$long,A=dat_Nov$NN.Mdist_today_1800, B=dat_Nov$NN.Mdist_1800_today, C=dat_Nov$NN.Mdist_2100_4.5_today, D=dat_Nov$NN.Mdist_today_2100_4.5, E=dat_Nov$NN.Mdist_2100_8.5_today, F=dat_Nov$NN.Mdist_today_2100_8.5)
write.csv(for_liqing_Md,  "../results/SigmaD_plot2.csv", row.names=FALSE)
  
  
plot(final_dat12$lat_today_2100_8.5, final_dat12$lat) #If A is today and B is future, 
#  where station No's climate in the future will come from today (or the closest similar climate today).
# Latitude of the station today (x-axis) where the latitude of the station in the future (y-axis) will come from


plot(final_dat12$lat, final_dat12$lat_2100_8.5_today) # If A is future and B is today, this 
# represents where station "No" will be found in the future (or the closest similar climate).
# Latitude of the station today (x-axis) and where it will be in the future (y-axis)


  
  
  
  
  
  
  