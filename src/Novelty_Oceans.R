### KE Lotterhos
### Oct 2018 - June 2019
### Northeastern University
### Mod. Áki - Nov. 2018 - June 2019

## Based on basic code for calculating and mapping climatic novelty using the sigma dissimilarity metric. 
## written by Colin Mahony
## Supp Material of "A closer look at novel climates: new methods and insights at continental to landscape scales"

# Current Plans
# Compare 1800 analog (A) to today surface (B), with today's surface variability (C)
# Compare today surface analog (A) to 2100-RCP8.5 (B), with today's surface variability (C)
# Compare today surface analog (A) to 2100-RCP4.5 (B), with today's surface variability (C)

# Future Plans
# Compare today all depths analog (A) to 2100 (B), with today's surface variability (C)

##################################
#### System setup ####
##################################
## Specify location of data ####
  setwd("/Users/katie/Desktop/Repos/OceanClimateNovelty/") 
  #setwd("~/Desktop/PostDoc/OceanClimateNovelty/")

  source("src/Novelty_Oceans_Functions.R")

##################################
#### Read in the input data ####
##################################

dat8.5 <- fread("data/large_files/T_Ar_Ca_pH_RCP85.txt", sep = ",")
dat4.5 <- fread("data/large_files/T_Ar_Ca_pH_RCP45.txt", sep = ",")

head(dat8.5)
head(dat4.5)
unique(dat8.5$Year)
unique(dat4.5$Year)

#boxplot(dat8.5$pH ~ dat8.5$Year)
#boxplot(dat4.5$pH ~ dat4.5$Year)
# note that these two datasets are exactly the same for years < 2000


#--------------------------------
#### 40-year climate normals ####
#--------------------------------
  # for example, 1930 represents from 1/1/1925 to 12/31/1934.  

  dat_1800 <- dat8.5 %>% filter(Year<1850)
  dim(dat_1800)
  dat_2000 <-   dat8.5 %>% filter(Year>1960 & Year<2010)
  dim(dat_2000)  
  dat_2100_8.5 <-   dat8.5 %>% filter(Year>2050)
  dim(dat_2100_8.5) 
  dat_2100_4.5 <-   dat4.5 %>% filter(Year>2050)
  dim(dat_2100_4.5) 


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
#--------------------------------
  length(norm_1800$No)
  length(norm_2000$No)
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
  (whichcols <- which(names(dat_2000) %in% c("SST", "Arag", "Calc", "pH")))
  C <- data.frame(dat_2000[,c(1, whichcols)], dat_2000[,c( whichcols)])
  head(C)
  
  NN.sigma_1800_today <- loop_sigma_D(A, B, C)
  dim(A)
  dim(B)
  head(NN.sigma_1800_today)
  
  final_dat <- full_join(NN.sigma_1800_today, stationInfo, by="No")
  head(final_dat)
  
  final_dat[which(is.infinite(final_dat$NN.sigma_1800_today)),]
  
  final_dat[which(!complete.cases(final_dat)),]
  
  #NN.sigma_1800_today[which(is.infinite(NN.sigma_1800_today))] <- NA
 
  head(final_dat)
  identical(nrow(final_dat), nrow(stationInfo)) 
    # should be true, same number as stations as stationInfo
  dim(B)
  
  # Visualize
  world <- map_data("world2")
  Plot_nonInt(final_dat$lat, final_dat$long, 
              final_dat$NN.sigma_1800_today, world, "sigma dis.")
  
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
  NN.sigma_today_2100_8.5 <- loop_sigma_D(A1, B1, C)
  names(NN.sigma_today_2100_8.5)[2:4] <- paste0(names(NN.sigma_today_2100_8.5)[2:4],"_today_2100_8.5")
  
  which(is.infinite(NN.sigma_today_2100_8.5$NN.sigma_today_2100_8.5))
  which(is.na(NN.sigma_today_2100_8.5$NN.sigma_today_2100_8.5))
  
  final_dat2 <- full_join(NN.sigma_today_2100_8.5, final_dat, by="No")
  head(final_dat2)
  
  dim(final_dat2)
  dim(stationInfo)
  
  # Visualize
  Plot_nonInt(final_dat2$lat, final_dat2$long, 
              final_dat2$NN.sigma_today_2100_8.5, world, "sigma dis.")
  
  
#--------------------------------  
### Today analog to 2100 RCP 4.5 ####
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
  NN.sigma_today_2100_4.5 <- loop_sigma_D(A1, B2, C)

  which(is.infinite(NN.sigma_today_2100_4.5$NN.sigma_today_2100_4.5))
  which(is.na(NN.sigma_today_2100_4.5$NN.sigma_today_2100_4.5))
  
  
  final_dat3 <- full_join(NN.sigma_today_2100_4.5, final_dat2, by="No")
  head(final_dat3)
  
  dim(final_dat3)
  dim(stationInfo)
  
  
  cond <- which(!complete.cases(final_dat3))
  length(cond)
  final_dat3[cond,]
  
  # Visualize
  Plot_nonInt(final_dat3$lat, final_dat3$long, 
              final_dat3$NN.sigma_today_2100_4.5, world, "sigma dis.")
  
  
### Write to file ####
write.csv(final_dat2, "results/SigmaD.csv", row.names = FALSE)
  

# Figure out NAs and NANs ###
  stations 671 NA
  station 773 NaN
  
  

  
  
  
  
  
  
  
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
