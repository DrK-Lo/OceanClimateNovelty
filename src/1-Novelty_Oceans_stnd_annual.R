### KE Lotterhos
### Oct 2018 - Oct 2019
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

# This code uses the annnual climate normals and monthly ICV for the novelty calculation

##################################
#### System setup ####
##################################
## Specify location of data ####
setwd("/Users/lotterhos/Documents/GitHub/OceanClimateNovelty") 
#setwd("~/Desktop/PostDoc/OceanClimateNovelty/")

source("src/0-Novelty_Oceans_Functions.R")

##################################
#### Read in the input data ####
##################################
dat1 <- fread("/Users/lotterhos/Google Drive/katie_research/2018-OceanClimateNovelty/Data/Lotterhos/Katie_Temp_Arag_1800_2000.txt", sep = ",")
dat4.5 <- fread("/Users/lotterhos/Google Drive/katie_research/2018-OceanClimateNovelty/Data/Lotterhos/Katie_Temp_Arag_2070_2100_RCP45.txt", sep = ",")
dat8.5 <- fread("/Users/lotterhos/Google Drive/katie_research/2018-OceanClimateNovelty/Data/Lotterhos/Katie_Temp_Arag_2070_2100_RCP85.txt", sep = ",")

ICV <- fread("/Users/lotterhos/Google Drive/katie_research/2018-OceanClimateNovelty/Data/ICV_Temp_Arag_pH_1960_2020.txt")

head(dat1)
head(dat8.5)
head(dat4.5)
head(ICV)
unique(dat1$Year)
unique(dat8.5$Year)
unique(dat4.5$Year)
unique(ICV$Year)

hist(c(dat1$Arag, dat8.5$Arag))
hist(log10(c(dat1$Arag, dat8.5$Arag)))

##################################
#### Log-transform Arag ####
##################################
dat1$Arag <- log10(dat1$Arag)
dat4.5$Arag <- log10(dat4.5$Arag)
dat8.5$Arag <- log10(dat8.5$Arag)
ICV$Arag <- log10(ICV$Arag)

#  Aragonite saturation is a ratio variable; it is limited at zero and 
# proportional changes are meaningful. The analysis uses raw values of 
# aragonite saturation state. Under raw scaling, the difference between 1 and 2 
# is equal to the difference between zero and one, which doesn’t reflect the 
# meaning of the variable. Instead, the difference between 1 and 2 should 
# be given the same significance as the difference between 0.5 and 1 (doubling vs halving). 
# Log-transforming aragonite saturation state will produce this more meaningful 
# proportional scaling.

# More generally on the same topic: The need for proportional scaling is the 
# reason why precipitation and other ratio variables (e.g. degree-days) are 
# typically log-transformed in climate space analysis. In theory, all variables 
# should be log-transformed to produce an environmental space where distances 
# are comparable for different locations; however, in practice temperature 
# doesn’t need to be log-transformed because it doesn’t vary across orders 
# of magnitude, and in our case pH is already a log-scaled variable.

#boxplot(dat8.5$pH ~ dat8.5$Year)
#boxplot(dat4.5$pH ~ dat4.5$Year)
# note that these two datasets are exactly the same for years < 2000


#--------------------------------
#### 40-year climate normals ####
#--------------------------------
# for example, 1930 represents from 1/1/1925 to 12/31/1934.  
dat_2000 <-   dat1 %>% filter(Year>1960 & Year<2010)
dim(dat_2000)
# Below is a standardization code that can be deleted becuase the PCA
# does the standardization
# x_SST <- mean(dat_2000$SST)
# s_SST <- sd(dat_2000$SST)
# x_Arag <- mean(dat_2000$Arag)
# s_Arag <- sd(dat_2000$Arag)
# x_pH <- mean(dat_2000$pH)
# s_pH <- sd(dat_2000$pH)
# 
# head(dat_2000)  
#   dat_2000$SST <- (dat_2000$SST - x_SST)/s_SST
#   dat_2000$Arag <- (dat_2000$Arag - x_Arag)/s_Arag
#   dat_2000$pH <- (dat_2000$pH - x_pH)/s_pH
#   dat_2000 %>% summarise_each(sd)
# 
dat_1800 <- dat1 %>% filter(Year<1850)
dim(dat_1800)
#   dat_1800$SST <- (dat_1800$SST - x_SST)/s_SST
#   dat_1800$Arag <- (dat_1800$Arag - x_Arag)/s_Arag
#   dat_1800$pH <- (dat_1800$pH - x_pH)/s_pH
#   
#   
dat_2100_8.5 <-   dat8.5 
#     dim(dat_2100_8.5) 
#     dat_2100_8.5$SST <- (dat_2100_8.5$SST - x_SST)/s_SST
#     dat_2100_8.5$Arag <- (dat_2100_8.5$Arag - x_Arag)/s_Arag
#     dat_2100_8.5$pH <- (dat_2100_8.5$pH - x_pH)/s_pH
#   
dat_2100_4.5 <-   dat4.5 
#   dim(dat_2100_4.5) 
#   dat_2100_4.5$SST <- (dat_2100_4.5$SST - x_SST)/s_SST
#   dat_2100_4.5$Arag <- (dat_2100_4.5$Arag - x_Arag)/s_Arag
#   dat_2100_4.5$pH <- (dat_2100_4.5$pH - x_pH)/s_pH

sum(!complete.cases(dat_1800))
sum(!complete.cases(dat_2000))
sum(!complete.cases(dat_2100_4.5))
sum(!complete.cases(dat_2100_8.5))
sum(!complete.cases(ICV))

norm_1800 <- calculate_normals_annual(dat_1800)
norm_2000 <- calculate_normals_annual(dat_2000)
norm_2100_8.5 <- calculate_normals_annual(dat_2100_8.5)
norm_2100_4.5 <- calculate_normals_annual(dat_2100_4.5)

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
### ICV matrix
head(ICV)  
C <- ICV %>% select(No, SST, Arag, pH)
head(C)
#--------------------------------   

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
B <- norm_2000
# 2000 climate normals 
head(B)
dim(B)
# same columns as A

# sanity check to make sure stations in right order
identical(A$No, B$No) # should be true

# check list 18906,18952, 29736, 29789, 8193, 29319,
# calc_sigma_D(A, B, C,  18906, "_1800_today")
# calc_sigma_D(A, B, C,  18952, "_1800_today")
# calc_sigma_D(A, B, C,  29736, "_1800_today")
# calc_sigma_D(A, B, C,  29789, "_1800_today")
# calc_sigma_D(A, B, C,  8193, "_1800_today")
# calc_sigma_D(A, B, C,  29319, "_1800_today")

NN.sigma_1800_today <- loop_sigma_D(A, B, C, "_1800_today")
head(NN.sigma_1800_today)
hist(NN.sigma_1800_today$numPCs_1800_today)

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
            final_dat$NN.sigma_1800_today, world, "sigD novel 2000")

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

B <- norm_1800
A <- norm_2000
identical(A$No, B$No) # should be true
identical(sort(A$No), A$No) # should be true


# check list 18906,18952, 29736, 29789, 8193, 29319,
# calc_sigma_D(A, B, C,  18906, "_today_1800")
# calc_sigma_D(A, B, C,  18952, "_today_1800")
# calc_sigma_D(A, B, C,  29736, "_today_1800")
# calc_sigma_D(A, B, C,  29789, "_today_1800")
# calc_sigma_D(A, B, C,  8193, "_today_1800")
# calc_sigma_D(A, B, C,  29319, "_today_1800")


NN.sigma_today_1800 <- loop_sigma_D(A, B, C, "_today_1800")

head(NN.sigma_today_1800)

final_dat2 <- full_join(final_dat, NN.sigma_today_1800)
head(final_dat2)
dim(final_dat2)
Plot_nonInt(final_dat2$lat, final_dat2$long, 
            final_dat2$NN.sigma_today_1800, world, "sigD dis 2000")

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
#hist(final_dat4$lat_1800_disappear, breaks=seq(-90,90, by=1))
#hist(final_dat4$lat_1800_today, breaks=seq(-90,90, by=1))
#hist(final_dat4$lat, breaks=seq(-90,90, by=1))

#--------------------------------  
### Today analog to 2100 RCP 8.5 (Novelty) ####
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

# check list 18906,18952, 29736, 29789, 8193, 29319,
# calc_sigma_D(A1, B1, C,  18906, "_today_2100_8.5")
# calc_sigma_D(A1, B1, C,  18952, "_today_2100_8.5")
# calc_sigma_D(A1, B1, C,  29736, "_today_2100_8.5")
# calc_sigma_D(A1, B1, C,  29789, "_today_2100_8.5")
# calc_sigma_D(A1, B1, C,  8193, "_today_2100_8.5")
# calc_sigma_D(A1, B1, C,  29319, "_today_2100_8.5")

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
            final_dat6$NN.sigma_today_2100_8.5, world, "sigma nov 2100 8.5")

#--------------------------------  
### 2100 RCP 8.5 to today, what are today's climates that ####
### will disappear in 2100 RCP 8.5?
#--------------------------------
# check list 18906,18952, 29736, 29789, 8193, 29319,
# calc_sigma_D( B1, A1, C,  18906, "_2100_8.5_today")
# calc_sigma_D( B1, A1, C,  18952, "_2100_8.5_today")
# calc_sigma_D( B1, A1, C,  29736, "_2100_8.5_today")
# calc_sigma_D( B1, A1, C,  29789, "_2100_8.5_today")
# calc_sigma_D( B1, A1, C,  8193, "_2100_8.5_today")
# calc_sigma_D( B1, A1, C,  29319, "_2100_8.5_today")

NN.sigma_2100_8.5_today <- loop_sigma_D( B1, A1, C, "_2100_8.5_today")
final_dat7 <- full_join(NN.sigma_2100_8.5_today, final_dat6, by="No")
head(final_dat7)
Plot_nonInt(final_dat7$lat, final_dat7$long, 
            final_dat7$NN.sigma_2100_8.5_today, world, "sigma dis. 2000 in 2100 8.5")

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
            final_dat9$NN.sigma_today_2100_4.5, world, "sig nov 2100 4.5")


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
            final_dat11$NN.sigma_2100_4.5_today, world, "sigma dis. 2000 in 2100 4.5")

stationInfo11 <- stationInfo
names(stationInfo11) <- paste0(names(stationInfo),"_2100_4.5_today")
names(stationInfo11)[1] <- "NN.station_2100_4.5_today"
head(stationInfo11)
final_dat12 <- left_join(final_dat11, stationInfo11)
dim(final_dat12)
head(final_dat12)

sum(!complete.cases(final_dat12))


### Final dataframe ####
dat_Nov <- final_dat12
### Write to file ####
write.csv(final_dat12, "results_annual/SigmaD.csv", row.names = FALSE)

head(dat_Nov)
for_liqing <- data.frame(No=dat_Nov$No, lat=dat_Nov$lat,long=dat_Nov$long,
                         A=dat_Nov$NN.sigma_today_1800, 
                         B=dat_Nov$NN.sigma_1800_today, 
                         C=dat_Nov$NN.sigma_2100_4.5_today, 
                         D=dat_Nov$NN.sigma_today_2100_4.5, 
                         E=dat_Nov$NN.sigma_2100_8.5_today, 
                         F=dat_Nov$NN.sigma_today_2100_8.5)
write.csv(for_liqing,  "results_annual/SigmaD_plot1.csv", row.names=FALSE)

for_liqing_Md <- data.frame(No=dat_Nov$No, lat=dat_Nov$lat,long=dat_Nov$long,A=dat_Nov$NN.Mdist_today_1800, B=dat_Nov$NN.Mdist_1800_today, C=dat_Nov$NN.Mdist_2100_4.5_today, D=dat_Nov$NN.Mdist_today_2100_4.5, E=dat_Nov$NN.Mdist_2100_8.5_today, F=dat_Nov$NN.Mdist_today_2100_8.5)
write.csv(for_liqing_Md,  "results_annual/SigmaD_plot2.csv", row.names=FALSE)


#plot(final_dat12$lat_today_2100_8.5, final_dat12$lat) #If A is today and B is future, 
#  where station No's climate in the future will come from today (or the closest similar climate today).
# Latitude of the station today (x-axis) where the latitude of the station in the future (y-axis) will come from


#plot(final_dat12$lat, final_dat12$lat_2100_8.5_today) # If A is future and B is today, this 
# represents where station "No" will be found in the future (or the closest similar climate).
# Latitude of the station today (x-axis) and where it will be in the future (y-axis)


### Calculate percentages ####  
N <- nrow(dat_Nov)  

calc_perc <- function(x, descrip){
  print(paste("Percent of cells with moderate dissimilarity for", descrip,":", 
              round(sum(x>2 & x<4)/N,3)*100, sep=" "))
  print(paste("Percent of cells with extreme dissimilarity for", descrip,":", 
              round(sum(x>4)/N, 3)*100, sep=" "))
}

calc_perc(dat_Nov$NN.sigma_today_1800, "1800 disappearing in 2000")
calc_perc(dat_Nov$NN.sigma_1800_today, "Novel climates in 2000 since 1800")
calc_perc(dat_Nov$NN.sigma_2100_4.5_today, "2000 disappearing in 2100 4.5")
calc_perc(dat_Nov$NN.sigma_today_2100_4.5, "Novel climates in 2100 4.5 since 2000")
calc_perc(dat_Nov$NN.sigma_2100_8.5_today, "2000 disappearing in 2100 8.5")
calc_perc(dat_Nov$NN.sigma_today_2100_8.5, "Novel climates in 2100 8.5 since 2000")

# Visualize relationship between sigma_D and M_D
#### Figure out what is going on ####
plot(dat_Nov$NN.sigma_today_1800, dat_Nov$NN.Mdist_today_1800)
plot(dat_Nov$NN.sigma_today_2100_8.5, dat_Nov$NN.Mdist_today_2100_8.5, col=dat_Nov$numPCs_today_2100_8.5)
legend(0, 15, 1:4, 1:4)
abline(h=12)
hist(dat_Nov$numPCs_today_2100_8.5)

A <- norm_2000
B <- norm_2100_8.5

weird_ind <- which(dat_Nov$NN.sigma_today_2100_8.5> 6 & dat_Nov$NN.sigma_today_2100_8.5 < 8 & 
                     dat_Nov$NN.Mdist_today_2100_8.5 > 12)
dat_Nov[weird_ind,]
hist(dat_Nov$lat[weird_ind])
plot(dat_Nov$NN.sigma_today_2100_8.5[weird_ind], dat_Nov$NN.Mdist_today_2100_8.5[weird_ind], 
     col=as.numeric(as.factor(dat_Nov$lat[weird_ind]>50)))
weird_stations <- dat_Nov$No[weird_ind]
length(weird_stations)     

dat_Nov_weird <- dat_Nov[weird_ind,]
dat_Nov_weird %>% filter(NN.Mdist_today_2100_8.5>15) %>% select(No)
dat_Nov_weird %>% filter(NN.Mdist_today_2100_8.5<13.05) %>% select(No, lat, long)
#stations to start with in Novelty functions

norm_2000_weird <- norm_2000[norm_2000$No %in% weird_stations,]
norm_2001_8.5_weird <- norm_2100_8.5[norm_2100_8.5$No %in% weird_stations,]

plot(norm_2000$Arag_sum, norm_2000$Arag_win)
abline(0,1, col="blue")
points(norm_2000_weird$Arag_sum, norm_2000_weird$Arag_win, col="red")
abline(0,1, col="blue")

plot(norm_2000$Arag_sum,  norm_2100_8.5$Arag_sum)
abline(0,1, col="blue")
points(norm_2000_weird$Arag_sum, norm_2001_8.5_weird$Arag_sum, col="red")

left_join(norm_2000_weird , dat_Nov$lat)
