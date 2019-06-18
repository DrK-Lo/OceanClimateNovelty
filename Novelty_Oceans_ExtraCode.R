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