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


##Interpolation for visualization
B2a<-B2[!is.na(B2$NN.sigma),]

B2a<-B2a[,c(4,3,2)]

for(i in 1:nrow(B2a)){
  if(B2a$long[i]>360){
    B2a$long[i]<-B2a$long[i]-360
  }
}

EB2 <- SpatialPoints(B2a) # this is your spatial points df

# Project sp object to WGS 84
proj4string(EB2) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

# Create an empty grid; n = number of cells
# Increase n to increase resolution
gr <- as.data.frame(spsample(EB2, 'regular', n  = 50000))
names(gr) <- c('X', 'Y')
coordinates(gr) <- c('X', 'Y')
gridded(gr) <- TRUE  
fullgrid(gr) <- TRUE  # Create SpatialGrid object
proj4string(gr) <- proj4string(EB2) 
#If this line throws an error, run this
# chunk again. It is a random grid.

# Interpolate the grid cells using inverse distance weighing power = 8
# NN.sigma ~ 1 = simple kriging
EB2.idw <- idw(NN.sigma ~ 1, EB2, newdata = gr, idp = 8)

# Convert to raster
r <- raster(EB2.idw)

world<-map("world2", fill=T,plot=F)
y<-map2SpatialPolygons(world, IDs = sapply(strsplit(world$names, ":"), function(x) x[1]), proj4string=CRS("+proj=longlat +datum=WGS84"))
z<-st_as_sf(y)
wr <- raster(z, res = 0.01)
wrld_r <- fasterize(z, wr)
gplot_wrld_r <- gplot_data(wrld_r)

gplot_r <- gplot_data(r)

#Change scale_fill_gradient value to name of variable
ggplot() +
  geom_tile(data = gplot_r, 
            aes(x = x, y = y, fill = value)) +
  geom_tile(data = dplyr::filter(gplot_wrld_r, !is.na(value)), 
            aes(x = x, y = y), fill = "grey20") +
  xlim(0,360) +
  ylim(-77,90) + 
  xlab("Long") +
  ylab("Lat") +
  ggtitle("Sigma dissimilarity: Today from 1800") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_gradient2(expression(paste(sigma," dis.")),
                       low = 'blue', mid = "yellow", high = 'red',
                       midpoint = 4,
                       limits=c(0,8.1),
                       na.value = NA) +
  coord_quickmap()



#Visualize with interpolation
B2a<-B2[!is.na(B2$NN.sigma),]

B2a<-B2a[,c(4,3,2)]

for(i in 1:nrow(B2a)){
  if(B2a$long[i]>360){
    B2a$long[i]<-B2a$long[i]-360
  }
}

EB2 <- SpatialPoints(B2a) # this is your spatial points df

# Project sp object to WGS 84
proj4string(EB2) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

# Function to create a random grid empty grid
gr<-makeGrid(EB2)

# Interpolate the grid cells using power value = 2
# NN.sigma ~ 1 = simple kriging

EB2.idw <- idw(NN.sigma ~ 1, EB2, newdata = gr, idp = 8)

# Convert to raster
r <- raster(EB2.idw)

world<-map("world2", fill=T,plot=F)
y<-map2SpatialPolygons(world, IDs = sapply(strsplit(world$names, ":"), function(x) x[1]), proj4string=CRS("+proj=longlat +datum=WGS84"))
z<-st_as_sf(y)
wr <- raster(z, res = 0.01)
wrld_r <- fasterize(z, wr)
gplot_wrld_r <- gplot_data(wrld_r)

gplot_r <- gplot_data(r)

#Change scale_fill_gradient value to name of variable
ggplot() +
  geom_tile(data = gplot_r, 
            aes(x = x, y = y, fill = value)) +
  geom_tile(data = dplyr::filter(gplot_wrld_r, !is.na(value)), 
            aes(x = x, y = y), fill = "grey20") +
  #xlim(2,358) +
  ylim(-78,90) + 
  xlab("Long") +
  ylab("Lat") +
  ggtitle("Sigma dissimilarity: Today from 2100") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_gradient2(expression(paste(sigma," dis.")),
                       low = 'blue', mid = "yellow", high = 'red',
                       midpoint = 4,
                       limits=c(0,8.1),
                       na.value = NA) +
  coord_quickmap()
