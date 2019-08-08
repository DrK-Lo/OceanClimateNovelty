Plot_interp<-function(data,column, B, title){ #Data is expected to be data.table format, column should be in the format 'data$column', title is a character string for the figure legend.
  FD12<-data[,c("long","lat")]
  FD12$sigma<-column
  FD12<-FD12[!is.na(FD12$sigma),]
  FD12[long>360,long:=long-360] #Correct format
  #FD12$long<-FD12$long-180 # Convert to WGS 1984 bounds
  FD12[long>180,long:=long-360]
  FD12_SP <- SpatialPoints(FD12) # Spatial points df
  proj4string(FD12_SP) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
  
  boolFALSE<-F
  MakeRast<-function(data){
    while(boolFALSE==F){
    tryCatch({
      gr <- as.data.frame(spsample(data, 'regular', n  = 10000000))
      names(gr) <- c('X', 'Y')
      coordinates(gr) <- c('X', 'Y')
      gridded(gr) <- TRUE 
      fullgrid(gr) <- TRUE
      proj4string(gr) <- proj4string(data)
      data.idw <- idw(sigma~1, data, newdata = gr, idp = 10)
      return(raster(data.idw))
      boolFalse<-T
      },
      error=function(e){
      },finally={})
      }
  }

  r<-MakeRast(FD12_SP)
  
  data('World', package = 'tmap')
  #get_projection(World)
  
  tm_shape(r) + 
    tm_raster(breaks = B, palette = 'plasma', # n = 10 may be better
              title = title) + 
    tm_legend(legend.outside = TRUE) +
    tm_shape(World, projection="longlat") +
    tm_fill()
}
