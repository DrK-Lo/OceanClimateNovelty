## 2019-07-09
## Áki Jarl - Ocean Climate Novelty script to view variance at existing oceanographic meaurement opperations
##
setwd("~/Desktop/PostDoc/OceanClimateNovelty/")
source("src/0-Novelty_Oceans_Functions.R")

require('bit64')
hots <- fread("data/real_data/HOT_bottle_data_mod.txt", sep = ",")
hots <- hots[hots$`date(mmddyy)` != -9] 
hots[hots == -9] <- NA
bats <- fread("data/real_data/BATS_bottle_data_mod.txt", sep = "\t")
bats[bats == -999] <- NA

years<-NULL
  for(i in as.numeric(substr(hots$`date(mmddyy)`,nchar(hots$`date(mmddyy)`)-1,nchar(hots$`date(mmddyy)`)))){
  if(i>80){
    years<-c(years,as.numeric(paste(19,i,sep="")))
  }
  else if(nchar(i)==2){
    years<-c(years,as.numeric(paste(20,i,sep="")))
  }
  else{
    years<-c(years,as.numeric(paste(200,i,sep="")))
  }
}

hots$Years<-years
bats$Years<-as.numeric(substr(bats$yyyymmdd,1,4))

summary(hots$Years)
hist(hots$Years)
summary(bats$Years)
hist(bats$Years)


hots_sd<-hots[hots$`press(dbar)`<10] %>% 
  group_by(Years) %>% 
  summarise(sd(`temp(ITS-90)`,na.rm=T))
colnames(hots_sd)<-c("Years","Temp_sd")

bats_sd<-bats[bats$Depth < 10]  %>% 
  group_by(Years) %>% 
  summarise(sd(Temp,na.rm=T))  
colnames(bats_sd)<-c("Years","Temp_sd")
        
        