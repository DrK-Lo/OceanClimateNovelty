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
bbh <- fread("data/real_data/BBH.csv", sep = ",")

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

bbh$Years<-as.numeric(substr(bbh$COLLECTION_DATE,nchar(bbh$COLLECTION_DATE)-3,nchar(bbh$COLLECTION_DATE)))

summary(hots$Years)
summary(as.factor(hots$Years))
hist(hots$Years)
summary(bats$Years)
summary(as.factor(bats$Years))
hist(bats$Years)
summary(as.factor(bbh$Years))
hist(bbh$Years)

hots_sd<-hots[hots$`press(dbar)`<10 & hots$Years>=1990] %>%
group_by(Years) %>%
summarise(sd(`temp(ITS-90)`,na.rm=T),sd(ph,na.rm=T))
colnames(hots_sd)<-c("Years","Temp_sd","pH_sd")

ggplot(data=hots_sd)+
geom_point(aes(Years,Temp_sd))+
ylab("Surface Temperature (standard deviation)")+
ggtitle("HOTS")+
geom_hline(yintercept=mean(hots_sd$Temp_sd),col="red")

ggplot(data=hots_sd)+
geom_point(aes(Years,pH_sd))+
ylab("pH (standard deviation)")+
ggtitle("HOTS")+
geom_hline(yintercept=mean(hots_sd$pH_sd, na.rm = T),col="red")

bats_sd<-bats[bats$Depth < 10 & bats$Years>=1990]  %>%
group_by(Years) %>%
summarise(sd(Temp,na.rm=T))
colnames(bats_sd)<-c("Years","Temp_sd")

ggplot(data=bats_sd)+
geom_point(aes(Years,Temp_sd))+
ylab("Surface Temperature (standard deviation)")+
ggtitle("BATS")+
geom_hline(yintercept=mean(bats_sd$Temp_sd, na.rm = T),col="red")

bbh_sd<-bbh[bbh$Years>=1990]  %>%
group_by(Years) %>%
summarise(sd(`Sea Surface Temp Ave C`,na.rm=T))
colnames(bbh_sd)<-c("Years","Temp_sd")

ggplot(data=bbh_sd)+
geom_point(aes(Years,Temp_sd))+
ylab("Surface Temperature (standard deviation)")+
ggtitle("BBH")+
geom_hline(yintercept=mean(bbh_sd$Temp_sd, na.rm = T),col="red")

ggplot()+
geom_point(data=bats_sd,aes(bats_sd$Years,bats_sd$Temp_sd),col="blue")+
geom_point(data=hots_sd, aes(hots_sd$Years,hots_sd$Temp_sd),col="red")+
geom_point(data=bbh_sd, aes(bbh_sd$Years,bbh_sd$Temp_sd),col="cyan")+
xlab("Years")+
ylab("Surface Temperature (standard deviation)")+
geom_hline(yintercept=mean(bats_sd$Temp_sd, na.rm = T),col="blue")+
geom_hline(yintercept=mean(hots_sd$Temp_sd, na.rm = T),col="red")+
geom_hline(yintercept=mean(bbh_sd$Temp_sd, na.rm = T),col="cyan")

dat_2100_4.5_new <- fread("data/new_data/Katie_Temp_Arag_2070_2100_RCP45.txt", sep = ",")
colnames(dat_2100_4.5_new)<-c("No", "longitude", "latitude", "year", "month", "SST", "Aragonite")

dat_2100_8.5_new <- fread("data/new_data/Katie_Temp_Arag_2070_2100_RCP85.txt", sep = ",")
colnames(dat_2100_8.5_new)<-c("No", "longitude", "latitude", "year", "month", "SST", "Aragonite")

dat_1800_new <- fread("data/new_data/Katie_Temp_Arag_1800_2000.txt", sep = ",")
colnames(dat_1800_new)<-c("No", "longitude", "latitude", "year", "month", "SST", "Aragonite")


