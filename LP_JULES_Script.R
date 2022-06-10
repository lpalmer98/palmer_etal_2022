rm(list=ls())

setwd("~/Desktop/new_analysis/")

yr1<-1979
yr2<-2016


## 1- Upload data from each site

loc_site<-read.csv("~/Desktop/new_analysis/Site_location.csv",header=T)
loc_site<-as.data.frame(loc_site)
loc_site$Elev<-as.numeric(as.character(loc_site$Elev))

#loc_site <-loc_site[order(loc_site$Site, decreasing = FALSE), ]  

data_obs<-read.table("~/Desktop/new_analysis/data_compiled_by_Lewis.txt", header=T)

##data_obs<-subset(data_obs,!(Site == "Oxford"))

data_obs$Site<-factor(data_obs$Site, levels = loc_site$Site)
data_obs<-data_obs[order(data_obs$Site, decreasing = FALSE),] 




## 2- Open dataset
library(ncdf4)

JULES_data<-nc_open("~/Desktop/new_analysis/WFDEI_JULES.UK.month.2d.1979.2016-1.nc")
light<-nc_open("~/Desktop/new_analysis/PPFD_fapar_PAR.UK.month.2d.1982.2016.nc")

lon<-ncvar_get(JULES_data,'lon')
lat<-ncvar_get(JULES_data,'lat')

D13C_jules<-ncvar_get(JULES_data,'D13Cleaf')#[(3*12+1):456,,]
GPP_jules<-ncvar_get(JULES_data,'gpp')#[(3*12+1):456,,]
Tair_jules<-ncvar_get(JULES_data,'Tair')#[(3*12+1):456,,]
co2_jules<-ncvar_get(JULES_data,'co2')#[(3*12+1):456,,]
fsmc_jules<-ncvar_get(JULES_data,'fsmc')#[(3*12+1):456,,]
smc_avail_top_jules<-ncvar_get(JULES_data,'smc_avail_top')#[(3*12+1):456,,]
vpd_jules<-ncvar_get(JULES_data,'vpd')#[(3*12+1):456,,]
par_jules<-array(NA,dim = c(length(lon),length(lat),456))
par_jules[,,(3*12+1):456]<-ncvar_get(light,'PAR')


## 3a- Calculate GPP-weighted data for each variable for each grid-point to integrate the whole growing season
D13C_jules_gw<-array(NA,dim = c(length(lon),length(lat),length(yr1:yr2)))
GPP_jules_gw<-array(NA,dim = c(length(lon),length(lat),length(yr1:yr2)))
Tair_jules_gw<-array(NA,dim = c(length(lon),length(lat),length(yr1:yr2)))
co2_jules_gw<-array(NA,dim = c(length(lon),length(lat),length(yr1:yr2)))
fsmc_jules_gw<-array(NA,dim = c(length(lon),length(lat),length(yr1:yr2)))
smc_avail_top_jules_gw<-array(NA,dim = c(length(lon),length(lat),length(yr1:yr2)))
vpd_jules_gw<-array(NA,dim = c(length(lon),length(lat),length(yr1:yr2)))
par_jules_gw<-array(NA,dim = c(length(lon),length(lat),length(yr1:yr2)))

for (z in 1:length(lat)){
  for (k in 1:length(lon)){
    i<-1
    for (j in 1:length(yr1:yr2)){
      D13C_jules_gw[k,z,j]<-sum(D13C_jules[(i:(i+11)),k,z]*GPP_jules[(i:(i+11)),k,z],na.rm=T)/sum(GPP_jules[(i:(i+11)),k,z],na.rm=T)
      GPP_jules_gw[k,z,j]<-mean(GPP_jules[(i:(i+11)),k,z],na.rm=T)
      Tair_jules_gw[k,z,j]<-sum(Tair_jules[(i:(i+11)),k,z]*GPP_jules[(i:(i+11)),k,z],na.rm=T)/sum(GPP_jules[(i:(i+11)),k,z],na.rm=T)
      co2_jules_gw[k,z,j]<-sum(co2_jules[(i:(i+11)),k,z]*GPP_jules[(i:(i+11)),k,z],na.rm=T)/sum(GPP_jules[(i:(i+11)),k,z],na.rm=T)
      fsmc_jules_gw[k,z,j]<-sum(fsmc_jules[(i:(i+11)),k,z]*GPP_jules[(i:(i+11)),k,z],na.rm=T)/sum(GPP_jules[(i:(i+11)),k,z],na.rm=T)
      smc_avail_top_jules_gw[k,z,j]<-sum(smc_avail_top_jules[(i:(i+11)),k,z]*GPP_jules[(i:(i+11)),k,z],na.rm=T)/sum(GPP_jules[(i:(i+11)),k,z],na.rm=T)
      vpd_jules_gw[k,z,j]<-sum(vpd_jules[(i:(i+11)),k,z]*GPP_jules[(i:(i+11)),k,z],na.rm=T)/sum(GPP_jules[(i:(i+11)),k,z],na.rm=T)
     par_jules_gw[k,z,j]<-sum(par_jules[k,z,(i:(i+11))]*GPP_jules[(i:(i+11)),k,z],na.rm=T)/sum(GPP_jules[(i:(i+11)),k,z],na.rm=T)
      
      i<-i+12
      
    }
  }
}


## 3b- Average data only for July-August
mth1<-7
mth2<-8

D13C_jules_ja<-array(NA,dim = c(length(lon),length(lat),length(yr1:yr2)))
GPP_jules_ja<-array(NA,dim = c(length(lon),length(lat),length(yr1:yr2)))
Tair_jules_ja<-array(NA,dim = c(length(lon),length(lat),length(yr1:yr2)))
co2_jules_ja<-array(NA,dim = c(length(lon),length(lat),length(yr1:yr2)))
fsmc_jules_ja<-array(NA,dim = c(length(lon),length(lat),length(yr1:yr2)))
smc_avail_top_jules_ja<-array(NA,dim = c(length(lon),length(lat),length(yr1:yr2)))
vpd_jules_ja<-array(NA,dim = c(length(lon),length(lat),length(yr1:yr2)))
par_jules_ja<-array(NA,dim = c(length(lon),length(lat),length(yr1:yr2)))

for (z in 1:length(lat)){
  for (k in 1:length(lon)){
    i<-1
    
    for (j in 1:length(yr1:yr2)){
      
      D13C_jules_ja[k,z,j]<-mean(D13C_jules[(i:(i+11)),k,z][mth1:mth2],na.rm=T)
      GPP_jules_ja[k,z,j]<-mean(GPP_jules[(i:(i+11)),k,z][mth1:mth2],na.rm=T)
      Tair_jules_ja[k,z,j]<-mean(Tair_jules[(i:(i+11)),k,z][mth1:mth2],na.rm=T)
      co2_jules_ja[k,z,j]<-mean(co2_jules[(i:(i+11)),k,z][mth1:mth2],na.rm=T)
      fsmc_jules_ja[k,z,j]<-mean(fsmc_jules[(i:(i+11)),k,z][mth1:mth2],na.rm=T)
      smc_avail_top_jules_ja[k,z,j]<-mean(smc_avail_top_jules[(i:(i+11)),k,z][mth1:mth2],na.rm=T)
      vpd_jules_ja[k,z,j]<-mean(vpd_jules[(i:(i+11)),k,z][mth1:mth2],na.rm=T)
      par_jules_ja[k,z,j]<-mean(par_jules[k,z,(i:(i+11))][mth1:mth2],na.rm=T)
      i<-i+12
      
    }
    
  }
}



## 4- Extract your datapoints for the different sites

for (i in 1:nrow(loc_site)){
ilat=which.min(abs(lat - loc_site$Lat[i]))
ilon=which.min(abs(lon - loc_site$Lon[i]))
  
## Extract 3 closest latitudes and longitudes gridpoints
ilat2=c((which.min(abs(lat - loc_site$Lat[i]))-1):(which.min(abs(lat - loc_site$Lat[i]))+1))
ilon2=c((which.min(abs(lon - loc_site$Lon[i]))-1):(which.min(abs(lon - loc_site$Lon[i]))+1))

if(i==1){
 data_gw <- cbind(as.character(loc_site[i,1]),loc_site[i,4]/1000,c(1979:2016),D13C_jules_gw[ilon,ilat,],Tair_jules_gw[ilon,ilat,],co2_jules_gw[ilon,ilat,],fsmc_jules_gw[ilon,ilat,],smc_avail_top_jules_gw[ilon,ilat,],vpd_jules_gw[ilon,ilat,]/1000,par_jules_gw[ilon,ilat,])
 data_ja <- cbind(as.character(loc_site[i,1]),loc_site[i,4]/1000,c(1979:2016),D13C_jules_ja[ilon,ilat,],Tair_jules_ja[ilon,ilat,],co2_jules_ja[ilon,ilat,],fsmc_jules_ja[ilon,ilat,],smc_avail_top_jules_ja[ilon,ilat,],vpd_jules_ja[ilon,ilat,]/1000,par_jules_ja[ilon,ilat,])

  data_gw2 <- cbind(as.character(loc_site[i,1]),loc_site[i,4]/1000,rep(1979:2016,each=3),D13C_jules_gw[ilon2,ilat2,],Tair_jules_gw[ilon2,ilat2,],co2_jules_gw[ilon2,ilat2,],fsmc_jules_gw[ilon2,ilat2,],smc_avail_top_jules_gw[ilon2,ilat2,],vpd_jules_gw[ilon2,ilat2,]/1000,par_jules_gw[ilon2,ilat2,])
 data_ja2 <- cbind(as.character(loc_site[i,1]),loc_site[i,4]/1000,rep(1979:2016,each=3),D13C_jules_ja[ilon2,ilat2,],Tair_jules_ja[ilon2,ilat2,],co2_jules_ja[ilon2,ilat2,],fsmc_jules_ja[ilon2,ilat2,],smc_avail_top_jules_ja[ilon2,ilat2,],vpd_jules_ja[ilon2,ilat2,]/1000,par_jules_ja[ilon2,ilat2,])
}else{
  data_gw <-rbind(data_gw,cbind(as.character(loc_site[i,1]),loc_site[i,4]/1000,c(1979:2016),D13C_jules_gw[ilon,ilat,],Tair_jules_gw[ilon,ilat,],co2_jules_gw[ilon,ilat,],fsmc_jules_gw[ilon,ilat,],smc_avail_top_jules_gw[ilon,ilat,],vpd_jules_gw[ilon,ilat,]/1000,par_jules_gw[ilon,ilat,]))
  data_ja <-rbind(data_ja,cbind(as.character(loc_site[i,1]),loc_site[i,4]/1000,c(1979:2016),D13C_jules_ja[ilon,ilat,],Tair_jules_ja[ilon,ilat,],co2_jules_ja[ilon,ilat,],fsmc_jules_ja[ilon,ilat,],smc_avail_top_jules_ja[ilon,ilat,],vpd_jules_ja[ilon,ilat,]/1000,par_jules_ja[ilon,ilat,]))
  
  data_gw2 <-rbind(data_gw2,cbind(as.character(loc_site[i,1]),loc_site[i,4]/1000,rep(1979:2016,each=3),D13C_jules_gw[ilon2,ilat2,],Tair_jules_gw[ilon2,ilat2,],co2_jules_gw[ilon2,ilat2,],fsmc_jules_gw[ilon2,ilat2,],smc_avail_top_jules_gw[ilon2,ilat2,],vpd_jules_gw[ilon2,ilat2,]/1000,par_jules_gw[ilon2,ilat2,]))
  data_ja2 <-rbind(data_ja2,cbind(as.character(loc_site[i,1]),loc_site[i,4]/1000,rep(1979:2016,each=3),D13C_jules_ja[ilon2,ilat2,],Tair_jules_ja[ilon2,ilat2,],co2_jules_ja[ilon2,ilat2,],fsmc_jules_ja[ilon2,ilat2,],smc_avail_top_jules_ja[ilon2,ilat2,],vpd_jules_ja[ilon2,ilat2,]/1000,par_jules_ja[ilon2,ilat2,]))
}

}

## Dataframe with growing season GPP-weighted variables
data_gw<-as.data.frame(data_gw)
colnames(data_gw) <-c("Site","Elev_km","Year","D13C_permil","Tair_degC","CO2_ppm","fsmc_unitless","SMC_mm","VPD_kPa","par_molm2mth")

data_gw$Year<-as.numeric(as.character(data_gw$Year))
data_gw$Elev_km<-as.numeric(as.character(data_gw$Elev_km))
data_gw$D13C_permil<-as.numeric(as.character(data_gw$D13C_permil))
data_gw$CO2_ppm<-as.numeric(as.character(data_gw$CO2_ppm))
data_gw$Tair_degC<-as.numeric(as.character(data_gw$Tair_degC))
data_gw$fsmc_unitless<-as.numeric(as.character(data_gw$fsmc_unitless))
data_gw$SMC_mm<-as.numeric(as.character(data_gw$SMC_mm))
data_gw$VPD_kPa<-as.numeric(as.character(data_gw$VPD_kPa))
data_gw$par_molm2mth<-as.numeric(as.character(data_gw$par_molm2mth))
data_gw$par_molm2mth[data_gw$par_molm2mth == 0] <- NA
head(data_gw)


## Dataframe with July-August average variables
data_ja<-as.data.frame(data_ja)
colnames(data_ja) <-c("Site","Elev_km","Year","D13C_permil","Tair_degC","CO2_ppm","fsmc_unitless","SMC_mm","VPD_kPa","par_molm2mth")

data_ja$Year<-as.numeric(as.character(data_ja$Year))
data_ja$Elev_km<-as.numeric(as.character(data_ja$Elev_km))
data_ja$D13C_permil<-as.numeric(as.character(data_ja$D13C_permil))
data_ja$CO2_ppm<-as.numeric(as.character(data_ja$CO2_ppm))
data_ja$Tair_degC<-as.numeric(as.character(data_ja$Tair_degC))
data_ja$fsmc_unitless<-as.numeric(as.character(data_ja$fsmc_unitless))
data_ja$SMC_mm<-as.numeric(as.character(data_ja$SMC_mm))
data_ja$VPD_kPa<-as.numeric(as.character(data_ja$VPD_kPa))
data_ja$par_molm2mth<-as.numeric(as.character(data_ja$par_molm2mth))
data_ja$par_molm2mth[data_ja$par_molm2mth == 0] <- NA

head(data_ja)


## Dataframe with growing season GPP-weighted variables
data_gw2<-as.data.frame(data_gw2)
colnames(data_gw2) <-c("Site","Elev_km","Year","D13C_permil","Tair_degC","CO2_ppm","fsmc_unitless","SMC_mm","VPD_kPa","par_molm2mth")

data_gw2$Year<-as.numeric(as.character(data_gw2$Year))
data_gw2$Elev_km<-as.numeric(as.character(data_gw2$Elev_km))
data_gw2$D13C_permil<-as.numeric(as.character(data_gw2$D13C_permil))
data_gw2$CO2_ppm<-as.numeric(as.character(data_gw2$CO2_ppm))
data_gw2$Tair_degC<-as.numeric(as.character(data_gw2$Tair_degC))
data_gw2$fsmc_unitless<-as.numeric(as.character(data_gw2$fsmc_unitless))
data_gw2$SMC_mm<-as.numeric(as.character(data_gw2$SMC_mm))
data_gw2$VPD_kPa<-as.numeric(as.character(data_gw2$VPD_kPa))
data_gw2$par_molm2mth<-as.numeric(as.character(data_gw2$par_molm2mth))
data_gw2$par_molm2mth[data_gw2$par_molm2mth == 0] <- NA



## Dataframe with July-August average variables
data_ja2<-as.data.frame(data_ja2)
colnames(data_ja2) <-c("Site","Elev_km","Year","D13C_permil","Tair_degC","CO2_ppm","fsmc_unitless","SMC_mm","VPD_kPa","par_molm2mth")

data_ja2$Year<-as.numeric(as.character(data_ja2$Year))
data_ja2$Elev_km<-as.numeric(as.character(data_ja2$Elev_km))
data_ja2$D13C_permil<-as.numeric(as.character(data_ja2$D13C_permil))
data_ja2$CO2_ppm<-as.numeric(as.character(data_ja2$CO2_ppm))
data_ja2$Tair_degC<-as.numeric(as.character(data_ja2$Tair_degC))
data_ja2$fsmc_unitless<-as.numeric(as.character(data_ja2$fsmc_unitless))
data_ja2$SMC_mm<-as.numeric(as.character(data_ja2$SMC_mm))
data_ja2$VPD_kPa<-as.numeric(as.character(data_ja2$VPD_kPa))
data_ja2$par_molm2mth<-as.numeric(as.character(data_ja2$par_molm2mth))
data_ja2$par_molm2mth[data_ja2$par_molm2mth == 0] <- NA



## Calculate average of grid-points for each year and site
  
for (i in unique(data_gw2$Site)){  
  for (j in yr1:yr2){
    if(i==unique(data_gw2$Site)[1] & j==yr1){
  data_gw_final <- c(i,apply(subset(data_gw2,Site ==i & Year == j)[,2:10],2,mean,na.rm=T))
  data_ja_final <- c(i,apply(subset(data_ja2,Site ==i & Year == j)[,2:10],2,mean,na.rm=T))
    }else{
  data_gw_final <- rbind(data_gw_final,c(i,apply(subset(data_gw2,Site ==i & Year == j)[,2:10],2,mean,na.rm=T)))
  data_ja_final <- rbind(data_ja_final,c(i,apply(subset(data_ja2,Site ==i & Year == j)[,2:10],2,mean,na.rm=T)))
  }  }
  }


## Dataframe with growing season GPP-weighted variables
data_gw_final<-as.data.frame(data_gw_final)
colnames(data_gw_final) <-c("Site","Elev_km","Year","D13C_permil","Tair_degC","CO2_ppm","fsmc_unitless","SMC_mm","VPD_kPa","par_molm2mth")

data_gw_final$Year<-as.numeric(as.character(data_gw_final$Year))
data_gw_final$Elev_km<-as.numeric(as.character(data_gw_final$Elev_km))
data_gw_final$D13C_permil<-as.numeric(as.character(data_gw_final$D13C_permil))
data_gw_final$CO2_ppm<-as.numeric(as.character(data_gw_final$CO2_ppm))
data_gw_final$Tair_degC<-as.numeric(as.character(data_gw_final$Tair_degC))
data_gw_final$fsmc_unitless<-as.numeric(as.character(data_gw_final$fsmc_unitless))
data_gw_final$SMC_mm<-as.numeric(as.character(data_gw_final$SMC_mm))
data_gw_final$VPD_kPa<-as.numeric(as.character(data_gw_final$VPD_kPa))
data_gw_final$par_molm2mth<-as.numeric(as.character(data_gw_final$par_molm2mth))
head(data_gw_final)


## Dataframe with growing season JA variables
data_ja_final<-as.data.frame(data_ja_final)
colnames(data_ja_final) <-c("Site","Elev_km","Year","D13C_permil","Tair_degC","CO2_ppm","fsmc_unitless","SMC_mm","VPD_kPa","par_molm2mth")

data_ja_final$Year<-as.numeric(as.character(data_ja_final$Year))
data_ja_final$Elev_km<-as.numeric(as.character(data_ja_final$Elev_km))
data_ja_final$D13C_permil<-as.numeric(as.character(data_ja_final$D13C_permil))
data_ja_final$CO2_ppm<-as.numeric(as.character(data_ja_final$CO2_ppm))
data_ja_final$Tair_degC<-as.numeric(as.character(data_ja_final$Tair_degC))
data_ja_final$fsmc_unitless<-as.numeric(as.character(data_ja_final$fsmc_unitless))
data_ja_final$SMC_mm<-as.numeric(as.character(data_ja_final$SMC_mm))
data_ja_final$VPD_kPa<-as.numeric(as.character(data_ja_final$VPD_kPa))
data_ja_final$par_molm2mth<-as.numeric(as.character(data_ja_final$par_molm2mth))
head(data_ja_final)




head(data_obs)

# Version when extracting the exact lat/lon
data_gw$Obs <- data_obs$Obs
data_ja$Obs <- data_obs$Obs

# Version when averaging over several lat/lon gridpoints
data_gw_final$Obs <- data_obs$Obs
data_ja_final$Obs <- data_obs$Obs


## Agreement between observations and predictions

## Weighted with GPP
summary(lm(data_gw$D13C_permil~data_gw$Obs))   ## R2 = 0.16
summary(lm(data_gw_final$D13C_permil~data_gw_final$Obs))  ## R2 = 0.19


## For July-August
summary(lm(data_ja$D13C_permil~data_ja$Obs))  ## R2 = 0.22
summary(lm(data_ja_final$D13C_permil~data_ja_final$Obs))  ## R2 = 0.22



### 5- Multiple regression models over 1982-2016

library(performance)
library(car)
library(MASS)

## selection data
data = data_ja  ## average for July-August

## a- Simple model for D13C values againt JA data
composite = subset(data_ja,!is.na(par_molm2mth) & !is.na(Obs))

mod1a <-lm(data=data_ja, Obs ~ Tair_degC + CO2_ppm + VPD_kPa + Elev_km)
summary(mod1a)  ## R2 = 0.22 (observations)

mod1b <-lm(data=data_ja, D13C_permil ~ Tair_degC + CO2_ppm + VPD_kPa + Elev_km)
summary(mod1b)  ## R2 = 0.78 (predictions)

step <- stepAIC(mod1a, direction="both")
step$anova  
## give you an indication of the model more suitable
## here final model include: Tair_degC + VPD_kPa

## So we consider only these environmental variables for comparison with predictions
mod1a <-lm(data=data, Obs ~ Tair_degC + VPD_kPa)
summary(mod1a)  ## R2 = 0.56 (observations)

mod2a <-lm(data=data, D13C_permil ~ Tair_degC + VPD_kPa)
summary(mod2a)  ## R2 = 0.57 (predictions)

compare_performance(mod1a,mod2a) ## RMSE values


## b- Individual sites

i=1;loc_site$Site[i]  ## change here number of site

site = subset(data,Site == loc_site$Site[i] & !is.na(par_molm2mth) & !is.na(Obs))
mod1a <-lm(data=site, Obs ~ Tair_degC + CO2_ppm + VPD_kPa + Elev_km)

step <- stepAIC(mod1a, direction="both")
step$anova 

as.character(step$call)[2]  ## model to use for both observations and predictions

## Manually change environmental drivers below and report them in Tables 4/5:
mod1a <-lm(data=site, Obs ~  Tair_degC + CO2_ppm + VPD_kPa) 
mod2a <-lm(data=site, D13C_permil ~  Tair_degC + CO2_ppm + VPD_kPa)

summary(mod1a)  
summary(mod2a)  
compare_performance(mod1a,mod2a)  ## extract RMSE value


## c- Figure with regression coefficient

library(broom)
library(GGally)
require(broom.helpers)
library(sjPlot)
library(ggpubr)
library(ggplot2)

## Reorder from high to low latitudes (north to south)
## loc_site <- loc_site[order(loc_site$Lat,decreasing=TRUE),]

model_obs = "Obs ~ Tair_degC + CO2_ppm + VPD_kPa"
model_jules = "D13C_permil ~ Tair_degC + CO2_ppm + VPD_kPa"

site1 <- glm(model_obs, data = subset(data,Site == loc_site$Site[1] & !is.na(Obs)))
site2 <- glm(model_obs, data = subset(data,Site == loc_site$Site[2] & !is.na(Obs)))
site3 <- glm(model_obs, data = subset(data,Site == loc_site$Site[3] & !is.na(Obs)))
site4 <- glm(model_obs, data = subset(data,Site == loc_site$Site[4] & !is.na(Obs)))
site5 <- glm(model_obs, data = subset(data,Site == loc_site$Site[5] & !is.na(Obs)))
site6 <- glm(model_obs, data = subset(data,Site == loc_site$Site[6] & !is.na(Obs)))
site7 <- glm(model_obs, data = subset(data,Site == loc_site$Site[7] & !is.na(Obs)))
site8 <- glm(model_obs, data = subset(data,Site == loc_site$Site[8] & !is.na(Obs)))
site9 <- glm(model_obs, data = subset(data,Site == loc_site$Site[9] & !is.na(Obs)))
site10 <- glm(model_obs, data = subset(data,Site == loc_site$Site[10] & !is.na(Obs)))
site11 <- glm(model_obs, data = subset(data,Site == loc_site$Site[11] & !is.na(Obs)))
site12 <- glm(model_obs, data = subset(data,Site == loc_site$Site[12] & !is.na(Obs)))


site1_jules <- glm(model_jules, data = subset(data,Site == loc_site$Site[1] & !is.na(Obs)))
site2_jules <- glm(model_jules, data = subset(data,Site == loc_site$Site[2] & !is.na(Obs)))
site3_jules <- glm(model_jules, data = subset(data,Site == loc_site$Site[3] & !is.na(Obs)))
site4_jules <- glm(model_jules, data = subset(data,Site == loc_site$Site[4] & !is.na(Obs)))
site5_jules <- glm(model_jules, data = subset(data,Site == loc_site$Site[5] & !is.na(Obs)))
site6_jules <- glm(model_jules, data = subset(data,Site == loc_site$Site[6] & !is.na(Obs)))
site7_jules <- glm(model_jules, data = subset(data,Site == loc_site$Site[7] & !is.na(Obs)))
site8_jules <- glm(model_jules, data = subset(data,Site == loc_site$Site[8] & !is.na(Obs)))
site9_jules <- glm(model_jules, data = subset(data,Site == loc_site$Site[9] & !is.na(Obs)))
site10_jules <- glm(model_jules, data = subset(data,Site == loc_site$Site[10] & !is.na(Obs)))
site11_jules <- glm(model_jules, data = subset(data,Site == loc_site$Site[11] & !is.na(Obs)))
site12_jules <- glm(model_jules, data = subset(data,Site == loc_site$Site[12] & !is.na(Obs)))


p1 <- plot_models(site1,site2,site3, site4,site5,site6,site7,site8,site9,site10,site11,site12,show.values = TRUE,colors = "Paired",show.intercept = FALSE,
            value.size = 3,show.p = TRUE, #,p.shape = TRUE,
            dot.size = 2,axis.lim=c(-15,40),
            line.size = 1.0,spacing = 0.85,
            legend.title = "Sites",std.est = TRUE,m.labels = loc_site$Site,title = " ",axis.labels = c("VPD (kPa)","CO2 (ppm)","Tair (°C)")) 
            
p1$layers[[4]]$aes_params$vjust<-(+0.2)
p1$layers[[4]]$aes_params$hjust<-(-1.8) 
           
p1 + scale_color_brewer(palette = "Paired", labels = loc_site$Site[12:1])
            




p2 <- plot_models(site1_jules,site2_jules,site3_jules, site4_jules,site5_jules,site6_jules,site7_jules,site8_jules,site9_jules,site10_jules,site11_jules,site12_jules,show.values = TRUE,colors = "Paired",show.intercept = FALSE,
                  value.size = 3,show.p = TRUE, #,p.shape = TRUE,
                  dot.size = 2,axis.lim=c(-15,40),
                  line.size = 1.0,spacing = 0.85,
                  legend.title = "Sites",std.est = TRUE, m.labels = loc_site$Site, title = " ",axis.labels = c("VPD (kPa)","CO2 (ppm)","Tair (°C)"))
         

p2$layers[[4]]$aes_params$vjust<-(+0.2)
p2$layers[[4]]$aes_params$hjust<-(-1.8)

p2 + scale_color_brewer(palette = "Paired", labels = loc_site$Site[12:1])  

combo_plot_D13C<-ggarrange(p1 + scale_color_brewer(palette = "Paired", labels = loc_site$Site[12:1]),p2 + scale_color_brewer(palette = "Paired", labels = loc_site$Site[12:1]),ncol = 2,nrow=1, align="v",widths = c(2,2), vjust = 1,hjust = 0.05,labels = c("(a) Tree rings","(b) JULES"),label.x = 0.035,label.y = 0.99,font.label = list(size = 10),common.legend = TRUE,legend = "right")
combo_plot_D13C

ggsave(paste0("Figure_reg_coef.jpeg"),combo_plot_D13C, width = 10,height=8, dpi = 1000)



## Figure 2
library(ggplot2)
library(ggpubr)
library(ggpmisc)
 
data_corr <-data_ja[,c("Site","Obs","D13C_permil")]
data_corr$Obs <- data_corr$Obs + 1960
 
plot1<-ggplot(data_ja) +
  theme_classic() +
  scale_x_continuous(limits = c(1979,2016),breaks = seq(1980,2016,10)) +
  scale_y_continuous(limits = c(18,24),breaks = seq(18,24,2)) +
  labs(x=" ",y= expression(Delta^{13}*"C (‰)")) +
  geom_line(aes(x=Year,y=Obs,group = Site,color = "Observations"),size=1,alpha=0.5,show.legend = TRUE) +
  geom_line(aes(x=Year,y=D13C_permil,group = Site,color = "Predictions"),size=1,alpha=0.5,lty=1,show.legend = TRUE) +
  geom_point(data=data_corr,aes(x=Obs,y=D13C_permil,group = Site),col="white",alpha=0.01) +
  stat_cor(data=data_corr,aes(x=Obs,y=D13C_permil,group = Site),method = "spearman",label.x = 1980, label.y = 23.5,digits = 2,r.digits =2,p.digits = 2,p.accuracy=0.001) +
  facet_wrap(~Site,ncol=4) +
  geom_smooth(aes(x=Year,y=Obs,group = Site,color = "Observations"),method = 'lm',formula=y ~ x,cex=0.75,alpha=0.2,show.legend = TRUE) +
  geom_smooth(aes(x=Year,y=D13C_permil,group = Site,color = "Predictions"),method = 'lm',formula=y ~ x,cex=0.75,alpha=0.2,show.legend = TRUE) +
  scale_color_manual(values=c("black","blue"),labels=c("Observations","Predictions")) +
  theme(
    axis.ticks.length = unit(.15, "cm"),
    axis.text.y = element_text(size=12),
    axis.text.x = element_text(size=12),
    axis.title.x = element_text( size=12),
    axis.title.y = element_text(size=12),
    legend.text = element_text(family="sans", size = 12),
    legend.background = element_rect(fill=NA),
    legend.title = element_blank(),
    legend.position = "bottom",
    text=element_text(family="Helvetica"))
plot1
 
ggsave("Figure2.jpeg",plot1, width = 10,height=8, dpi = 1000)


## Figure S1 - 3 gridsquares
library(ggplot2)
library(ggpubr)
library(ggpmisc)
 
data_corr <-data_ja_final[,c("Site","Obs","D13C_permil")]
data_corr$Obs <- data_corr$Obs + 1960
 
plot2<-ggplot(data_ja_final) +
  theme_classic() +
  scale_x_continuous(limits = c(1979,2016),breaks = seq(1980,2016,10)) +
  scale_y_continuous(limits = c(18,24),breaks = seq(18,24,2)) +
  labs(x=" ",y= expression(Delta^{13}*"C (‰)")) +
  geom_line(aes(x=Year,y=Obs,group = Site,color = "Observations"),size=1,alpha=0.5,show.legend = TRUE) +
  geom_line(aes(x=Year,y=D13C_permil,group = Site,color = "Predictions"),size=1,alpha=0.5,lty=1,show.legend = TRUE) +
  geom_point(data=data_corr,aes(x=Obs,y=D13C_permil,group = Site),col="white",alpha=0.01) +
  stat_cor(data=data_corr,aes(x=Obs,y=D13C_permil,group = Site),method = "spearman",label.x = 1980, label.y = 23.5,digits = 2,r.digits =2,p.digits = 2,p.accuracy=0.001) +
  facet_wrap(~Site,ncol=4) +
  geom_smooth(aes(x=Year,y=Obs,group = Site,color = "Observations"),method = 'lm',formula=y ~ x,cex=0.75,alpha=0.2,show.legend = TRUE) +
  geom_smooth(aes(x=Year,y=D13C_permil,group = Site,color = "Predictions"),method = 'lm',formula=y ~ x,cex=0.75,alpha=0.2,show.legend = TRUE) +
  scale_color_manual(values=c("black","blue"),labels=c("Observations","Predictions")) +
  theme(
    axis.ticks.length = unit(.15, "cm"),
    axis.text.y = element_text(size=12),
    axis.text.x = element_text(size=12),
    axis.title.x = element_text( size=12),
    axis.title.y = element_text(size=12),
    legend.text = element_text(family="sans", size = 12),
    legend.background = element_rect(fill=NA),
    legend.title = element_blank(),
    legend.position = "bottom",
    text=element_text(family="Helvetica"))
plot2
 
ggsave("FigureS1.jpeg",plot2, width = 10,height=8, dpi = 1000)


