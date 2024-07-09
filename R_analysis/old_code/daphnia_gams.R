#==============================================================================
# Model species distributions based on abundance, using a GAMM framework.
#==============================================================================
library(tidyverse)
library(mgcv)
#==============================================================================
#Load data sets:
#==============================================================================
#Current environment from worldclim
currentEnv=raster::getData('worldclim', var="bio", res=2.5, download=T)

#Future environment based on the 8.5 concentration scenario for 2070 from HADGEM2-ES model
futureEnv=raster::getData('CMIP5', var='bio', res=2.5, rcp=85, model='HE', year=70, download=T)
names(futureEnv)=names(currentEnv)

###This set of data comes from Jones and 
#Contemporary counts (2011)
zabundC<-read.csv("contempA.csv")
colnames(zabundC)

#Historical counts
zabundH<-read.csv("historicalA.csv")
colnames(zabundH)

#Environmental data set
zenv<-read.csv("vars.csv")
colnames(zenv)

names(zabundC)[names(zabundC)=="lake"] <- "Lake.name"
colnames(zabundC)
zabundC = data.frame(year = c(matrix(2011,nrow(zabundC),1) ), zabundC )

names(zabundH)[names(zabundH)=="lake"] <- "Lake.name"
colnames(zabundH)

names(zenv)[names(zenv)=="Longitude..Â.."]<-"Longitude"
names(zenv)[names(zenv)=="Latitude.Â.."]<-"Latitude"
colnames(zenv)

#Merge the two Zooplankton sets
zabund = full_join (zabundC,zabundH)

abundance_environment <- merge(zabundC, zenv, by=c('Lake.name'))
head(abundance_environment)
view(abundance_environment)

#z1_gam = gam ( D.pulex ~ te(Latitude,Longitude)+ s(Lake.size..km2.) +  s(pH)+ s (year, bs = "re", k=3), data =abundance_environment  )
z1_gam = gam ( D.pulex ~ te(Latitude,Longitude)+ s(pH)+s (year, bs = "re"), data =abundance_environment  )
summary(z1_gam)
AIC(z1_gam)
plot(z1_gam)
coef(z1_gam)
evaluate(abundance_environment, abundance_environment, z1_gam)
?evaluate
?maxent
?gam
