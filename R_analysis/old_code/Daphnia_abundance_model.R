library(tidyverse)
library(raster)
library(rgdal)
library(maps)
library(mapdata)
library(dismo)  
library(rJava)  
library(maptools)
library(jsonlite)
library(rgbif) 
library("readxl")
install.packages("xlsx")
library("xlsx")

#Current environment from worldclim
currentEnv=getData('worldclim', var="bio", res=2.5, download=T)
view(currentEnv)
#Future environment based on the 8.5 concentration scenario for 2070 from HADGEM2-ES model
futureEnv=getData('CMIP5', var='bio', res=2.5, rcp=85, model='HE', year=70, download=T)
names(futureEnv)=names(currentEnv)

#World map:
data(wrld_simpl)

dir<-'C:\\Users\\larac\\Desktop\\BMSC DS - GitHub\\BMSC-DS'
setwd(dir)

zooplankton_abundance_data<-read.csv("Zooplankton abundance data - historical and contemporary.csv")
colnames(zooplankton_abundance_data)

zooplankton_environmental_data<-read.csv("Zooplankton abundance data - envrionmental data.csv")
colnames(zooplankton_environmental_data)

names(zooplankton_abundance_data)[names(zooplankton_abundance_data)=="lake"] <- "Lake.name"
colnames(zooplankton_abundance_data)
names(zooplankton_environmental_data)[names(zooplankton_environmental_data)=="Longitude...."]<-"Longitude"
names(zooplankton_environmental_data)[names(zooplankton_environmental_data)=="Latitude..."]<-"Latitude"
colnames(zooplankton_environmental_data)

abundance_environment <- merge(zooplankton_abundance_data, zooplankton_environmental_data, by=c('Lake.name'))
head(abundance_environment)
view(abundance_environment)

daphnia_pulex_abundance<-filter(abundance_environment, period=="contemporary") %>% 
  select(Lake.name, year, period, date, abundance, D.pulex, Latitude:pH) %>% 
  filter(D.pulex > 0)
daphnia_longiremis_abundance<-filter(abundance_environment, period=="contemporary") %>% 
  select(Lake.name, year, period, date, abundance, D.longiremis, Latitude:pH) %>% 
  filter(D.longiremis > 0)

#Initial plot
plot(wrld_simpl, xlim=c(min(daphnia_pulex_abundance$Longitude)-1,max(daphnia_pulex_abundance$Longitude)+1), ylim=c(min(daphnia_abundance$Latitude)-1,max(daphnia_abundance$Latitude)+1), axes=TRUE, col="light yellow")
points(daphnia_pulex_abundance$Longitude, daphnia_pulex_abundance$Latitude, col="orange", pch=20, cex=0.75)
points(daphnia_longiremis_abundance$Longitude, daphnia_longiremis_abundance$Latitude, col="red", pch=20, cex=0.75)


########Not sure how to adapt MaxEnt to accound for abundance data...
library(lubridate)
library(sf)
library(raster)
install.packages("dggridR")
library(dggridR)
install.packages("pdp")
library(pdp)
library(mgcv)
install.packages("fitdistrplus")
library(fitdistrplus)
library(viridis)
install.packages("fields")
library(fields)
library(tidyverse)
library(ggplot2)
install.packages(sim_sad)
# resolve namespace conflicts
select <- dplyr::select
map <- purrr::map
projection <- raster::projection

sad_lnorm1 <- sim_sad(s_pool = 100, n_sim = 10000, sad_type = "lnorm",
                      sad_coef = list("meanlog" = 5, "sdlog" = 0.5))
plot(sad_lnorm1, method = "octave")
plot(sad_lnorm1, method = "rank")



model.extent<-extent(min(daphnia_abundance$Longitude)-10,max(daphnia_abundance$Longitude)+10,min(daphnia_abundance$Latitude)-10,max(daphnia_abundance$Latitude)+10)
modelEnv=crop(currentEnv, model.extent)
modelFutureEnv=crop(futureEnv, model.extent)

m1 <- glm(daphnia_abundance ~ bio1, family="binomial", data=daphnia_abundance)

?glm
# map mean annual temperature
plot(modelEnv[["bio1"]]/10, main="Annual Mean Temperature")
plot(wrld_simpl,xlim=c(min(daphnia_abundance$Longitude....)-10,max(daphnia_abundance$Longitude....)+10), ylim=c(min(daphnia_abundance$Latitude)-10,max(daphnia_abundance$Latitude)+10), fill=FALSE, add=TRUE)
points(daphnia_abundance$Longitude...., daphnia_abundance$Latitude, col="red", pch=20, cex=0.75)

# map future mean annual temperature
plot(modelFutureEnv[["bio1"]]/10, main="Future Annual Mean Temperature")
plot(wrld_simpl,xlim=c(min(daphnia_abundance$Longitude....)-10,max(daphnia_abundance$Longitude....)+10), ylim=c(min(daphnia_abundance$Latitude)-10,max(daphnia_abundance$Latitude)+10), fill=FALSE, add=TRUE)
points(daphnia_abundance$Longitude...., daphnia_abundance$Latitude, col="red", pch=20, cex=0.75)

#SDM assessment
daphnia_abundance_coordinates=cbind.data.frame(daphnia_abundance$Longitude....,daphnia_abundance$Latitude) #first, just make a data frame of latitudes and longitudes for the model
fold <- kfold(daphnia_abundance_coordinates, k=5) # add an index that makes five random groups of observations
daphnia_abundance_test <- daphnia_abundance_coordinates[fold == 1, ] # hold out one fifth as test data
daphnia_abundance_train <- daphnia_abundance_coordinates[fold != 1, ] # the other four fifths are training data

#Fit SDM using the Maxent algorithm
daphnia_abundance_model <- maxent(modelEnv, daphnia_abundance_train) #note we just using the training data
plot(daphnia_abundance_model)

#How does the likelihood of species occurrence respond to variation in these climatic conditions?
response(daphnia_abundance_train)

#Predicted values for every cell in our region
daphnia_abundance_pred <- predict(daphnia_abundance_train, modelEnv)
daphnia_abundance_pred
plot(wrld_simpl, fill=FALSE, add=TRUE)
points(dmC$lon, dmC$lat, pch="+", cex=0.2)

dpC_pred <- predict(dpC_model, modelEnv)
plot(dpC_pred, main="Predicted Suitability")
plot(wrld_simpl, fill=FALSE, add=TRUE)
points(dpC$lon, dpC$lat, pch="+", cex=0.2)

#Generate background points for pseudoabsences
#To evaluate the predictive accuracy of the model, we turn back to our test data, and use cross-validation to test the model. Our evaluation will use an approach called the Area Under the Receiver Operator Curve (AUC). Basically, the process is to set thresholds on the prediction to generate different levels of the false postive rate (predicting presence for places where the species is absent), then assess the true postive rate (successfully predicting presence) as a function of the false positive rate. The area under this curve, which varies from zero to one, provides an assessment of the model. An AUC value of 0.5 is the same as random guessing of presence/absence, while values towards one mean our predictions are more reliable.
bg <- randomPoints(modelEnv, 1000)
e1 <- evaluate(dmC_model, p=dmC_test, a=bg, x=modelEnv)
plot(e1, 'ROC')

e2<-evaluate(dpC_model, p=dpC_test, a=bg, x=modelEnv)
plot(e2, 'ROC')

#Projecting responses to climate change
dmC_2070 <- predict(dmC_model, modelFutureEnv)
plot(dmC_2070, main="Predicted Future Suitability")
plot(wrld_simpl, fill=FALSE, add=TRUE)
points(dmC$lon, dmC$lat, pch="+", cex=0.2)

dpC_2070 <- predict(dpC_model, modelFutureEnv)
plot(dpC_2070, main="Predicted Future Suitability")
plot(wrld_simpl, fill=FALSE, add=TRUE)
points(dpC$lon, dpC$lat, pch="+", cex=0.2)

#difference in habitat suitability between current conditions and those in 2070
dmC_change <- dmC_2070-dmC_pred
plot(dmC_2070, main="Predicted Change in Suitability")
plot(wrld_simpl, fill=FALSE, add=TRUE)
points(dmC$lon, dmC$lat, pch="+", cex=0.2)

dpC_change <- dpC_2070-dpC_pred
plot(dpC_2070, main="Predicted Change in Suitability")
plot(wrld_simpl, fill=FALSE, add=TRUE)
points(dpC$lon, dpC$lat, pch="+", cex=0.2)

#Visualize changes in histogram
dmCChangePoints <- extract(dmC_change, dmC_coordinates)
hist(dmCChangePoints, main="")
abline(v=0, col="red")

dpCChangePoints <- extract(dpC_change, dpC_coordinates)
hist(dpCChangePoints, main="")
abline(v=0, col="red")




