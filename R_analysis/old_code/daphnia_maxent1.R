#==============================================================================
#Adapt the code from kerkhoff_maxent.R example to Daphnia to take a first look
#
#
#==============================================================================

#################################Presence only data

library(raster)
library(rgdal)
library(maps)
library(mapdata)
library(dismo)  
library(rJava)  
library(maptools)
library(jsonlite)
library(rgbif) 

#system.file(package="dismo")
#system.file("java", package="dismo")

#BC boundary map
#library(ggplot2)
#install.packages("bcmaps")
#install.packages('bcmapsdata', repos='https://bcgov.github.io/drat/')
#library(bcmaps)
#available_layers()
#library(sf)
#bc <- bc_bound()
#plot(st_geometry(bc))

#Canada <- getData('GADM', country="CAN", level=1)
#plot(Canada)
#map("worldHires","Canada", xlim=c(-130,-90), ylim=c(30,60), col="gray95", fill=TRUE, add=TRUE)
#install.packages("ggExtra")
#library(ggExtra)
#install.packages("spatialEco")
#library(spatialEco)
#chelsa<- raster("Chelsa_tw.tif")

#get_chelsa <- function(type = "bioclim", layer = 1:19, period, model_string, scenario_string, future_years, output_dir)

#current <- getData('get_chelsa', type = "bioclim", layer = 1:19, period, model_string, scenario_string, future_years, output_dir)

             
#get_chelsa <- function(type = "bioclim", layer = 1:19, period, model_string, scenario_string, future_years, output_dir)

 
# Please note that you have to set download=T if you haven't downloaded the data before:
#bio_fut <- getData('CMIP5', var='bio', res=0.5, lon=5.5, lat=45.5, rcp=45, model='NO', year=50, download=T)[[c(2,5,14)]]

#==============================================================================
###1.1 Data Access
#==============================================================================
#Current environment from worldclim
currentEnv=getData('worldclim', var="bio", res=2.5, download=T)
#Future environment based on the 8.5 concentration scenario for 2070 from HADGEM2-ES model
futureEnv=getData('CMIP5', var='bio', res=2.5, rcp=85, model='HE', year=70, download=T)
names(futureEnv)=names(currentEnv)
#Species locations from gbif - for example, the timber rattlesnake
#dm = gbif('daphnia','magna*') #Daphnia magna - 1897 records
dp = gbif('daphnia', 'pulex*') #Daphnia pulex - 2022 records
#dr = gbif('daphnia','rosea*') #Daphnia rosea- 145 records
#dpu = gbif('daphnia', 'pulicaria*') #Daphnia pulicaria - 293 records
#dl = gbif('daphnia', 'lumholtzi*') #Daphnia lumholtzi - 600 records
#dmd = gbif('daphnia', 'middendorffiana*') #Daphnia middendorffiana - 140 records
dlong = gbif('daphnia', 'longiremis*') #Daphnia longiremis - 228 records 

# Please note that you have to set download=T if you haven't downloaded the data before:
#bio_curr <- getData('worldclim', var='bio', res=0.5, lon=5.5, lat=45.5, path='', download=T)[[c(2,5,14)]]

# Please note that you have to set download=T if you haven't downloaded the data before:
#bio_fut <- getData('CMIP5', var='bio', res=0.5, lon=5.5, lat=45.5, rcp=45, model='NO', year=50, download=TRUE)[[c(2,5,14)]]  
#names(bio_curr)=names(bio_fut)

#World map:
data(wrld_simpl)

#==============================================================================
###1.2 Data Processing
#==============================================================================
# limit number of predictors just a bit
#currentEnv=dropLayer(currentEnv, c("bio2", "bio3", "bio4", "bio10", "bio11", "bio13", "bio14", "bio15"))
#futureEnv=dropLayer(future_bio, c("bio2", "bio3", "bio4", "bio10", "bio11", "bio13", "bio14", "bio15"))

#bio_curr=dropLayer(bio_curr, c("bio2", "bio3", "bio4", "bio10", "bio11", "bio13", "bio14", "bio15"))
#bio_fut=dropLayer(bio_fut, c("bio2", "bio3", "bio4", "bio10", "bio11", "bio13", "bio14", "bio15"))


# get rid of occurences without location information
#dm=subset(dm, !is.na(lon) & !is.na(lat)) #now Daphnia magna has 905 records
dp=subset(dp, !is.na(lon) & !is.na(lat)) #now Daphnia pulex has 1599 records
#dr=subset(dr, !is.na(lon) & !is.na(lat)) #now Daphnia rosea has 94 records
#dpu=subset(dpu, !is.na(lon) & !is.na(lat)) #now Daphnia pulicaria has 212 records
#dl=subset(dl, !is.na(lon) & !is.na(lat)) #now Daphnia lumhotzi has 598 records
#dmd=subset(dmd, !is.na(lon) & !is.na(lat)) #now Daphnia middendorffiana has 58 records
dlong=subset(dlong, !is.na(lon)&!is.na(lat)) #204 records

# find and eliminate duplicate locations

dpdups = duplicated(dp[, c("lon", "lat")])
dp = dp[!dpdups, ] #now Daphnia pulex has 764 records

#dmdups = duplicated(dm[, c("lon", "lat")])
#dm = dm[!dmdups, ] #now Daphnia magna has 479 records

#dpudups = duplicated(dpu[, c("lon", "lat")])
#dpu = dpu[!dpudups, ] #now Daphnia pulicaria has 79 records

#drdups = duplicated(dr[, c("lon", "lat")])
#dr = dr[!drdups, ] #now Daphnia rosea has 74 records

#dldups = duplicated(dl[, c("lon", "lat")])
#dl = dl[!dldups, ] #now Daphnia lumholzi has 307 records

#dmddups = duplicated(dmd[, c("lon", "lat")])
#dmd = dmd[!dmddups, ] #now Daphnia middendorffiana has 58 records

dlongdups = duplicated(dlong[, c("lon", "lat")])
dlong = dlong[!dlongdups, ] #117 records

#==============================================================================
###1.3 Intial plots
#==============================================================================
# make initial plot (worldwide) for diagnostic purposes

plot(wrld_simpl, xlim=c(min(dp$lon)-1,max(dp$lon)+1), ylim=c(min(dp$lat)-1,max(dp$lat)+1), axes=TRUE, col="light yellow")
#points(dm$lon, dm$lat, col="orange", pch=20, cex=0.75)
points(dp$lon, dp$lat, col="blue", pch=20, cex=0.75)
#points(dr$lon, dr$lat, col="green", pch=20, cex=0.75)
#points(dpu$lon, dpu$lat, col="red", pch=20, cex=0.75)
#points(dl$lon, dl$lat, col="purple", pch=20, cex=0.75)
#points(dmd$lon, dmd$lat, col="pink", pch=20, cex=0.75)
points(dlong$lon, dlong$lat, col="red", pch=20, cex=0.75)

#eliminate questionable points likely outside the native range
# rattler = rattler[rattler$lon < -40 & rattler$lat > 25 , ]
#Focus in on Canada: 
#dmC = subset(dm, country == "Canada") #Daphnia magna - 29 records in Canada
dpC = subset(dp, country == "Canada") #Daphnia pulex - 94 observations in Canada
#drC = subset(dr, country == "Canada") #Daphnia rosea - 48 records in Canada
#dpuC = subset(dpu, country == "Canada") #Daphnia pulicaria - 38 observations in Canada
#dlC = subset(dl, country == "Canada") #Daphnia lumholzi - 4 observations in Canada
#dmdC = subset(dmd, country == "Canada") #Daphnia middendorffiana - 49 observations in Canada
dlongC = subset(dlong, country == "Canada") #22 observations in Canada

# make initial plot (Canada) for diagnostic purposes
plot(wrld_simpl, xlim=c(min(dpC$lon)-1,max(dpC$lon)+1), ylim=c(min(dpC$lat)-1,max(dpC$lat)+1), axes=TRUE, col="light yellow")

#points(dmC$lon, dmC$lat, col="orange", pch=20, cex=0.75)
points(dpC$lon, dpC$lat, col="blue", pch=20, cex=0.75)
#points(drC$lon, drC$lat, col="green", pch=20, cex=0.75)
#points(dpuC$lon, dpuC$lat, col="red", pch=20, cex=0.75)
#points(dlC$lon, dlC$lat, col="purple", pch=20, cex=0.75)
#points(dmdC$lon, dmdC$lat, col="pink", pch=20, cex=0.75)
points(dlongC$lon, dlongC$lat, col="red", pch=20, cex=0.75)


# make another plot, with higher resolution background
#map('worldHires', xlim=c(min(dm_BC$lon)-1,max(dm_BC$lon)+1), ylim=c(min(dm_BC$lat)-1,max(dm_BC$lat)+1), fill=TRUE, col="light yellow")
#points(dm_BC$lon, dm_BC$lat, col="orange", pch=20, cex=0.75)
#points(dpC$lon, dpC$lat, col="blue", pch=20, cex=0.75)
#points(drC$lon, drC$lat, col="green", pch=20, cex=0.75)
#points(dpuC$lon, dpuC$lat, col="red", pch=20, cex=0.75)
#points(dlC$lon, dlC$lat, col="red", pch=20, cex=0.75)
#points(dmC$lon, dmC$lat, col="red", pch=20, cex=0.75)

#==============================================================================
###2.0 Create MaxEnt model
#==============================================================================

model.extent<-extent(min(dpC$lon)-10,max(dpC$lon)+10,min(dpC$lat)-10,max(dpC$lat)+10)
modelEnv=crop(currentEnv, model.extent)
modelFutureEnv=crop(futureEnv, model.extent)

# map mean annual temperature - example
#plot(modelEnv[["bio1"]]/10, main="Annual Mean Temperature")
#plot(wrld_simpl,xlim=c(min(dpC$lon)-10,max(dpC$lon)+10), ylim=c(min(dpC$lat)-10,max(dpC$lat)+10), fill=FALSE, add=TRUE)
#points(dmC$lon, dmC$lat, col="orange", pch=20, cex=0.75)
#points(dpC$lon, dpC$lat, col="blue", pch=20, cex=0.75)
#points(drC$lon, drC$lat, col="green", pch=20, cex=0.75)
#points(dpuC$lon, dpuC$lat, col="red", pch=20, cex=0.75)
#points(dlC$lon, dlC$lat, col="purple", pch=20, cex=0.75)
#points(dmdC$lon, dmdC$lat, col="pink", pch=20, cex=0.75)
#points(dlongC$lon, dlongC$lat, col="red", pch=20, cex=0.75)


# map future mean annual temperature - example
#plot(modelFutureEnv[["bio1"]]/10, main="Future Annual Mean Temperature")
#plot(wrld_simpl,xlim=c(min(dpC$lon)-10,max(dpC$lon)+10), ylim=c(min(dpC$lat)-10,max(dpC$lat)+10), fill=FALSE, add=TRUE)
#points(dm$lon, dm$lat, col="orange", pch=20, cex=0.75)
#points(dpC$lon, dpC$lat, col="blue", pch=20, cex=0.75)
#points(drC$lon, drC$lat, col="green", pch=20, cex=0.75)
#points(dpuC$lon, dpuC$lat, col="red", pch=20, cex=0.75)
#points(dlC$lon, dlC$lat, col="red", pch=20, cex=0.75)
#points(dmC$lon, dmC$lat, col="red", pch=20, cex=0.75)
#points(dlongC$lon, dlongC$lat, col="red", pch=20, cex=0.75)

#==============================================================================
###2.1 Model assessment
#==============================================================================
set.seed(1234)
#SDM assessment
#dmC_coordinates=cbind.data.frame(dmC$lon,dmC$lat) #first, just make a data frame of latitudes and longitudes for the model
#fold <- kfold(dmC_coordinates, k=5) # add an index that makes five random groups of observations
#dmC_test <- dmC_coordinates[fold == 1, ] # hold out one fifth as test data
#dmC_train <- dmC_coordinates[fold != 1, ] # the other four fifths are training data


 dpC_coordinates=cbind.data.frame(dpC$lon,dpC$lat) #first, just make a data frame of latitudes and longitudes for the model
 fold <- kfold(dpC_coordinates, k=5) # add an index that makes five random groups of observations
 dpC_test <- dpC_coordinates[fold == 1, ] # hold out one fifth as test data
 dpC_train <- dpC_coordinates[fold != 1, ] # the other four fifths are training data
 
 #drC_coordinates=cbind.data.frame(drC$lon,drC$lat) #first, just make a data frame of latitudes and longitudes for the model
 #fold <- kfold(drC_coordinates, k=5) # add an index that makes five random groups of observations
 #drC_test <- drC_coordinates[fold == 1, ] # hold out one fifth as test data
 #drC_train <- drC_coordinates[fold != 1, ] # the other four fifths are training data

 dlongC_coordinates=cbind.data.frame(dlongC$lon,dlongC$lat) #first, just make a data frame of latitudes and longitudes for the model
 fold <- kfold(dlongC_coordinates, k=5) # add an index that makes five random groups of observations
 dlongC_test <- dlongC_coordinates[fold == 1, ] # hold out one fifth as test data
 dlongC_train <- dlongC_coordinates[fold != 1, ] # the other four fifths are training data
 
 #==============================================================================
 ###2.2 Fit SDM using Maxent algorithm
 #============================================================================== 

 #dmC_model <- maxent(modelEnv, dmC_train) #note we just using the training data
 #plot(dmC_model)

dpC_model <- maxent(modelEnv, dpC_train) #note we just using the training data
plot(dpC_model)

#drC_model <- maxent(modelEnv, drC_train) #note we just using the training data
#plot(drC_model)

dlongC_model <- maxent(modelEnv, dlongC_train) #note we just using the training data
plot(dlongC_model)

#system.file(package="dismo")

#==============================================================================
###2.3 Examine sensitivity to climatic variables
#============================================================================== 

#How does the likelihood of species occurrence respond to variation in these climatic conditions?
#response(dmC_model)
response(dpC_model)
#response(drC_model)
response(dlongC_model)

#==============================================================================
###2.4 Predict current suitability
#============================================================================== 

#Predicted values for every cell in our region
#dmC_pred <- predict(dmC_model, modelEnv)
#plot(dmC_pred, main="Predicted Suitability")
#plot(wrld_simpl, fill=FALSE, add=TRUE)
#points(dmC$lon, dmC$lat, pch="+", cex=0.2)


dpC_pred <- predict(dpC_model, modelEnv)
plot(dpC_pred, main="Predicted Suitability")
plot(wrld_simpl, fill=FALSE, add=TRUE)
points(dpC$lon, dpC$lat, pch="+", cex=0.2)

#drC_pred <- predict(drC_model, modelEnv)
#plot(drC_pred, main="Predicted Suitability")
#plot(wrld_simpl, fill=FALSE, add=TRUE)
#points(drC$lon, drC$lat, pch="+", cex=0.2)

dlongC_pred <- predict(dlongC_model, modelEnv)
plot(dlongC_pred, main="Predicted Suitability")
plot(wrld_simpl, fill=FALSE, add=TRUE)
points(dlongC$lon, dlongC$lat, pch="+", cex=0.2)

#==============================================================================
###2.5 Pseudoabsences and model accuracy
#==============================================================================

#Generate background points for pseudoabsences
#To evaluate the predictive accuracy of the model, we turn back to our test data, and use cross-validation to test the model. Our evaluation will use an approach called the Area Under the Receiver Operator Curve (AUC). Basically, the process is to set thresholds on the prediction to generate different levels of the false postive rate (predicting presence for places where the species is absent), then assess the true postive rate (successfully predicting presence) as a function of the false positive rate. The area under this curve, which varies from zero to one, provides an assessment of the model. An AUC value of 0.5 is the same as random guessing of presence/absence, while values towards one mean our predictions are more reliable.
bg <- randomPoints(modelEnv, 1000)
#e1 <- evaluate(dmC_model, p=dmC_test, a=bg, x=modelEnv)
#plot(e1, 'ROC')

e2<-evaluate(dpC_model, p=dpC_test, a=bg, x=modelEnv)
plot(e2, 'ROC')

#e3<-evaluate(drC_model, p=drC_test, a=bg, x=modelEnv)
#plot(e3, 'ROC')

e4<-evaluate(dlongC_model, p=dlongC_test, a=bg, x=modelEnv)
plot(e4, 'ROC')

#==============================================================================
###2.6 Predict future (2070) suitability
#==============================================================================

#Projecting responses to climate change
#dmC_2070 <- predict(dmC_model, modelFutureEnv)
#plot(dmC_2070, main="Predicted Future Suitability")
#plot(wrld_simpl, fill=FALSE, add=TRUE)
#points(dmC$lon, dmC$lat, pch="+", cex=0.2)

dpC_2070 <- predict(dpC_model, modelFutureEnv)
plot(dpC_2070, main="Predicted Future Suitability")
plot(wrld_simpl, fill=FALSE, add=TRUE)
points(dpC$lon, dpC$lat, pch="+", cex=0.2)

#drC_2070 <- predict(drC_model, modelFutureEnv)
#plot(drC_2070, main="Predicted Future Suitability")
#plot(wrld_simpl, fill=FALSE, add=TRUE)
#points(drC$lon, drC$lat, pch="+", cex=0.2)

dlongC_2070 <- predict(dlongC_model, modelFutureEnv)
plot(dlongC_2070, main="Predicted Future Suitability")
plot(wrld_simpl, fill=FALSE, add=TRUE)
points(dlongC$lon, dlongC$lat, pch="+", cex=0.2)

#==============================================================================
###3.0 Difference in habitat suitability between current and future (2070) conditions
#==============================================================================

#dmC_change <- dmC_2070-dmC_pred
#plot(dmC_2070, main="Predicted Change in Suitability")
#plot(wrld_simpl, fill=FALSE, add=TRUE)
#points(dmC$lon, dmC$lat, pch="+", cex=0.2)

dpC_change <- dpC_2070-dpC_pred
plot(dpC_2070, main="Predicted Change in Suitability")
plot(wrld_simpl, fill=FALSE, add=TRUE)
points(dpC$lon, dpC$lat, pch="+", cex=0.2)

#drC_change <- drC_2070-dpC_pred
#plot(drC_2070, main="Predicted Change in Suitability")
#plot(wrld_simpl, fill=FALSE, add=TRUE)
#points(drC$lon, drC$lat, pch="+", cex=0.2)

dlongC_change <- dlongC_2070-dlongC_pred
plot(dlongC_2070, main="Predicted Change in Suitability")
plot(wrld_simpl, fill=FALSE, add=TRUE)
points(dlongC$lon, dlongC$lat, pch="+", cex=0.2)

#==============================================================================
###3.1 Visual change in habitat suitability between current and future (2070) conditions
#==============================================================================
#Visualize changes in histogram
#dmCChangePoints <- extract(dmC_change, dmC_coordinates)
#hist(dmCChangePoints, main="")
#abline(v=0, col="red")

dpCChangePoints <- extract(dpC_change, dpC_coordinates)
hist(dpCChangePoints, main ="")
abline(v=0, col="red")

#drCChangePoints <- extract(drC_change, drC_coordinates)
#hist(drCChangePoints, main ="")
#abline(v=0, col="red")

dlongCChangePoints <- extract(dlongC_change, dlongC_coordinates)
hist(dlongCChangePoints, main ="")
abline(v=0, col="red")

###################################################################################
####END OF CODE