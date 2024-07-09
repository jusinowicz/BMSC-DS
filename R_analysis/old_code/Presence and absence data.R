#################################Presence and absence data

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

#==============================================================================
###1.1 Data Access
#==============================================================================

#Current environment from worldclim
currentEnv=getData('worldclim', var="bio", res=2.5, download=T)
#Future environment based on the 8.5 concentration scenario for 2070 from HADGEM2-ES model
futureEnv=getData('CMIP5', var='bio', res=2.5, rcp=85, model='HE', year=70, download=T)
names(futureEnv)=names(currentEnv)

#World map:
data(wrld_simpl)

zooplankton_data<-read.csv("/Users/larac/Desktop/BMSC FP 2020/Directed Studies/zooplankton_data/Loewen et al. 2018. Ecography_Data_Dryad.csv")
colnames(zooplankton_data)


#==============================================================================
###1.2 Data Processing
#==============================================================================

#Get daphnia data
daphnia_data<-select(zooplankton_data, ID, NAME, Daphnia.ambigua:Daphnia.schodleri, LATITUDE, LONGITUDE, RECENT_DATE)
colnames(daphnia_data)

# get rid of occurences without location information or years
daphnia_data=subset(daphnia_data, !is.na(LONGITUDE) & !is.na(LATITUDE) & !is.na(RECENT_DATE)) 
daphnia_data

#daphnia_data$Daphnia.pulex<-as.numeric(daphnia_data$Daphnia.pulex)

#D.pulex records
daphnia_pulex<-select(daphnia_data, NAME, Daphnia.pulex, LONGITUDE, LATITUDE, RECENT_DATE)

#Select D. pulex, presence points
daphnia_pulex_presence<- select(daphnia_data, NAME, Daphnia.pulex, LONGITUDE, LATITUDE, RECENT_DATE) %>% 
  filter(Daphnia.pulex=="1")

#find and eliminate duplicate locations
daphnia_pulex_presence_duplicated = duplicated(daphnia_pulex_presence[, c("LONGITUDE", "LATITUDE")])
daphnia_pulex_presence = daphnia_pulex_presence[!daphnia_pulex_presence_duplicated, ] #Daphnia pulexhas 290 presence records

#Select D. pulex, absence points
daphnia_pulex_absence<- select(daphnia_data, NAME, Daphnia.pulex, LONGITUDE, LATITUDE, RECENT_DATE) %>% 
  filter(Daphnia.pulex=="0")

#find and eliminate duplicate locations
daphnia_pulex_absence_duplicated = duplicated(daphnia_pulex_absence[, c("LONGITUDE", "LATITUDE")])
daphnia_pulex_absence = daphnia_pulex_absence[!daphnia_pulex_absence_duplicated, ] #Daphnia pulex has 950 absence records

#Select D. longiremis presence points
daphnia_longiremis_presence<- select(daphnia_data, NAME, Daphnia.longiremis, LONGITUDE, LATITUDE, RECENT_DATE) %>% 
  filter(Daphnia.longiremis=="1")

#find and eliminate duplicate locations
daphnia_longiremis_presence_duplicated = duplicated(daphnia_longiremis_presence[, c("LONGITUDE", "LATITUDE")])
daphnia_longiremis_presence = daphnia_longiremis_presence[!daphnia_longiremis_presence_duplicated, ] #Daphnia longiremis has 88 presence records


#Select D.longiremis, absence points
daphnia_longiremis_absence<- select(daphnia_data, NAME, Daphnia.longiremis, LONGITUDE, LATITUDE, RECENT_DATE) %>% 
  filter(Daphnia.longiremis=="0")

#find and eliminate duplicate locations
daphnia_longiremis_absence_duplicated = duplicated(daphnia_longiremis_absence[, c("LONGITUDE", "LATITUDE")])
daphnia_longiremis_absence = daphnia_longiremis_absence[!daphnia_longiremis_absence_duplicated, ] #Daphnia pulex has 950 absence records

#==============================================================================
###1.3 Intial plots
#==============================================================================

#Initial plot, presence and absence
plot(wrld_simpl, xlim=c(min(daphnia_data$LONGITUDE)-1,max(daphnia_data$LONGITUDE)+1), ylim=c(min(daphnia_data$LATITUDE)-1,max(daphnia_data$LATITUDE)+1), axes=TRUE, col="light yellow")
points(daphnia_pulex_presence$LONGITUDE, daphnia_pulex_presence$LATITUDE, col="orange", pch=20, cex=0.75)
#points(daphnia_pulex_absence$LONGITUDE, daphnia_pulex_absence$LATITUDE, col="blue", pch=20, cex=0.75)
points(daphnia_longiremis_presence$LONGITUDE, daphnia_longiremis_presence$LATITUDE, col="green", pch=20, cex=0.75)
#points(daphnia_longiremis_absence$LONGITUDE, daphnia_longiremis_absence$LATITUDE, col="purple", pch=20, cex=0.75)

#Focus in on Canada:
#figure out extent of BC

#==============================================================================
###2.0 Create MaxEnt model
#==============================================================================

model.extent<-extent(min(daphnia_data$LONGITUDE)-10,max(daphnia_data$LONGITUDE)+10,min(daphnia_data$LATITUDE)-10,max(daphnia_data$LATITUDE)+10)
modelEnv=crop(currentEnv, model.extent)
modelFutureEnv=crop(futureEnv, model.extent)

# map mean annual temperature - example
# plot(modelEnv[["bio1"]]/10, main="Annual Mean Temperature")
# plot(wrld_simpl,xlim=c(min(daphnia_data$LONGITUDE)-10,max(daphnia_data$LONGITUDE)+10), ylim=c(min(daphnia_data$LATITUDE)-10,max(daphnia_data$LATITUDE)+10), fill=FALSE, add=TRUE)
# points(daphnia_pulex_presence$LONGITUDE, daphnia_pulex_presence$LATITUDE, col="red", pch=20, cex=0.75)
# points(daphnia_longiremis_presence$LONGITUDE, daphnia_longiremis_presence$LATITUDE, col="blue", pch=20, cex=0.75)

# map future mean annual temperature - example
# plot(modelFutureEnv[["bio1"]]/10, main="Future Annual Mean Temperature")
# plot(wrld_simpl,xlim=c(min(daphnia_data$LONGITUDE)-10,max(daphnia_data$LONGITUDE)+10), ylim=c(min(daphnia_data$LATITUDE)-10,max(daphnia_data$LATITUDE)+10), fill=FALSE, add=TRUE)
# points(daphnia_pulex_presence$LONGITUDE, daphnia_pulex_presence$LATITUDE, col="red", pch=20, cex=0.75)
# points(daphnia_longiremis_presence$LONGITUDE, daphnia_longiremis_presence$LATITUDE, col="blue", pch=20, cex=0.75)

#==============================================================================
###2.1 Model assessment
#==============================================================================
set.seed(1234)

#SDM assessment
daphnia_pulex_coordinates=cbind.data.frame(daphnia_pulex_presence$LONGITUDE,daphnia_pulex_presence$LATITUDE) #first, just make a data frame of latitudes and longitudes for the model
fold_pulex <- kfold(daphnia_pulex_coordinates, k=5) # add an index that makes five random groups of observations
daphnia_pulex_test <- daphnia_pulex_coordinates[fold_pulex == 1, ] # hold out one fifth as test data
daphnia_pulex_train <- daphnia_pulex_coordinates[fold_pulex != 1, ] # the other four fifths are training data

daphnia_longiremis_coordinates=cbind.data.frame(daphnia_longiremis_presence$LONGITUDE,daphnia_longiremis_presence$LATITUDE) #first, just make a data frame of latitudes and longitudes for the model
fold_longiremis <- kfold(daphnia_longiremis_coordinates, k=5) # add an index that makes five random groups of observations
daphnia_longiremis_test <- daphnia_longiremis_coordinates[fold_longiremis == 1, ] # hold out one fifth as test data
daphnia_longiremis_train <- daphnia_longiremis_coordinates[fold_longiremis != 1, ] # the other four fifths are training data

#==============================================================================
###2.2 Fit SDM using Maxent algorithm
#============================================================================== 

#Fit SDM using the Maxent algorithm
daphnia_pulex_model <- maxent(modelEnv, daphnia_pulex_train) #note we just using the training data
plot(daphnia_pulex_model)

daphnia_longiremis_model <- maxent(modelEnv, daphnia_longiremis_train) #note we just using the training data
plot(daphnia_longiremis_model)

#==============================================================================
###2.3 Examine sensitivity to climatic variables
#============================================================================== 

#How does the likelihood of species occurrence respond to variation in these climatic conditions?
response(daphnia_pulex_model)
response(daphnia_longiremis_model)

#==============================================================================
###2.4 Predict current suitability
#============================================================================== 

#Predicted values for every cell in our region
daphnia_pulex_pred <- predict(daphnia_pulex_model, modelEnv)
plot(daphnia_pulex_pred, main="Predicted Suitability")
plot(wrld_simpl, fill=FALSE, add=TRUE)
points(daphnia_data$LONGITUDE, daphnia_data$LATITUDE, pch="+", cex=0.2)

daphnia_longiremis_pred <- predict(daphnia_longiremis_model, modelEnv)
plot(daphnia_longiremis_pred, main="Predicted Suitability")
plot(wrld_simpl, fill=FALSE, add=TRUE)
points(daphnia_data$LONGITUDE, daphnia_data$LATITUDE, pch="+", cex=0.2)

#==============================================================================
###2.5 Absence records and model accuracy
#==============================================================================
#Generate background points for pseudoabsences
#Already have absence points so using those instead.

daphnia_pulex_absence_xy <- select(daphnia_pulex_absence, LONGITUDE:LATITUDE)
head(daphnia_pulex_absence_xy)
e1 <- evaluate(daphnia_pulex_model, p=daphnia_pulex_test, a=daphnia_pulex_absence_xy, x=modelEnv)
plot(e1, 'ROC')
#very low accuracy...

daphnia_longiremis_absence_xy <- select(daphnia_longiremis_absence, LONGITUDE:LATITUDE)
e2 <- evaluate(daphnia_longiremis_model, p=daphnia_longiremis_test, a=daphnia_longiremis_absence_xy, x=modelEnv)
plot(e2, 'ROC')

#==============================================================================
###2.6 Predict future (2070) suitability
#==============================================================================

#Projecting responses to climate change
daphnia_pulex_2070 <- predict(daphnia_pulex_model, modelFutureEnv)
plot(daphnia_pulex_2070, main="Predicted Future Suitability")
plot(wrld_simpl, fill=FALSE, add=TRUE)
points(daphnia_pulex_presence$LONGITUDE, daphnia_pulex_presence$LATITUDE, pch="+", cex=0.2)

daphnia_longiremis_2070 <- predict(daphnia_longiremis_model, modelFutureEnv)
plot(daphnia_longiremis_2070, main="Predicted Future Suitability")
plot(wrld_simpl, fill=FALSE, add=TRUE)
points(daphnia_longiremis_presence$LONGITUDE, daphnia_longiremis_presence$LATITUDE, pch="+", cex=0.2)

#==============================================================================
###3.0 Difference in habitat suitability between current and future (2070) conditions
#==============================================================================
#difference in habitat suitability between current conditions and those in 2070

daphnia_pulex_change <- daphnia_pulex_2070-daphnia_pulex_pred
plot(daphnia_pulex_2070, main="Predicted Change in Suitability")
plot(wrld_simpl, fill=FALSE, add=TRUE)
points(daphnia_pulex_presence$LONGITUDE, daphnia_pulex_presence$LATITUDE, pch="+", cex=0.2)

daphnia_longiremis_change <- daphnia_longiremis_2070-daphnia_longiremis_pred
plot(daphnia_longiremis_2070, main="Predicted Change in Suitability")
plot(wrld_simpl, fill=FALSE, add=TRUE)
points(daphnia_longiremis_presence$LONGITUDE, daphnia_longiremis_presence$LATITUDE, pch="+", cex=0.2)

#==============================================================================
###3.1 Visual change in habitat suitability between current and future (2070) conditions
#==============================================================================

#Visualize changes in histogram
daphnia_pulex_ChangePoints <- extract(daphnia_pulex_change, daphnia_pulex_coordinates)
hist(daphnia_pulex_ChangePoints, main="")
abline(v=0, col="red")

daphnia_longiremis_ChangePoints <- extract(daphnia_longiremis_change, daphnia_longiremis_coordinates)
hist(daphnia_longiremis_ChangePoints, main ="")
abline(v=0, col="red")

###################################################################################
####END OF CODE



