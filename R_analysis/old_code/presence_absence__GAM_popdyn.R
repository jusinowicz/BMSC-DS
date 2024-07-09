#==============================================================================
# Model species distributions based on presence/absence, using a GAM framework and competition
#==============================================================================
library(raster)
library(mgcv) # gam
library(dismo)
library(rgdal)
library(rJava)
library(XML)
library(rgdal)
library(maps)
library(mapdata)
library(dismo)  
library(maptools)
library(jsonlite)
library(rgbif) 
library(rlist)


library(rms) 
library(randomForest)
library(ellipse)
library(XML)

#==============================================================================
# Choose which species to run: 
#==============================================================================
d_species = c('Daphnia.pulex', 'Daphnia.longiremis')
#'Daphnia.middendorffiana','Daphnia.ambigua','Daphnia.dentifera', 'Daphnia.schodleri'
nspp = length(d_species)
daph_loewen = vector("list", nspp) #Vector of lists to store Daphnia data

#==============================================================================
#Data Access:
#==============================================================================
#Current environment from worldclim
currentEnv=raster::getData('worldclim', var="bio", res=2.5, download=T)
#Future environment based on the 8.5 concentration scenario for 2070 from HADGEM2-ES model
futureEnv=raster::getData('CMIP5', var='bio', res=2.5, rcp=85, model='HE', year=70, download=T)
names(futureEnv)=names(currentEnv)

#World map:
data(wrld_simpl)

#Daphnia data from Loewen et al. (2018)
for (i in 1:nspp ) {
  daph_loewen[[i]]<-read.csv("./Loewen et al. 2018. Ecography_Data_Dryad.csv") %>% 
    select(d_species[i], 'LATITUDE', 'LONGITUDE', 'RECENT_DATE', 'NAME')
}

# plot(raster(currentEnv, 1))
# plot(wrld_simpl, add=TRUE)

#===============================================================
### Drop variables and visual inspection of collinearity ###
#===========================================================================

#drop precipipation variables
currentEnv=dropLayer(currentEnv, c("bio16", "bio17", "bio18", "bio19", "bio12", "bio13", "bio14", "bio15"))
futureEnv=dropLayer(futureEnv, c("bio16", "bio17", "bio18", "bio19", "bio12", "bio13", "bio14", "bio15"))

current_col <- cor(getValues(currentEnv), use = "complete.obs")
#corrplot(current_col, method="ellipse")
#ellipse going from lower left to upper right is positive correlation
#ellipse going from bottom right to upper left is negative correlation
#width indicates strength of correlation
#straight line = perfect correlation
corrplot(current_col, order="AOE", method="color", addCoef.col="gray")
#dark red is strong negative correlation, dark blue is strong positive correlation

### Select an uncorrelated subset of environmental variables ###
#currentEnv<-subset(currentEnv, c("bio3", "bio4", "bio6", "bio11", "bio9", "bio13", "bio14", "bio15"))

future_col <- cor(getValues(futureEnv), use = "complete.obs")
corrplot(future_col, order="AOE", method="color", addCoef.col="gray")
### Select an uncorrelated subset of environmental variables ###
#futureEnv<-subset(futureEnv, c("bio2", "bio3", "bio4", "bio10", "bio11", "bio13", "bio14", "bio15"))


#==============================================================================
# Processing points:
#==============================================================================
daph_spp = vector("list",nspp) #Vector of lists to store tidy daphnia data

# get rid of occurences without location information
for (i in 1:nspp ) {
  
  daph_spp[[i]] = subset(daph_loewen[[i]], !is.na(LONGITUDE) & !is.na(LATITUDE) & !is.na(RECENT_DATE))
  
  # find and eliminate duplicate locations
  dpdups = duplicated(daph_spp[[i]][, c("LONGITUDE", "LATITUDE")])
  daph_spp[[i]]  = daph_spp[[i]] [!dpdups, ]
  
}

###Select species records for only where environmental info is available?
# for (i in 1:nspp ) {
#   
#   daph_spp[[i]]<- daph_spp[[i]][complete.cases(raster::extract(currentEnv, daph_spp)),]
#   
# }

###Add projection info?
#proj4string(daph_spp[i]) <- CRS("+proj=longlat +datum=WGS84")

###Separate present and absent records
#Present:
daph_present = vector("list",nspp) #Vector of lists to store daphnia presence records

for (i in 1:nspp ) {
  
  daph_present[[i]] = subset(daph_spp[[i]], daph_spp[[i]][1] == 1)
}

#Absent:
daph_absent = vector("list",nspp) #Vector of lists to store daphnia absence records

for (i in 1:nspp ) {
  
  daph_absent[[i]] = subset(daph_spp[[i]], daph_spp[[i]][1] == 0)
  
}

#==============================================================================
###1.3 Initial plots
#==============================================================================
# make initial plots (worldwide) for diagnostic purposes
#use presence records
col_use = c("red","blue")
legend_use = c("D.pulex", "D.longiremis")

plot(wrld_simpl, xlim=c(min(daph_spp[[1]]$LONGITUDE)-1,max(daph_spp[[1]]$LONGITUDE)+1), 
     ylim=c(min(daph_spp[[1]]$LATITUDE)-1,max(daph_spp[[1]]$LATITUDE)+1), axes=TRUE, col="light yellow")
for (i in 1:nspp ) {
  
  points(daph_present[[i]]$LONGITUDE, daph_present[[i]]$LATITUDE, col=col_use[i], pch=20, cex=0.75)
  
}

legend("topright", 
       legend = legend_use, text.font=c(3),
       col = col_use, 
       pch = c(19,19), 
       bty = "n", 
       cex = 1.2, 
       text.col = "black", 
       inset = c(0.01, 0.01))

#==============================================================================
### Training and testing sets
#==============================================================================
set.seed(2)

###Select only one presence record in each cell of the environmental layer?

# for (i in 1:nspp ) {
# daph_present[[i]]= gridSample(daph_present[[i]], currentEnv, n =1)
# }
#error

#Combine  presence and background points
for (i in 1:nspp ) {
  fulldaph[[i]]  = SpatialPointsDataFrame(rbind(daph_present[[i]], daph_absent[[i]]),
                                  data=data.frame("species"=rep(c(1,0),
                                                               c(nrow(presence), nrow(background)))),
                                  match.ID=FALSE, 
                                  proj4string = CRS(projection(env)))
}
#Error

###Add environmental data at point locations

for (i in 1:nspp ) {
  
  fulldaph[[i]]@data=cbind(fulldaph[[i]]@data, extract(currentEnv, fulldaph[[i]]))
}
#Error

### Split data set into training and test data 
set.seed(2)
daph_coordinates = vector("list",nspp) #Vector of lists to store data
fold = vector("list",nspp) #Vector of lists to store data
daph_test = vector("list",nspp) #Vector of lists to store data
daph_train = vector("list",nspp) #Vector of lists to store data

for (i in 1:nspp ) {
  
  #first, just make a data frame of latitudes and longitudes for the model
  daph_coordinates[[i]] = cbind.data.frame(fulldaph[[i]]$LONGITUDE,fulldaph[[i]]$LATITUDE)
  # add an index that makes five random groups of observations
  fold[[i]] = kfold(daph_coordinates[[i]], k=5) 
  # hold out one fifth as test data
  daph_test[[i]] = daph_coordinates[[i]][fold[[i]] == 1, ] 
  # the other four fifths are training data
  daph_train[[i]] = daph_coordinates[[i]][fold[[i]] != 1, ] 
}

variable_names = c("bio2", "bio3", "bio4", "bio10", "bio11", "bio13", "bio14", "bio15")

#==============================================================================
###Make GAM 
#=======================================================================
##Generalized additive models
daph_model = gam (daph_present ~ + s(bio2)+s (bio14) +s(other predictors i'm using'), data =daph_train, family="binomial")
summary(daph_model)
plot(daph_model)

#Evaluate model on test data
#Predict to test data
daph_gam_test = predict(daph_model, newdata=daph_test, type ="response")
#calculate performace
val.prob(daph_gam_test, daph_test[["daph_present"]])

#variable importance
var_imp= varImpBiomod(daph_model, variable_names,
                     daph_test)
barplot()

#======================================================
##Map predictive map
#======================================================
gam_map<- predict(currentEnv, daph_model, type="response")
plot(gam_map)


#==============================================================================
###3 Apply population dynamics (including competition), using the 
#	 maxent suitability maps to describe species intrinsic reproductive fitness
#==============================================================================
#=============================================================================
# 3.1 Discrete-time Leslie-Gower model
#=============================================================================
#Time 
tend = 100
delta1 = 1
times  = seq(from = 0, to = tend, by = delta1)
tl = length(times)

#Use the version of as.matrix from "raster" to convert raster to a matrix. 
#Do this here initially to make it easier to get spatial extents: 
l1_tmp = raster::as.matrix(daphC_pred[[1]])
nx = ncol(l1_tmp)
ny = nrow(l1_tmp)

#Get each species' lambda: 
l1_tmp = array( matrix(0,ny,nx), dim = c(ny,nx,nspp )  )
co_suit = matrix(1,ny,nx)
for ( i in 1:nspp) {
  
  l1_tmp[,,i] = raster::as.matrix(daphC_pred[[i]])
  co_suit = co_suit*l1_tmp[,,i] 
}
co_suit[is.na(co_suit)] = 0

#Set the reproduction rate: 
l1_a = 5

#Constant lambda across time
lambda_i = l1_tmp *l1_a

#Turn this into a projection across time
# lambda_i = array( l1_tmp, dim = c(ny,nx,nspp,tend) )
# rm(l1_tmp) #Just get rid of this since it is a large object. 

#Create the population matrix with the same times as continuous-time model: 
Nt = array( 0, dim = c(ny,nx,nspp,tend ))

#Initial conditions (applied to every spatial point):
Nt[,,,1] = c(0.1,0.1)

### For now, assume that competition and survival are independent of space
#Competition: 
#This creates a square matrix for aij with intraspecific competition = 1 on 
#the diagonal. 
aii1 = 1 
aij1 = 0.8
aij = matrix ( aij1, nspp, nspp )
diag(aij) = matrix(aii1,nspp,1)

#Survival (should be less than 1)
s_i = c(matrix(0.9,nspp,1))

#Loop over all of time: 
for (t in 1:(tend-1) ) {
  
  #Loop over all species: 
  for (s in 1:nspp) {
    
    #Multiply each competitor by aij using sweep( )
    intra_comp = sweep( array(Nt[, , -s, t], dim =c( dim( Nt[, , -s, t])[1], dim( Nt[, , -s, t] )[2], 
      (nspp-1)) ), 3, t(aij[s, (-s)]), FUN = "*")
    #A handy trick with rowSums to sum over 3rd dimension
    intra_comp = rowSums(intra_comp, dim=2) 
    
    Nt[,,s,t+1] = ( lambda_i[,,s] / (1+ aij[s,s] *Nt[,,s,t]+ intra_comp) + 
                      s_i[s] ) * Nt[,,s,t]
    Nt[,,s,t+1][is.na(Nt[,,s,t+1])] = 0
    
  }
  
}


###Total population 
pop_tot = matrix(0,tend,nspp)
for (i in 1:nspp) {
  for (t in 1:tend) {	
    pop_tot[t,i] = sum( colSums (Nt[,,i,t]) )
  }
}

###Make a suitability metric out of the population densities. Use the max to scale
###Max population 
Nt_suit = Nt 
pop_max = matrix(0,tend,nspp)
for (i in 1:nspp) {
  for (t in 1:tend) {	
    pop_max[t,i] = max( max (Nt[,,i,t]) )
    Nt_suit[,,i,t] =Nt[,,i,t] / pop_max[t,i] 
  }
}

suit_final = Nt_suit[,,,tend]
suit_change = suit_final - l1_tmp 

###Convert these to rasters with same qualities as original rasters:
daphC_suit_final = vector("list",nspp)
daphC_suit_change = vector("list",nspp)
daphC_cosuit = vector("list",nspp)
daphC_comp = vector("list",nspp)

for ( i in 1:nspp) {
  
  daphC_suit_final[[i]] =daphC_pred[[i]]
  daphC_suit_change[[i]] = daphC_pred[[i]]
  daphC_cosuit[[i]] =daphC_pred[[i]]
  #daphC_comp[[i]] = daphC_pred[[i]]
  
  daphC_suit_final[[i]]@data = raster(suit_final[,,i])@data
  daphC_suit_change[[i]]@data = raster(suit_change[,,i])@data
  daphC_cosuit[[i]]@data = raster(log(10^3*co_suit+1))@data
  #daphC_comp[[i]]@data = raster(suit_change[,,i])@data
  
  
}


#==============================================================================
# Compare and Plot population density-suitability .
#==============================================================================
#MaxEnt suitability: 
fig.name = paste("daphnia_maxent_suit.pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  
  plot(daphC_pred[[i]], main="Predicted Suitability")
  plot(wrld_simpl, fill=FALSE, add=TRUE)
  points(daph_sppC[[i]]$lon, daph_sppC[[i]]$lat, pch="+", cex=0.2)
  
}

dev.off()

#MaxEnt co-suitability: 
fig.name = paste("daphnia_maxent_cosuit.pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  
  plot(daphC_cosuit[[i]], main="Predicted Suitability")
  plot(wrld_simpl, fill=FALSE, add=TRUE)
  points(daph_sppC[[i]]$lon, daph_sppC[[i]]$lat, pch="+", cex=0.2)
  
}

dev.off()

#Final population-driven suitability
fig.name = paste("daphnia_pop_driven_suit.pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  
  plot(daphC_suit_final[[i]], main="Predicted Suitability")
  plot(wrld_simpl, fill=FALSE, add=TRUE)
  points(daph_sppC[[i]]$lon, daph_sppC[[i]]$lat, pch="+", cex=0.2)
  
}

dev.off()


#Final change in suitability
fig.name = paste("daphnia_pop_driven_change.pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  
  plot(daphC_suit_change[[i]], main="Predicted Suitability")
  plot(wrld_simpl, fill=FALSE, add=TRUE)
  points(daph_sppC[[i]]$lon, daph_sppC[[i]]$lat, pch="+", cex=0.2)
  
}

dev.off()
