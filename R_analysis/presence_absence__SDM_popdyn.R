#==============================================================================
# Fit an SDM with maxent and then project population dynamics (with competition) 
# onto the range projections. 
# 
# 1. Fit an SDM with maxent in the usual way. Get the suitability map. Code
#	 adapted from kerkhoff_maxent.R example
# 2. Transform the maxent suitability map into a fitness metric by scaling 
#	 relative to a known fitness measurement for these species. 
# 3. Run competitive population dynamics in each grid cell. 
#
#
#==============================================================================

#################################Presence and absence data
options(java.parameters = "-Xmx4g")

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
library(rlist)

#For database management
library(sf)
library(RSQLite)

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


daph_present = vector("list",nspp) #Vector of lists to store daphnia presence records

for (i in 1:nspp ) {
  
  daph_present[[i]] = subset(daph_spp[[i]], daph_spp[[i]][1] == 1)
}
head(daph_present)

daph_absent = vector("list",nspp) #Vector of lists to store daphnia absence records

for (i in 1:nspp ) {
  
  daph_absent[[i]] = subset(daph_spp[[i]], daph_spp[[i]][1] == 0)

}



# limit number of predictors just a bit
#currentEnv=dropLayer(currentEnv, c("bio2", "bio3", "bio4", "bio10", "bio11", "bio13", "bio14", "bio15"))
#futureEnv=dropLayer(future_bio, c("bio2", "bio3", "bio4", "bio10", "bio11", "bio13", "bio14", "bio15"))

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
#Focus in on Canada: 
#==============================================================================
# daph_sppC = vector("list",nspp) #Vector of lists to store data
# for (i in 1:nspp ) {
#   daph_sppC[[i]] = subset(daph_spp[[i]], country == "Canada")
# }
# 
# # make initial plot (Canada) for diagnostic purposes
# plot(wrld_simpl, xlim=c(min(daph_present[[i]]$lon)-1,max(daph_present[[i]]$lon)+1), 
#      ylim=c(min(daph_present[[i]]$lat)-1,max(daph_present[[i]]$lat)+1), axes=TRUE, col="light yellow")
# 
# for (i in 1:nspp ) {
#   
#   points(daph_sppC[[i]]$lon, daph_sppC[[i]]$lat, col=col_use[i], pch=20, cex=0.75)
#   
# }
# 
# legend("topright", 
#        legend = legend_use, text.font=c(3),
#        col = col_use, 
#        pch = c(19,19), 
#        bty = "n", 
#        cex = 1.2, 
#        text.col = "black", 
#        inset = c(0.03, 0.01))

#==============================================================================
###2.0 Create MaxEnt model
#==============================================================================
model.extent = vector("list",nspp) #Vector of lists to store data
modelEnv = vector("list",nspp) #Vector of lists to store data
modelFutureEnv = vector("list",nspp) #Vector of lists to store data

#max and mins for d.pulex in model.extent list
 x_lims = c(	-147.7167,	 -103.746)
 y_lims = c(	30.43043,			73.75)
 
 #max and mins for lat/long of d.pulux
# x_lims = c(	-138.7512, -113.6869)
# y_lims = c(	40.43043, 66.17637)

for (i in 1:nspp ) {
  
  model.extent[[i]] = extent(min(daph_present[[i]]$LONGITUDE)-10,max(daph_present[[i]]$LONGITUDE)+10,
                             min(daph_present[[i]]$LATITUDE)-10,max(daph_present[[i]]$LATITUDE)+10)
  modelEnv[[i]] = crop(currentEnv, model.extent[[i]])
  modelFutureEnv[[i]] =  crop(futureEnv, model.extent[[i]])
  
  #Get the overall min/max extents. This is to create one common extent 
   x_lims = c( min( c( x_lims[1], xmin(model.extent[[i]] ) ) ) , 
               max(c( x_lims[2], xmax(model.extent[[i]] ) ) )
   )
   y_lims = c( min( c( y_lims[1], ymin(model.extent[[i]] ) ) ) , 
               max(c( y_lims[2], ymax(model.extent[[i]] ) ) )
   )
  
}


#==============================================================================
###2.1 Model assessment
#==============================================================================
set.seed(1234)

daph_coordinates = vector("list",nspp) #Vector of lists to store data
fold = vector("list",nspp) #Vector of lists to store data
daph_test = vector("list",nspp) #Vector of lists to store data
daph_train = vector("list",nspp) #Vector of lists to store data

for (i in 1:nspp ) {
  
  #first, just make a data frame of latitudes and longitudes for the model
  daph_coordinates[[i]] = cbind.data.frame(daph_present[[i]]$LONGITUDE,daph_present[[i]]$LATITUDE)
  # add an index that makes five random groups of observations
  fold[[i]] = kfold(daph_coordinates[[i]], k=5) 
  # hold out one fifth as test data
  daph_test[[i]] = daph_coordinates[[i]][fold[[i]] == 1, ] 
  # the other four fifths are training data
  daph_train[[i]] = daph_coordinates[[i]][fold[[i]] != 1, ] 
}

#==============================================================================
###2.2 Fit SDM using Maxent algorithm
#============================================================================== 

daph_model = vector("list",nspp)

par (mfrow = c(ceiling(nspp/2),2))
for (i in 1:nspp ) {
  
  daph_model[[i]] = maxent(modelEnv[[i]], daph_train[[i]])
  plot(daph_model[[i]],
       xlab="Variable contribution (%)", 
       ylab="Bioclimatic variables", 
       main=" ")
  
}

#==============================================================================
###2.3 Examine sensitivity to climatic variables
#============================================================================== 

#How does the likelihood of species occurrence respond to variation 
#in these climatic conditions?
response (daph_model[[i]])

#==============================================================================
###2.4 Predict current suitability
#============================================================================== 
daph_pred = vector("list",nspp)

par (mfrow = c(ceiling(nspp/2),2))
for (i in 1:nspp){
  
  #Note: I'm projecting these to the same modelEnv to get a common
  #resolution and extent. 
  daph_pred[[i]]  = predict(daph_model[[i]], modelEnv[[1]])
  plot(daph_pred[[i]], main="Predicted Suitability")
  plot(wrld_simpl, fill=FALSE, add=TRUE)
  points(daph_present[[i]]$LONGITUDE, daph_present[[i]]$LATITUDE, pch="+", cex=0.2)
  
  #Project current suitability to the same extent: 
  # extent( daphC_pred[[i]]) = c(x_lims,y_lims)
  # res( daphC_pred[[i]] ) = res (daphC_pred[[1]])
  
}

#==============================================================================
###2.5 Pseudoabsences and model accuracy
#==============================================================================
#Generate background points for pseudoabsences
#Already have absence points so using those instead.

#bg = vector("list",nspp)
ev_bg = vector("list",nspp)

#daph_absent can only have 2 columns (lat/long) to be used in evaluate function
for (i in 1:nspp ) {
  
  daph_absent[[1]]$Daphnia.pulex=NULL
  daph_absent[[2]]$Daphnia.longiremis=NULL
  daph_absent[[i]]$RECENT_DATE=NULL
  daph_absent[[i]]$NAME=NULL
  
}


par (mfrow = c(ceiling(nspp),2))
for (i in 1:nspp) {
  
  #bg[[i]] = randomPoints(modelEnv[[i]], 1000)
  ev_bg[[i]] = evaluate (daph_model[[i]], p = daph_test[[i]], 
                         a = daph_absent[[i]], x = modelEnv[[i]] )
  plot (ev_bg[[i]], 'ROC' )
  
}

######Stuck

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
l1_tmp = raster::as.matrix(daph_pred[[1]])
nx = ncol(l1_tmp)
ny = nrow(l1_tmp)

#Get each species' lambda: 
l1_tmp = array( matrix(0,ny,nx), dim = c(ny,nx,nspp )  )
co_suit = matrix(1,ny,nx)
for ( i in 1:nspp) {
  
  l1_tmp[,,i] = raster::as.matrix(daph_pred[[i]])
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
aij1 = 1
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
    #Checks on NAs and negatives: 
    Nt[,,s,t+1][is.na(Nt[,,s,t+1])] = 0
    
    #The model with no competition
    Nt_base [,,s,t+1] =  ( lambda_i[,,s]/(1+ aij[s,s] *Nt_base[,,s,t]) + s_i[s] ) * Nt_base[,,s,t]
    Nt_base[,,s,t+1][is.na(Nt_base[,,s,t+1])] = 0
    
    # Nt[,,s,t+1] = ( lambda_i[,,s] / tot_comp + s_i[s] ) * Nt[,,s,t]
    # Nt[,,s,t+1][is.na(Nt[,,s,t+1])] = 0
    
  }
  
}

###Total population 
pop_tot = matrix(0,tend,nspp)
for (i in 1:nspp) {
  for (t in 1:tend) {	
    pop_tot[t,i] = sum( colSums (Nt[,,i,t]) )
  }
}

Nt_base - no competition
#Nt - competition
max(Nt_base[,,,tend],na.rm=T)
max(Nt[,,,tend],na.rm=T)

#No competition - equivalent to maxent suitability
Nt_base_pop = Nt_base[,,,tend]
dim(Nt_base_pop)

#Nt_base_pop = Nt_base

#Interspecific competition
Nt_pop = Nt[,,,tend]
dim(Nt_pop)

#change in population between no competition and competition
#delta_pop = (daphC_Nt_base - Nt_base_pop)/daphC_Nt_base
delta_pop = (Nt_base[,,,tend] - Nt[,,,tend])/Nt_base[,,,tend]
dim(delta_pop)


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
suit_change = suit_final -l1_tmp 

###Convert these to rasters with same qualities as original rasters:
daphC_suit_final = vector("list",nspp)
daphC_suit_change = vector("list",nspp)
daphC_delta_pop = vector("list", nspp)
daphC_Nt_base = vector("list", nspp)
daphC_Nt = vector("list", nspp)

for ( i in 1:nspp) {
  
  daphC_suit_final[[i]] =daph_pred[[i]]
  daphC_suit_change[[i]] = daph_pred[[i]]
  daphC_delta_pop[[i]]=daph_pred[[i]]
  daphC_Nt_base[[i]]=daph_pred[[i]]
  daphC_Nt[[i]]=daph_pred[[i]]

  daphC_suit_final[[i]]@data = raster(suit_final[,,i])@data
  daphC_suit_change[[i]]@data = raster(suit_change[,,i])@data
  daphC_delta_pop[[i]]@data=raster(delta_pop[,,i])@data
  daphC_Nt_base[[i]]@data=raster(Nt_base_pop[,,i])@data
  daphC_Nt[[i]]@data=raster(Nt_pop[,,i])@data
  
  
}

#==============================================================================
# Save these as raster files for GIS.
#==============================================================================

for (n in 1:nspp){ 
  output_suit_final = paste("./../QGIS/",d_species[n],"_suit_final.tif",sep="")
  writeRaster(daphC_suit_final[[i]], filename = output_suit_final, format = "GTiff", overwrite = TRUE)

  output_suit_change = paste("./../QGIS/",d_species[n],"_suit_change.tif",sep="")
  writeRaster(daphC_suit_change[[i]] , filename = output_suit_change, format = "GTiff", overwrite = TRUE)

  output_delta_pop = paste("./../QGIS/",d_species[n],"_delta_pop.tif",sep="")
  writeRaster(daphC_delta_pop[[i]], filename =  output_delta_pop, format = "GTiff", overwrite = TRUE)

  output_Nt_base = paste("./../QGIS/",d_species[n],"_suit_final.tif",sep="")
  writeRaster(daphC_Nt_base[[i]], filename = output_Nt_base, format = "GTiff", overwrite = TRUE)

  output_NT = paste("./../QGIS/",d_species[n],"_NT.tif",sep="")
  writeRaster(daphC_Nt[[i]], filename =  output_NT , format = "GTiff", overwrite = TRUE)
}

#==============================================================================
# Save in SQLite database 
#==============================================================================


#==============================================================================
# Compare and Plot population density-suitability .
#==============================================================================
#MaxEnt suitability: 
fig.name = paste("daphnia_maxent_suit.pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  
  plot(daph_pred[[i]], main="Predicted Suitability")
  plot(wrld_simpl, fill=FALSE, add=TRUE)
  points(daph_present[[i]]$LONGITUDE, daph_present[[i]]$LATITUDE, pch="+", cex=0.2)
  
}

dev.off()

#MaxEnt co-suitability: 
fig.name = paste("daphnia_maxent_cosuit.pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  
  plot(daph_cosuit[[i]], main="Predicted Suitability")
  plot(wrld_simpl, fill=FALSE, add=TRUE)
  points(daph_present[[i]]$LONGITUDE, daph_present[[i]]$LATITUDE, pch="+", cex=0.2)
  
}

dev.off()

#Final population-driven suitability
fig.name = paste("daphnia_pop_driven_suit.pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  
  plot(daphC_suit_final[[i]], main="Predicted Suitability")
  plot(wrld_simpl, fill=FALSE, add=TRUE)
  points(daph_present[[i]]$LONGITUDE, daph_present[[i]]$LATITUDE, pch="+", cex=0.2)
  
}

dev.off()


#Final change in suitability
fig.name = paste("daphnia_pop_driven_change.pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  
  plot(daph_suit_change[[i]], main="Predicted Suitability")
  plot(wrld_simpl, fill=FALSE, add=TRUE)
  points(daph_present[[i]]$LONGITUDE, daph_present[[i]]$LATITUDE, pch="+", cex=0.2)
  
}

dev.off()


#==============================================================================
### Predict future (2070) suitability
#==============================================================================

daph_pred_future = vector("list",nspp)

par (mfrow = c(ceiling(nspp/2),2))
for (i in 1:nspp){
  
  #Note: I'm projecting these to the same modelFutureEnv to get a common
  #resolution and extent. 
  daph_pred_future[[i]]  = predict(daph_model[[i]], modelFutureEnv[[1]])
  plot(daph_pred_future[[i]], main="Predicted Future Suitability")
  plot(wrld_simpl, fill=FALSE, add=TRUE)
  points(daph_present[[i]]$LONGITUDE, daph_present[[i]]$LATITUDE, pch="+", cex=0.2)
  
  #Project current suitability to the same extent: 
  # extent( daphC_pred_future[[i]]) = c(x_lims,y_lims)
  # res( daphC_pred_future[[i]] ) = res (daphC_pred_future[[1]])
  
}

#==============================================================================
### Apply population dynamics (including competition), using the 
#	 maxent suitability maps to describe species intrinsic reproductive fitness
#==============================================================================
#=============================================================================
# Discrete-time Leslie-Gower model
#=============================================================================
#Time 
tend = 100
delta1 = 1
times  = seq(from = 0, to = tend, by = delta1)
tl = length(times)

#Use the version of as.matrix from "raster" to convert raster to a matrix. 
#Do this here initially to make it easier to get spatial extents: 
l1_tmp = raster::as.matrix(daph_pred_future[[1]])
nx = ncol(l1_tmp)
ny = nrow(l1_tmp)

#Get each species' lambda: 
l1_tmp = array( matrix(0,ny,nx), dim = c(ny,nx,nspp )  )
co_suit_fut = matrix(1,ny,nx)
for ( i in 1:nspp) {
  
  l1_tmp[,,i] = raster::as.matrix(daph_pred_future[[i]])
  co_suit_fut = co_suit_fut*l1_tmp[,,i] 
}
co_suit_fut[is.na(co_suit_fut)] = 0

#Set the reproduction rate: 
l1_a = 5

#Constant lambda across time
lambda_i = l1_tmp *l1_a

#Turn this into a projection across time
# lambda_i = array( l1_tmp, dim = c(ny,nx,nspp,tend) )
# rm(l1_tmp) #Just get rid of this since it is a large object. 

#Create the population matrix with the same times: 
Nt = array( 0, dim = c(ny,nx,nspp,tend )) 

#Initial conditions (applied to every spatial point):
Nt[,,,1] = array(0.1, dim = c(ny,nx,nspp) )

#Added Nt_base to provide a baseline. This is the density without any competition. 
Nt_base = Nt

### For now, assume that competition and survival are independent of space
#Competition: 
#This creates a square matrix for aij with intraspecific competition = 1 on 
#the diagonal. 
aii1 = 1 
aij1 = 0.5
aij = matrix ( aij1, nspp, nspp )
diag(aij) = matrix(aii1,nspp,1)

#Survival (should be less than 1)
s_i = c(matrix(0.5,nspp,1))

#Loop over all of time: 
for (t in 1:(tend-1) ) {
  
  #Total competition at each for lottery model
  tot_comp = (lambda_i*Nt[,,,t])
  tot_comp = rowSums(tot_comp, dim=2) 
  
  #Loop over all species: 
  for (s in 1:nspp) {
    
    #Multiply each competitor by aij using sweep( )
    intra_comp = sweep( Nt[, , -s, t], 3, aij[s, (-s)], FUN = "*")
    #A handy trick with rowSums to sum over 3rd dimension
    intra_comp = rowSums(intra_comp, dim=2) 
    
    Nt[,,s,t+1] = ( lambda_i[,,s] / (1+ aij[s,s] *Nt[,,s,t]+ intra_comp) + 
                      s_i[s] ) * Nt[,,s,t]
    #Checks on NAs and negatives: 
    Nt[,,s,t+1][is.na(Nt[,,s,t+1])] = 0
    
    #The model with no competition
    Nt_base [,,s,t+1] =  ( lambda_i[,,s]/(1+ aij[s,s] *Nt_base[,,s,t]) + s_i[s] ) * Nt_base[,,s,t]
    Nt_base[,,s,t+1][is.na(Nt_base[,,s,t+1])] = 0
    
    # Nt[,,s,t+1] = ( lambda_i[,,s] / tot_comp + s_i[s] ) * Nt[,,s,t]
    # Nt[,,s,t+1][is.na(Nt[,,s,t+1])] = 0
    
  }
  
}


#Nt_base - no competition
#Nt - competition
max(Nt_base[,,,tend])
max(Nt[,,,tend])

#change in population between no competition and competition
#delta_pop = (daphC_Nt_base - Nt_base_pop)/daphC_Nt_base
delta_pop_fut = (Nt_base[,,,tend] - Nt[,,,tend])/Nt_base[,,,tend]
dim(delta_pop_fut)

#No competition - equivalent to maxent suitability
Nt_base_pop_fut = Nt_base[,,,tend]
dim(Nt_base_pop_fut)

#Interspecific competition
Nt_pop_fut = Nt[,,,tend]
dim(Nt_pop_fut)

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

suit_final_fut = Nt_suit[,,,tend]
suit_change_fut = suit_final -l1_tmp 

###Convert these to rasters with same qualities as original rasters:
daph_suit_final_fut = vector("list",nspp)
daph_suit_change_fut = vector("list",nspp)
daphC_delta_pop_fut = vector("list", nspp)
daphC_Nt_base_fut=vector("list", nspp)
daphC_Nt_fut=vector("list", nspp)

for ( i in 1:nspp) {
  
  daph_suit_final_fut[[i]] =daph_pred_future[[i]]
  daph_suit_change_fut[[i]] = daph_pred_future[[i]]
  daph_cosuit_fut[[i]] =daph_pred_future[[i]]
  daphC_delta_pop_fut[[i]]=daph_pred_future[[i]]
  daphC_Nt_base_fut[[i]]=daph_pred_future[[i]]
  daphC_Nt_fut[[i]]=daph_pred_future[[i]]
  
  daph_suit_final_fut[[i]]@data = raster(suit_final_fut[,,i])@data
  daph_suit_change_fut[[i]]@data = raster(suit_change_fut[,,i])@data
  daph_cosuit_fut[[i]]@data = raster(log(10^3*co_suit_fut+1))@data
  daphC_delta_pop_fut[[i]]@data=raster(delta_pop_fut[,,i])@data
  daphC_Nt_base_fut[[i]]@data=raster(Nt_base_pop_fut[,,i])@data
  daphC_Nt_fut[[i]]@data=raster(Nt_pop_fut[,,i])@data
  
  
}

#==============================================================================
# Save these as raster files for GIS.
#==============================================================================

for (n in 1:nspp){ 
  output_suit_final_fut = paste("./../QGIS/",d_species[n],"_suit_final_fut.tif",sep="")
  writeRaster(daphC_suit_final_fut[[i]], filename = output_suit_final_fut, format = "GTiff", overwrite = TRUE)

  output_suit_change_fut = paste("./../QGIS/",d_species[n],"_suit_change_fut.tif",sep="")
  writeRaster(daphC_suit_change_fut[[i]] , filename = output_suit_change_fut, format = "GTiff", overwrite = TRUE)

  output_delta_pop_fut = paste("./../QGIS/",d_species[n],"_delta_pop_fut.tif",sep="")
  writeRaster(daphC_delta_pop_fut[[i]], filename =  output_delta_pop_fut, format = "GTiff", overwrite = TRUE)

  output_Nt_base_fut = paste("./../QGIS/",d_species[n],"_suit_final_fut.tif",sep="")
  writeRaster(daphC_Nt_base_fut[[i]], filename = output_Nt_base_fut, format = "GTiff", overwrite = TRUE)

  output_NT_fut = paste("./../QGIS/",d_species[n],"_NT_fut.tif",sep="")
  writeRaster(daphC_Nt_fut[[i]], filename =  output_NT_fut , format = "GTiff", overwrite = TRUE)
}

#==============================================================================
# Compare and Plot population density-suitability .
#==============================================================================
#MaxEnt suitability: 
fig.name = paste("daphnia_maxent_suit_fut.pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  
  plot(daph_pred_future[[i]], main="Predicted Future Suitability")
  plot(wrld_simpl, fill=FALSE, add=TRUE)
  points(daph_present[[i]]$LONGITUDE, daph_present[[i]]$LATITUDE, pch="+", cex=0.2)
  
}

dev.off()

#MaxEnt co-suitability: 
fig.name = paste("daphnia_maxent_cosuit_fut.pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  
  plot(daph_cosuit_fut[[i]], main="Predicted Future Suitability")
  plot(wrld_simpl, fill=FALSE, add=TRUE)
  points(daph_present[[i]]$LONGITUDE, daph_present[[i]]$LATITUDE, pch="+", cex=0.2)
  
}

dev.off()

#Final population-driven suitability
fig.name = paste("daphnia_pop_driven_suit_fut.pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  
  plot(daph_suit_final_fut[[i]], main="Predicted Future Suitability")
  plot(wrld_simpl, fill=FALSE, add=TRUE)
  points(daph_present[[i]]$LONGITUDE, daph_present[[i]]$LATITUDE, pch="+", cex=0.2)
  
}

dev.off()


#Final change in suitability
fig.name = paste("daphnia_pop_driven_change_fut.pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  
  plot(daph_suit_change_fut[[i]], main="Predicted Suitability")
  plot(wrld_simpl, fill=FALSE, add=TRUE)
  points(daph_present[[i]]$LONGITUDE, daph_present[[i]]$LATITUDE, pch="+", cex=0.2)
  
}

dev.off()

#==============================================================================
###3.0 Difference in habitat suitability between current and future (2070) conditions
#==============================================================================
#Change between future and present maxEnt suitability: 
suit_change_pre_fut<-daphC_pred_future[[i]]-daphC_pred[[i]]
fig.name = paste("daphnia_maxent_suit_change_pre_fut.pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  
  plot(suit_change_pre_fut[[i]], main="Predicted Change in Suitability")
  plot(wrld_simpl, fill=FALSE, add=TRUE)
  points(daph_sppC[[i]]$lon, daph_sppC[[i]]$lat, pch="+", cex=0.2)
  
}

dev.off()


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

#drCChangePoints <- extract(d rC_change, drC_coordinates)
#hist(drCChangePoints, main ="")
#abline(v=0, col="red")

dlongCChangePoints <- extract(dlongC_change, dlongC_coordinates)
hist(dlongCChangePoints, main ="")
abline(v=0, col="red")

###################################################################################
####END OF CODE