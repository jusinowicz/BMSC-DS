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

#################################Presence only data
options(java.parameters = "-Xmx4g")

library(raster)
library(rgdal)
library(maps)
library(mapdata)
library(dismo)  
library(rJava)  
library(maptools)
library(jsonlite)
library(rgbif) 
library(proxy)
library(bcmaps)
library(corrplot)

#==============================================================================
# Choose which species to run: 
#==============================================================================

d_genus = "daphnia"
d_species = c("pulex*", 'longiremis*','pulicaria*')#, 'rosea*' )
# 'lumholtzi*', 'middendorffiana*','magna*', 
nspp = length(d_species)
spp_gbif = vector("list", nspp) #Vector of lists to store gbif data

#==============================================================================
### Data Access
#==============================================================================
#Current environment from worldclim
#Note that you have to set download=T if you haven't downloaded the data before:
currentEnv=raster::getData('worldclim', var="bio", res=2.5, download=T)

#Future environment based on the 8.5 concentration scenario for 2070 from HADGEM2-ES model
futureEnv=raster::getData('CMIP5', var='bio', res=2.5, rcp=85, model='HE', year=70, download=T)
names(futureEnv)=names(currentEnv)

#World map:
data(wrld_simpl)

#Species locations from gbif - for example, the timber rattlesnake
for (i in 1:nspp ) {

	spp_gbif[[i]] = gbif(paste(d_genus), paste(d_species[i]) )

}

#==============================================================================
### Create map of BC
#==============================================================================
 provinces <- "British Columbia"
  canada <- raster::getData("GADM",country="CAN",level=1)
  bc <- canada[canada$NAME_1 %in% provinces,]
  bc_bbox <- bbox(bc)
  xlim_bc <- c(min(bc_bbox[1,1]),max(bc_bbox[1,2]))
  ylim_bc <- c(min(bc_bbox[2,1]),max(bc_bbox[2,2]))
 #plot(bc, xlim=xlim_bc, ylim=ylim_bc)

#==============================================================================
### Processing GBIF points
#==============================================================================
daph_spp = vector("list",nspp) #Vector of lists to store gbif data


# get rid of occurences without location information
for (i in 1:nspp ) {

	daph_spp[[i]] = subset(spp_gbif[[i]], !is.na(lon) & !is.na(lat) & !is.na(year))

	# find and eliminate duplicate locations
	dpdups = duplicated(daph_spp[[i]][, c("lon", "lat")])
	daph_spp[[i]]  = daph_spp[[i]] [!dpdups, ] 

}
  
#Daphnia species records before 1969
  for (i in 1:nspp ) {
    
    daph_spp[[i]] = subset(daph_spp[[i]], year >= 1970)
  
    
  }

#==============================================================================
### Look at climate data
#==============================================================================

# Remove precipitation bioclim variables
currentEnv=dropLayer(currentEnv, c("bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19"))
futureEnv=dropLayer(futureEnv, c("bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19"))

current_col <- cor(getValues(currentEnv), use = "complete.obs")
#corrplot(current_col, method="ellipse")
#ellipse going from lower left to upper right is positive correlation
#ellipse going from bottom right to upper left is negative correlation
#width indicates strength of correlation
#straight line = perfect correlation
corrplot(current_col, order="AOE", method="color", addCoef.col="gray")

# Remove bioclim variables highest collinearlity 
currentEnv=dropLayer(currentEnv, c("bio10", "bio11", "bio6"))
futureEnv=dropLayer(futureEnv, c("bio10", "bio11", "bio6"))

#==============================================================================
###1.3 Initial plots
#==============================================================================
# make initial plots (worldwide) for diagnostic purposes
col_use = c("red","blue", "green", "purple")
legend_use = c("D.pulex", "D.longiremis", "D.pulicaria", "D.rosea" )

plot(wrld_simpl, xlim=c(min(daph_spp[[1]]$lon)-1,max(daph_spp[[1]]$lon)+1), 
		ylim=c(min(daph_spp[[1]]$lat)-1,max(daph_spp[[1]]$lat)+1), axes=TRUE, col="light yellow")
for (i in 1:nspp ) {

		points(daph_spp[[i]]$lon, daph_spp[[i]]$lat, col=col_use[i], pch=20, cex=0.75)

}

legend("bottomleft", 
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
daph_sppC = vector("list",nspp) #Vector of lists to store data
for (i in 1:nspp ) {
	daph_sppC[[i]] = subset(daph_spp[[i]], country == "Canada")
}

# make initial plot (Canada) for diagnostic purposes
plot(wrld_simpl, xlim=c(min(daph_sppC[[i]]$lon)-10,max(daph_sppC[[i]]$lon)+25), 
     ylim=c(min(daph_sppC[[i]]$lat)-2,max(daph_sppC[[i]]$lat)+10), axes=TRUE, col="light yellow")

for (i in 1:nspp ) {
  
  points(daph_sppC[[i]]$lon, daph_sppC[[i]]$lat, col=col_use[i], pch=20, cex=0.75)
  
}

legend("topright", 
       legend = legend_use, text.font=c(3),
       col = col_use, 
       pch = c(19,19), 
       bty = "n", 
       cex = 1.2, 
       text.col = "black", 
       inset = c(0.03, 0.01))

#==============================================================================
### Focus in on BC: 
#==============================================================================
 bc_daph = vector("list", nspp)
 for (i in 1:nspp ) {
    bc_daph[[i]]=subset(daph_sppC[[i]], adm1== "British Columbia")
 
  }
 
  plot(bc, col="light yellow")
  for (i in 1:nspp ) {
 
    points(bc_daph[[i]]$lon, bc_daph[[i]]$lat, col=col_use[i], pch=20, cex=0.75)
 
  }
 
  legend("bottomleft",
         legend = legend_use, text.font=c(3),
         col = col_use,
         pch = c(19,19),
         bty = "n",
         cex = 1.2,
         text.col = "black",
         inset = c(0.01, 0.01))
  
  #d. magna and d. middendorffiana have less than 5 records


# old_daph = vector("list", nspp)
#  for (i in 1:nspp ) {
#    
#    old_daph[[i]] = subset(bc_daph[[i]], year < 1970)
#    
#  }


#==============================================================================
### Create MaxEnt model
#==============================================================================
set.seed(1234)
model.extent = vector("list",nspp) #Vector of lists to store data
modelEnv = vector("list",nspp) #Vector of lists to store data
modelFutureEnv = vector("list",nspp) #Vector of lists to store data

x_lims = c(-140.1667,-47.58333)
y_lims = c(33.75 ,69.125 )

for (i in 1:nspp ) {

	model.extent[[i]] = extent(min(daph_sppC[[i]]$lon)-10,max(daph_sppC[[i]]$lon)+10,
		min(daph_sppC[[i]]$lat)-10,max(daph_sppC[[i]]$lat)+10)
	modelEnv[[i]] = crop(currentEnv, model.extent[[i]])
	#modelFutureEnv[[i]] =  crop(futureEnv, model.extent[[i]])

	#Get the overall min/max extents. This is to create one common extent 
	x_lims = c( min( c( x_lims[1], xmin(model.extent[[i]] ) ) ) , 
				max(c( x_lims[2], xmax(model.extent[[i]] ) ) )
				)
	y_lims = c( min( c( y_lims[1], ymin(model.extent[[i]] ) ) ) , 
			max(c( y_lims[2], ymax(model.extent[[i]] ) ) )
			)

}

#==============================================================================
### Model assessment
#==============================================================================
set.seed(1234)

daphC_coordinates = vector("list",nspp) #Vector of lists to store data
fold = vector("list",nspp) #Vector of lists to store data
daphC_test = vector("list",nspp) #Vector of lists to store data
daphC_train = vector("list",nspp) #Vector of lists to store data

for (i in 1:nspp ) {

	#first, just make a data frame of latitudes and longitudes for the model
	daphC_coordinates[[i]] = cbind.data.frame(daph_sppC[[i]]$lon,daph_sppC[[i]]$lat)
	# add an index that makes five random groups of observations
	fold[[i]] = kfold(daphC_coordinates[[i]], k=5) 
	# hold out one fifth as test data
	daphC_test[[i]] = daphC_coordinates[[i]][fold[[i]] == 1, ] 
	# the other four fifths are training data
	daphC_train[[i]] = daphC_coordinates[[i]][fold[[i]] != 1, ] 
}

#==============================================================================
### Fit SDM using Maxent algorithm
#============================================================================== 
set.seed(1234)

daphC_model = vector("list",nspp)

par (mfrow = c(ceiling(nspp/2),2))
for (i in 1:nspp ) {

	daphC_model[[i]] = maxent(modelEnv[[i]], daphC_train[[i]])
	plot(daphC_model[[i]],
	xlab="Variable contribution (%)", 
	ylab="Bioclimatic variables", 
	main=" ")

}

#==============================================================================
### Examine sensitivity to climatic variables
#============================================================================== 
set.seed(1234)
#How does the likelihood of species occurrence respond to variation 
#in these climatic conditions?
response (daphC_model[[i]])

##############################################################################
#CURRENT SUITABILITY
#==============================================================================
### Predict current suitability
#============================================================================== 
set.seed(1234)
daphC_pred = vector("list",nspp)

par (mfrow = c(ceiling(nspp/2),2))
for (i in 1:nspp){

	#Note: I'm projecting these to the same modelEnv to get a common
	#resolution and extent. 
	daphC_pred[[i]]  = predict(daphC_model[[i]], modelEnv[[1]])
	plot(daphC_pred[[i]], main="Predicted Suitability")
	plot(wrld_simpl, fill=FALSE, add=TRUE)
	points(daph_sppC[[i]]$lon, daph_sppC[[i]]$lat, pch="+", cex=0.2)

	#Project current suitability to the same extent: 
	# extent( daphC_pred[[i]]) = c(x_lims,y_lims)
	# res( daphC_pred[[i]] ) = res (daphC_pred[[1]])

}

#==============================================================================
### Create pseudoabsences and model accuracy
#==============================================================================
set.seed(1234)
bg = vector("list",nspp)
ev_bg = vector("list",nspp)

par (mfrow = c(ceiling(nspp),2))
for (i in 1:nspp) {

	#Generate background points for pseudoabsences
	bg[[i]] = randomPoints(modelEnv[[i]], 1000)
	ev_bg[[i]] = evaluate (daphC_model[[i]], p = daphC_test[[i]], 
		a = bg [[i]], x = modelEnv[[i]] )
	plot (ev_bg[[i]], 'ROC' )

}

#==============================================================================
### Apply population dynamics (including competition), using the 
#	 maxent suitability maps to describe species intrinsic reproductive fitness
#==============================================================================
#=============================================================================
### Lottery model
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
co_suit = matrix(0,ny,nx)
#x = l1_tmp[,,i] - co_suit

for ( i in 1:nspp) {

	l1_tmp[,,i] = raster::as.matrix(daphC_pred[[i]]) #maxent suitability
	co_suit=proxy::dist(rbind(l1_tmp[,,i], 0))
	#co_suit = stats::dist(rbind(co_suit, method="euclidean"))
	#calculating euclidean distance between maxent suitability and co-occurrance 
	#co_suit = proxy::dist(rbind(co_suit, 0))
	#calculation eclidean distance between 0 (origin) and co-occurrance I think
	#co_suit = co_suit*l1_tmp[,,i]
}
co_suit[is.na(co_suit)] = 0

#Set the reproduction rate: 
l1_a = 5

#Constant lambda across time
#no variation in time
#equivalent to converting suitability to a per-capita reproduction
lambda_i = l1_tmp *l1_a

#Turn this into a projection across time
#lambda_i = array( l1_tmp, dim = c(ny,nx,nspp,tend) )
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
aij1 = 1
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
		Nt_base [,,s,t+1] =  ( lambda_i[,,s]/(1+ aij[s,s] *Nt[,,s,t]) + s_i[s] ) * Nt_base[,,s,t]
		Nt_base[,,s,t+1][is.na(Nt_base[,,s,t+1])] = 0

		# Nt[,,s,t+1] = ( lambda_i[,,s] / tot_comp + s_i[s] ) * Nt[,,s,t]
		# Nt[,,s,t+1][is.na(Nt[,,s,t+1])] = 0

	}

}

###Total population 
# pop_tot = matrix(0,tend,nspp)
# for (i in 1:nspp) {
# 	for (t in 1:tend) {	
# 		pop_tot[t,i] = sum( colSums (Nt[,,i,t]) )
# 	}
# }

delta_pop = (lambda_i-Nt[,,,tend])/lambda_i


###Convert suitability of maxent to population metric
#to only see competition effect
# Nt_suit = Nt #population suitability dataframe
# pop_max = matrix(0,tend,nspp)
# for (i in 1:nspp) {
# 	for (t in 1:tend) {	
# 	  Nt_suit = (lambda_i - Nt)/lambda_i
# 		pop_max[t,i] = max( max (Nt[,,i,t]) )
# 		Nt_suit[,,i,t] = Nt[,,i,t] / pop_max[t,i] #gets you population density
# 		Nt_suit[,,i,t] =  lambda_i/ Nt
# 	}
# }

#suit_final = Nt_suit[,,,tend] #population-driven suitability 
#suit_change = suit_final - l1_tmp 
#suit_change = lambda_i - 	Nt_suit #left with population suitability



###Convert these to rasters with same qualities as original rasters:
# daphC_suit_final = vector("list",nspp)
# daphC_suit_change = vector("list",nspp)
# daphC_cosuit = vector("list",nspp)
# daphC_comp = vector("list",nspp)
daphC_change = vector("list", nspp)


for ( i in 1:nspp) {

	# daphC_suit_final[[i]] =daphC_pred[[i]]
	# daphC_suit_change[[i]] = daphC_pred[[i]]
	# daphC_cosuit[[i]] =daphC_pred[[i]]
	#daphC_comp[[i]] = daphC_pred[[i]]
	daphC_change[[i]]=daphC_pred[[i]]

	# daphC_suit_final[[i]]@data = raster(suit_final[,,i])@data
	# daphC_suit_change[[i]]@data = raster(suit_change[,,i])@data
	# #daphC_cosuit[[i]]@data = raster::raster(co_suit)@data
	# daphC_cosuit[[i]]@data = raster(log(10^3*co_suit+1))@data
	# #daphC_comp[[i]]@data = raster(suit_change[,,i])@data
	daphC_change[[i]]@data=raster(delta_pop[,,i])@data
}
#==============================================================================
### Compare and Plot population density-suitability .
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

#MaxEnt co-suitability: co-occurrance
# fig.name = paste("daphnia_maxent_cosuit.pdf",sep="")
# pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)
# 
# par (mfrow = c(1,ceiling(nspp)))
# for (i in 1:nspp){
# 
# 	plot(daphC_cosuit[[i]], main="Predicted Suitability")
# 	plot(wrld_simpl, fill=FALSE, add=TRUE)
# 	points(daph_sppC[[i]]$lon, daph_sppC[[i]]$lat, pch="+", cex=0.2)
# 
# }
# 
# dev.off()

#Final change in suitability
fig.name = paste("daphnia_suit_change.pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

par (mfrow = c(1,ceiling(nspp)))
 for (i in 1:nspp){
 
 	plot(daphC_change[[i]], main="Predicted Suitability")
 	plot(wrld_simpl, fill=FALSE, add=TRUE)
 	points(daph_sppC[[i]]$lon, daph_sppC[[i]]$lat, pch="+", cex=0.2)
 
 }
 
 dev.off()


daphC_change[[i]]@data=raster(delta_pop[,,i])@data



#Final population-driven suitability
# fig.name = paste("daphnia_pop_driven_suit.pdf",sep="")
# pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)
# 
# par (mfrow = c(1,ceiling(nspp)))
# for (i in 1:nspp){
# 
# 	plot(daphC_suit_final[[i]], main="Predicted Suitability")
# 	plot(wrld_simpl, fill=FALSE, add=TRUE)
# 	points(daph_sppC[[i]]$lon, daph_sppC[[i]]$lat, pch="+", cex=0.2)
# 
# }
# 
# dev.off()
# 
# 
# #Final change in suitability
# fig.name = paste("daphnia_pop_driven_change.pdf",sep="")
# pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)
# 
# par (mfrow = c(1,ceiling(nspp)))
# for (i in 1:nspp){
# 
# 	plot(daphC_suit_change[[i]], main="Predicted Suitability")
# 	plot(wrld_simpl, fill=FALSE, add=TRUE)
# 	points(daph_sppC[[i]]$lon, daph_sppC[[i]]$lat, pch="+", cex=0.2)
# 
# }
# 
# dev.off()

#==============================================================================
### Compare difference in population-driven suitability and maxent suitability
#==============================================================================

maxent_suit = vector("list",nspp)
par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  #daphC_change[[i]]= extract(daphC_change[[i]], daphC_coordinates[[i]])
  maxent_suit[[i]]=hist(daphC_pred[[i]], breaks=seq(from = 0, to = 1, by = 0.1))
}

par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  #daphC_change[[i]]= extract(daphC_change[[i]], daphC_coordinates[[i]])
  daphC_change[[i]]=hist(daphC_change[[i]], breaks=seq(from = 0, to = 1, by = 0.1))
}

# cell_area=maxent_suit[[i]]$counts*

ar = vector("list",nspp)
for (i in 1:nspp){
ar[[i]] =  raster::area(daphC_pred[[1]], na.rm=TRUE, weights=FALSE)
}

# cell_size<-area(daphC_pred, na.rm=TRUE, weights=FALSE)
# ?area


#Suitability change between population-drive and maxent suitability: 

par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){

	plot(daphC_suit_change[[i]], main="Predicted change in suitability")
	plot(wrld_simpl, fill=FALSE, add=TRUE)
	points(daph_sppC[[i]]$lon, daph_sppC[[i]]$lat, pch="+", cex=0.2)

}

#Visualize change in a histogram
# change_points = vector("list",nspp)
# 
# par (mfrow = c(1,ceiling(nspp)))
#   for (i in 1:nspp){
# 
#   change_points[[i]] = extract(daphC_suit_change[[i]], daphC_coordinates[[i]])
#   hist(change_points[[i]], main="")
#   #lines(change_points[[i]], main="")
#   abline(v=0, col="red")
#   }


combined = vector("list",nspp)

for (i in 1:nspp){
daphC_pred[[i]]=as.vector(daphC_pred[[i]])
}

for (i in 1:nspp){
  daphC_suit_final[[i]]=as.vector(daphC_suit_final[[i]])
}


for (i in 1:nspp){
combined=rbind(daphC_pred[[i]], daphC_suit_final[[i]])
}

par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
    hist(combined, main="")
  abline(v=0, col="red")
}

#########################################################################################
#Calculating area of most suitable raster habitat

#### Start with maxent present suitability
#change raster to dataframe
daphC_pred_present = vector("list",nspp)
for (i in 1:nspp){
daphC_pred_present[[i]] = as.data.frame(daphC_pred[[i]], xy=TRUE)
}
#get rid of NAs
for (i in 1:nspp){
  daphC_pred_present[[i]] = na.omit(daphC_pred_present[[i]])
}



#keep cells that have a suitability score between 0-0.1 (scores range from 0 to 1)
daphC_pred_present_1 = vector("list",nspp)
for (i in 1:nspp){
  daphC_pred_present_1[[i]] = daphC_pred_present[[i]][daphC_pred_present[[i]]$layer<=0.1,]
  #re-rasterize just the suitable area
  daphC_pred_present_1[[i]]= rasterFromXYZ(daphC_pred_present_1[[i]])
  #Get sizes of all cells in current distribution raster
    daphC_pred_present_1[[i]]=area(daphC_pred_present_1[[i]], na.rm=TRUE, weights=FALSE)
  #delete NAs from all raster cells. NAs come back when going from dataframe to raster??
    daphC_pred_present_1[[i]]=daphC_pred_present_1[[i]][!is.na(daphC_pred_present_1[[i]])]
  #compute area [km2] of all cells in raster
    daphC_pred_present_1[[i]]=length(daphC_pred_present_1[[i]])*median(daphC_pred_present_1[[i]])
  }

#keep cells that have a suitability score between 0.1-0.2 (scores range from 0 to 1)
daphC_pred_present_2 = vector("list",nspp)
for (i in 1:nspp){
  daphC_pred_present_2[[i]] = daphC_pred_present[[i]][daphC_pred_present[[i]]$layer>0.1,]
  daphC_pred_present_2[[i]] = daphC_pred_present[[i]][daphC_pred_present[[i]]$layer<=0.2,]
  #re-rasterize just the suitable area
  daphC_pred_present_2[[i]]= rasterFromXYZ(daphC_pred_present_2[[i]])
  #Get sizes of all cells in current distribution raster
  daphC_pred_present_2[[i]]=area(daphC_pred_present_2[[i]], na.rm=TRUE, weights=FALSE)
  #delete NAs from all raster cells. NAs come back when going from dataframe to raster??
  daphC_pred_present_2[[i]]=daphC_pred_present_2[[i]][!is.na(daphC_pred_present_2[[i]])]
  #compute area [km2] of all cells in raster
  daphC_pred_present_2[[i]]=length(daphC_pred_present_2[[i]])*median(daphC_pred_present_2[[i]])
  
}

#keep cells that have a suitability score between 0.2-0.3 (scores range from 0 to 1)
daphC_pred_present_3 = vector("list",nspp)
for (i in 1:nspp){
  daphC_pred_present_3[[i]] = daphC_pred_present[[i]][daphC_pred_present[[i]]$layer>0.2,]
  daphC_pred_present_3[[i]] = daphC_pred_present[[i]][daphC_pred_present[[i]]$layer<=0.3,]
  #re-rasterize just the suitable area
  daphC_pred_present_3[[i]]= rasterFromXYZ(daphC_pred_present_3[[i]])
  #Get sizes of all cells in current distribution raster
  daphC_pred_present_3[[i]]=area(daphC_pred_present_3[[i]], na.rm=TRUE, weights=FALSE)
  #delete NAs from all raster cells. NAs come back when going from dataframe to raster??
  daphC_pred_present_3[[i]]=daphC_pred_present_3[[i]][!is.na(daphC_pred_present_3[[i]])]
  #compute area [km2] of all cells in raster
  daphC_pred_present_3[[i]]=length(daphC_pred_present_3[[i]])*median(daphC_pred_present_3[[i]])
  
}

#keep cells that have a suitability score between 0.3-0.4 (scores range from 0 to 1)
daphC_pred_present_4 = vector("list",nspp)
for (i in 1:nspp){
  daphC_pred_present_4[[i]] = daphC_pred_present[[i]][daphC_pred_present[[i]]$layer>0.3,]
  daphC_pred_present_4[[i]] = daphC_pred_present[[i]][daphC_pred_present[[i]]$layer<=0.4,]
  #re-rasterize just the suitable area
  daphC_pred_present_4[[i]]= rasterFromXYZ(daphC_pred_present_4[[i]])
  #Get sizes of all cells in current distribution raster
  daphC_pred_present_4[[i]]=area(daphC_pred_present_4[[i]], na.rm=TRUE, weights=FALSE)
  #delete NAs from all raster cells. NAs come back when going from dataframe to raster??
  daphC_pred_present_4[[i]]=daphC_pred_present_4[[i]][!is.na(daphC_pred_present_4[[i]])]
  #compute area [km2] of all cells in raster
  daphC_pred_present_4[[i]]=length(daphC_pred_present_4[[i]])*median(daphC_pred_present_4[[i]])
  
}

#keep cells that have a suitability score between 0.4-0.5 (scores range from 0 to 1)
daphC_pred_present_5 = vector("list",nspp)
for (i in 1:nspp){
  daphC_pred_present_5[[i]] = daphC_pred_present[[i]][daphC_pred_present[[i]]$layer>0.4,]
  daphC_pred_present_5[[i]] = daphC_pred_present[[i]][daphC_pred_present[[i]]$layer<=0.5,]
  #re-rasterize just the suitable area
  daphC_pred_present_5[[i]]= rasterFromXYZ(daphC_pred_present_5[[i]])
  #Get sizes of all cells in current distribution raster
  daphC_pred_present_5[[i]]=area(daphC_pred_present_5[[i]], na.rm=TRUE, weights=FALSE)
  #delete NAs from all raster cells. NAs come back when going from dataframe to raster??
  daphC_pred_present_5[[i]]=daphC_pred_present_5[[i]][!is.na(daphC_pred_present_5[[i]])]
  #compute area [km2] of all cells in raster
  daphC_pred_present_5[[i]]=length(daphC_pred_present_5[[i]])*median(daphC_pred_present_5[[i]])
  
}

#keep cells that have a suitability score between 0.5-0.6 (scores range from 0 to 1)
daphC_pred_present_6 = vector("list",nspp)
for (i in 1:nspp){
  daphC_pred_present_6[[i]] = daphC_pred_present[[i]][daphC_pred_present[[i]]$layer>0.5,]
  daphC_pred_present_6[[i]] = daphC_pred_present[[i]][daphC_pred_present[[i]]$layer<=0.6,]
  #re-rasterize just the suitable area
  daphC_pred_present_6[[i]]= rasterFromXYZ(daphC_pred_present_6[[i]])
  #Get sizes of all cells in current distribution raster
  daphC_pred_present_6[[i]]=area(daphC_pred_present_6[[i]], na.rm=TRUE, weights=FALSE)
  #delete NAs from all raster cells. NAs come back when going from dataframe to raster??
  daphC_pred_present_6[[i]]=daphC_pred_present_6[[i]][!is.na(daphC_pred_present_6[[i]])]
  #compute area [km2] of all cells in raster
  daphC_pred_present_6[[i]]=length(daphC_pred_present_6[[i]])*median(daphC_pred_present_6[[i]])
  
}

#keep cells that have a suitability score between 0.6-0.7 (scores range from 0 to 1)
daphC_pred_present_7 = vector("list",nspp)
for (i in 1:nspp){
  daphC_pred_present_7[[i]] = daphC_pred_present[[i]][daphC_pred_present[[i]]$layer>0.6,]
  daphC_pred_present_7[[i]] = daphC_pred_present[[i]][daphC_pred_present[[i]]$layer<=0.7,]
  #re-rasterize just the suitable area
  daphC_pred_present_7[[i]]= rasterFromXYZ(daphC_pred_present_7[[i]])
  #Get sizes of all cells in current distribution raster
  daphC_pred_present_7[[i]]=area(daphC_pred_present_7[[i]], na.rm=TRUE, weights=FALSE)
  #delete NAs from all raster cells. NAs come back when going from dataframe to raster??
  daphC_pred_present_7[[i]]=daphC_pred_present_7[[i]][!is.na(daphC_pred_present_7[[i]])]
  #compute area [km2] of all cells in raster
  daphC_pred_present_7[[i]]=length(daphC_pred_present_7[[i]])*median(daphC_pred_present_7[[i]])
  
}

#keep cells that have a suitability score between 0.7-0.8 (scores range from 0 to 1)
daphC_pred_present_8 = vector("list",nspp)
for (i in 1:nspp){
  daphC_pred_present_8[[i]] = daphC_pred_present[[i]][daphC_pred_present[[i]]$layer>0.7,]
  daphC_pred_present_8[[i]] = daphC_pred_present[[i]][daphC_pred_present[[i]]$layer<=0.8,]
  #re-rasterize just the suitable area
  daphC_pred_present_8[[i]]= rasterFromXYZ(daphC_pred_present_8[[i]])
  #Get sizes of all cells in current distribution raster
  daphC_pred_present_8[[i]]=area(daphC_pred_present_8[[i]], na.rm=TRUE, weights=FALSE)
  #delete NAs from all raster cells. NAs come back when going from dataframe to raster??
  daphC_pred_present_8[[i]]=daphC_pred_present_8[[i]][!is.na(daphC_pred_present_8[[i]])]
  #compute area [km2] of all cells in raster
  daphC_pred_present_8[[i]]=length(daphC_pred_present_8[[i]])*median(daphC_pred_present_8[[i]])
  
}


#keep cells that have a suitability score between 0.8-0.9 (scores range from 0 to 1)
daphC_pred_present_9 = vector("list",nspp)
for (i in 1:nspp){
  daphC_pred_present_9[[i]] = daphC_pred_present[[i]][daphC_pred_present[[i]]$layer>0.8,]
  daphC_pred_present_9[[i]] = daphC_pred_present[[i]][daphC_pred_present[[i]]$layer<=0.9,]
  #re-rasterize just the suitable area
  daphC_pred_present_9[[i]]= rasterFromXYZ(daphC_pred_present_9[[i]])
  #Get sizes of all cells in current distribution raster
  daphC_pred_present_9[[i]]=area(daphC_pred_present_9[[i]], na.rm=TRUE, weights=FALSE)
  #delete NAs from all raster cells. NAs come back when going from dataframe to raster??
  daphC_pred_present_9[[i]]=daphC_pred_present_9[[i]][!is.na(daphC_pred_present_9[[i]])]
  #compute area [km2] of all cells in raster
  daphC_pred_present_9[[i]]=length(daphC_pred_present_9[[i]])*median(daphC_pred_present_9[[i]])
  
}

#keep cells that have a suitability score between 0.9-1.0 (scores range from 0 to 1)
daphC_pred_present_10 = vector("list",nspp)
for (i in 1:nspp){
  daphC_pred_present_10[[i]] = daphC_pred_present[[i]][daphC_pred_present[[i]]$layer>0.9,]
  daphC_pred_present_10[[i]] = daphC_pred_present[[i]][daphC_pred_present[[i]]$layer<=1.0,]
  #re-rasterize just the suitable area
  daphC_pred_present_10[[i]]= rasterFromXYZ(daphC_pred_present_10[[i]])
  #Get sizes of all cells in current distribution raster
  daphC_pred_present_10[[i]]=area(daphC_pred_present_10[[i]], na.rm=TRUE, weights=FALSE)
  #delete NAs from all raster cells. NAs come back when going from dataframe to raster??
  daphC_pred_present_10[[i]]=daphC_pred_present_10[[i]][!is.na(daphC_pred_present_10[[i]])]
  #compute area [km2] of all cells in raster
  daphC_pred_present_10[[i]]=length(daphC_pred_present_10[[i]])*median(daphC_pred_present_10[[i]])
  
}

all = vector("list",nspp)
for (i in 1:nspp){

all[[i]] = do.call(c, list(daphC_pred_present_1, daphC_pred_present_2, daphC_pred_present_3, daphC_pred_present_4, daphC_pred_present_5, daphC_pred_present_6, daphC_pred_present_7, daphC_pred_present_8, daphC_pred_present_9, daphC_pred_present_10))
}
#Visualize change in a histogram

Habitat_suitability  = seq(from = 0, to = 1, by = 0.1)
Area=c()



par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  
  hist(all[[i]], main="")
  #lines(change_points[[i]], main="")
  abline(v=0, col="red")
}



####Do the same with population suitability
#change raster to dataframe
# daphC_suit_final_df = vector("list",nspp)
# for (i in 1:nspp){
#   daphC_suit_final_df[[i]] = as.data.frame(daphC_suit_final[[i]], xy=TRUE)
# }
# #get rid of NAs
# for (i in 1:nspp){
#   daphC_suit_final_df[[i]] = na.omit(daphC_suit_final_df[[i]])
# }
# 
# #keep cells that have a suitability score above 0.5 (scores range from 0 to 1)
# for (i in 1:nspp){
#   daphC_suit_final_df[[i]] = daphC_suit_final_df[[i]][daphC_suit_final_df[[i]]$layer>0.5,]
# }
# 
# #re-rasterize just the suitable area
# daphC_suit_final_ras = vector("list",nspp)
# for (i in 1:nspp){
#   daphC_suit_final_ras[[i]] = rasterFromXYZ(daphC_suit_final_df[[i]])
# }

####Get sizes of all cells in current distribution raster
daphC_pred_cell_size = vector("list",nspp)
for (i in 1:nspp){
  daphC_pred_cell_size[[i]]<-area(daphC_pred_ras[[i]], na.rm=TRUE, weights=FALSE)
}

#delete NAs from all raster cells. NAs come back when going from dataframe to raster??
for (i in 1:nspp){
  daphC_pred_cell_size[[i]]<-daphC_pred_cell_size[[i]][!is.na(daphC_pred_cell_size[[i]])]
}

#compute area [km2] of all cells in raster
daphC_pred_area = vector("list",nspp)
for (i in 1:nspp){
  daphC_pred_area[[i]]<-length(daphC_pred_cell_size[[i]])*median(daphC_pred_cell_size[[i]])
}

####Get sizes of all cells in population distribution raster
daphC_suit_final_cell_size=vector("list",nspp)
for (i in 1:nspp){
  daphC_suit_final_cell_size[[i]]<-area(daphC_suit_final_ras[[i]], na.rm=TRUE, weights=FALSE)
}

#delete NAs from all raster cells. NAs come back when going from dataframe to raster??
for (i in 1:nspp){
  daphC_suit_final_cell_size[[i]]<-daphC_suit_final_ras[[i]][!is.na(daphC_suit_final_ras[[i]])]
}

#compute area [km2] of all cells in raster
daphC_suit_final_area = vector("list",nspp)
for (i in 1:nspp){
  daphC_suit_final_area[[i]]<-length(daphC_suit_final_cell_size[[i]])*median(daphC_suit_final_cell_size[[i]])
}

###calculate change in area
delta_area = vector("list",nspp)
for (i in 1:nspp){
  delta_area[[i]] <- daphC_pred_area[[i]] - daphC_suit_final_area[[i]]
}


daphC_pred_df = vector("list",nspp)
for (i in 1:nspp){
  
  daphC_pred_df[[i]]=as.(daphC_pred[[i]], main="Predicted Suitability")
  plot(wrld_simpl, fill=FALSE, add=TRUE)
  points(daph_sppC[[i]]$lon, daph_sppC[[i]]$lat, pch="+", cex=0.2)
  
}






#Visualize change in a line plot
maxent_points = vector("list",nspp)
par (mfrow = c(1,ceiling(nspp)))
  for (i in 1:nspp){

  maxent_points[[i]] = extract(daphC_pred[[i]], daphC_coordinates[[i]])
  maxent_points[[i]] = subset(maxent_points[[i]], !is.na(maxent_points[[i]] ))
  }

pop_points=vector("list", nspp)
par (mfrow = c(1,ceiling(nspp)))
  for (i in 1:nspp){

  pop_points[[i]] = extract(daphC_suit_final[[i]], daphC_coordinates[[i]])
  }

#Suitability 

suit  = seq(from = 0, to = 1, by = 0.1)

suit
par (mfrow = c(1,ceiling(nspp)))
  for (i in 1:nspp){

  lines(suit, maxent_points[[i]])
  }



par (mfrow = c(1,ceiling(nspp)))
  for (i in 1:nspp){

  change_points[[i]] = extract(daphC_suit_change[[i]], daphC_coordinates[[i]])
  plot(change_points[[i]], type="1")
}

?plot

#####################################################################################
###FUTURE SUITABILITY
#==============================================================================
### Predict future (2070) suitability
#==============================================================================

daphC_pred_future = vector("list",nspp)

par (mfrow = c(ceiling(nspp/2),2))
for (i in 1:nspp){
  
  #Note: I'm projecting these to the same modelFutureEnv to get a common
  #resolution and extent. 
  daphC_pred_future[[i]]  = predict(daphC_model[[i]], modelFutureEnv[[1]])
  plot(daphC_pred_future[[i]], main="Predicted Future Suitability")
  plot(wrld_simpl, fill=FALSE, add=TRUE)
 points(daph_sppC[[i]]$lon, daph_sppC[[i]]$lat, pch="+", cex=0.2)
  
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
l1_tmp = raster::as.matrix(daphC_pred_future[[1]])
nx = ncol(l1_tmp)
ny = nrow(l1_tmp)

#Get each species' lambda: 
l1_tmp = array( matrix(0,ny,nx), dim = c(ny,nx,nspp )  )
co_suit_fut = matrix(1,ny,nx)
for ( i in 1:nspp) {
  
  l1_tmp[,,i] = raster::as.matrix(daphC_pred_future[[i]])
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

#Create the population matrix with time
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

	#Total competition at each for lottery model
	tot_comp = (lambda_i*Nt[,,,t])
	tot_comp = rowSums(tot_comp, dim=2) 

  
  #Loop over all species: 
  for (s in 1:nspp) {
    
    #Multiply each competitor by aij using sweep( )
    # intra_comp = sweep( Nt[, , -s, t], 3, aij[s, (-s)], FUN = "*")
    # #A handy trick with rowSums to sum over 3rd dimension
    # intra_comp = rowSums(intra_comp, dim=2) 
    # Nt[,,s,t+1] = ( lambda_i[,,s] / (1+ aij[s,s] *Nt[,,s,t]+ intra_comp) + 
    #                   s_i[s] ) * Nt[,,s,t]
    # Nt[,,s,t+1][is.na(Nt[,,s,t+1])] = 0

    Nt[,,s,t+1] = ( lambda_i[,,s] / tot_comp + s_i[s] ) * Nt[,,s,t]
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
### population 
Nt_suit = Nt 
pop_max = matrix(0,tend,nspp)
for (i in 1:nspp) {
  for (t in 1:tend) {	
    pop_max[t,i] = max( max (Nt[,,i,t]) )
    Nt_suit[,,i,t] =Nt[,,i,t] / pop_max[t,i] #population density
  }
}

suit_final_fut = Nt_suit[,,,tend]
suit_change_fut = suit_final -l1_tmp 

###Convert these to rasters with same qualities as original rasters:
daphC_suit_final_fut = vector("list",nspp)
daphC_suit_change_fut = vector("list",nspp)
daphC_cosuit_fut = vector("list",nspp)
daphC_comp_fut = vector("list",nspp)

for ( i in 1:nspp) {
  
  daphC_suit_final_fut[[i]] =daphC_pred_future[[i]]
  daphC_suit_change_fut[[i]] = daphC_pred_future[[i]]
  daphC_cosuit_fut[[i]] =daphC_pred_future[[i]]
  #daphC_comp_fut[[i]] = daphC_pred_future[[i]]
  
  daphC_suit_final_fut[[i]]@data = raster(suit_final_fut[,,i])@data
  daphC_suit_change_fut[[i]]@data = raster(suit_change_fut[,,i])@data
  daphC_cosuit_fut[[i]]@data = raster(log(10^3*co_suit_fut+1))@data
  #daphC_comp_fut[[i]]@data = raster(suit_change[,,i])@data
  
  
}


#==============================================================================
# Compare and Plot population density-suitability .
#==============================================================================
#MaxEnt suitability: 
fig.name = paste("daphnia_maxent_suit_fut.pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  
  plot(daphC_pred_future[[i]], main="Predicted Future Suitability")
  plot(wrld_simpl, fill=FALSE, add=TRUE)
  points(daph_sppC[[i]]$lon, daph_sppC[[i]]$lat, pch="+", cex=0.2)
  
}

dev.off()

#MaxEnt co-suitability: 
fig.name = paste("daphnia_maxent_cosuit_fut.pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  
  plot(daphC_cosuit_fut[[i]], main="Predicted Future Suitability")
  plot(wrld_simpl, fill=FALSE, add=TRUE)
  points(daph_sppC[[i]]$lon, daph_sppC[[i]]$lat, pch="+", cex=0.2)
  
}

dev.off()

#Final population-driven suitability
fig.name = paste("daphnia_pop_driven_suit_fut.pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  
  plot(daphC_suit_final_fut[[i]], main="Predicted Future Suitability")
  plot(wrld_simpl, fill=FALSE, add=TRUE)
  points(daph_sppC[[i]]$lon, daph_sppC[[i]]$lat, pch="+", cex=0.2)
  
}

dev.off()


#Final change in suitability
fig.name = paste("daphnia_pop_driven_change_fut.pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  
  plot(daphC_suit_change_fut[[i]], main="Predicted Suitability")
  plot(wrld_simpl, fill=FALSE, add=TRUE)
  points(daph_sppC[[i]]$lon, daph_sppC[[i]]$lat, pch="+", cex=0.2)
  
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