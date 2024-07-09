#==============================================================================
#Lara Calvo 
#BMSC Fall Program 2020
#Directed Studies in Marine Science - R-script
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
library(corrplot)
library(prettymapr)



#==============================================================================
# Choose which species to run: 
#==============================================================================

d_genus = "daphnia"
d_species = c("pulex*", 'longiremis*','pulicaria*', 'rosea*','lumholtzi*', 'middendorffiana*','magna*') 
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
 # provinces <- "British Columbia"
 #  canada <- raster::getData("GADM",country="CAN",level=1)
 #  bc <- canada[canada$NAME_1 %in% provinces,]
 #  bc_bbox <- bbox(bc)
 #  
 #  xlim_bc <- c(min(bc_bbox[1,1]),max(bc_bbox[1,2]))
 #  ylim_bc <- c(min(bc_bbox[2,1]),max(bc_bbox[2,2]))
 # plot(bc, xlim=xlim_bc, ylim=ylim_bc)

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
  
#Daphnia species records in/after 1970
  for (i in 1:nspp ) {
    
    daph_spp[[i]] = subset(daph_spp[[i]], year >= 1970)
  
    
  }

#==============================================================================
### Look at climate data
#==============================================================================

# Remove precipitation bioclim variables
currentEnv=dropLayer(currentEnv, c("bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19"))
futureEnv=dropLayer(futureEnv, c("bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19"))

#collinearity analysis
current_col <- cor(getValues(currentEnv), use = "complete.obs")
corrplot(current_col, order="AOE", method="color", addCoef.col="gray")

# Remove bioclim variables highest collinearity
currentEnv=dropLayer(currentEnv, c("bio10", "bio11", "bio6"))
futureEnv=dropLayer(futureEnv, c("bio10", "bio11", "bio6"))

#==============================================================================
###1.3 Initial plots
#==============================================================================
# make initial plots (worldwide) for diagnostic purposes
col_use = c("red","blue", "green", "purple", "black", "gray", "pink")
legend_use = c("D.pulex", "D.longiremis", "D.pulicaria","D.rosea", "lumholtzi", "middendorffiana","magna")

plot(wrld_simpl, xlim=c(min(daph_spp[[1]]$lon)-1,max(daph_spp[[1]]$lon)+1), 
		ylim=c(min(daph_spp[[1]]$lat)-1,max(daph_spp[[1]]$lat)+1), axes=TRUE, col="light yellow")
for (i in 1:nspp ) {

		points(daph_spp[[i]]$lon, daph_spp[[i]]$lat, col=col_use[i], pch=20, cex=0.75)

}

addscalebar()
addnortharrow(pos="bottomright")


legend("bottomleft", 
       legend = legend_use, text.font=c(3),
       col = col_use, 
       pch = c(19,19), 
       bty = "n", 
       cex = 1.2, 
       text.col = "black", 
       inset = c(0.01, 0.05))

#==============================================================================
#Focus in on North America: 
#==============================================================================
daph_sppC = vector("list",nspp) #Vector of lists to store data
for (i in 1:nspp ) {
	daph_sppC[[i]] = subset(daph_spp[[i]], country %in% c("Canada", "United States", "Mexico"))
}

# make initial plot (Canada) for diagnostic purposes
plot(wrld_simpl, xlim=c(min(daph_sppC[[i]]$lon)-10,max(daph_sppC[[i]]$lon)+25), 
     ylim=c(min(daph_sppC[[i]]$lat)-2,max(daph_sppC[[i]]$lat)+10), axes=TRUE, col="light yellow")

for (i in 1:nspp ) {
  
  points(daph_sppC[[i]]$lon, daph_sppC[[i]]$lat, col=col_use[i], pch=20, cex=0.75)
  
}

addscalebar()
addnortharrow(pos="bottomright")

legend("bottomleft", 
       legend = legend_use, text.font=c(3),
       col = col_use, 
       pch = c(19,19), 
       bty = "n", 
       cex = 1.2, 
       text.col = "black", 
       inset = c(0.03, 0.05))

#==============================================================================
### Focus in on BC: 
#==============================================================================
# bc_daph = vector("list", nspp)
#  
# for (i in 1:nspp ) {
#     bc_daph[[i]]=subset(daph_sppC[[i]], adm1== "British Columbia")
#  
#   }
#  
# plot(bc, col="light yellow")
#   for (i in 1:nspp ) {
#  
#     points(bc_daph[[i]]$lon, bc_daph[[i]]$lat, col=col_use[i], pch=20, cex=0.75)
#  
#   }
#   
#   addscalebar(padin = c(2, 0.01))
#   addnortharrow(padin = c(1, 0.2), pos="bottomright")
#  
#   legend("bottomleft",
#          legend = legend_use, text.font=c(3),
#          col = col_use,
#          pch = c(19,19),
#          bty = "n",
#          cex = 1.2,
#          text.col = "black",
#          inset = c(0.2, 0.04))


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
# set.seed(1234)
# #How does the likelihood of species occurrence respond to variation 
# #in these climatic conditions?
# response (daphC_model[[i]])

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
	addscalebar(pos="bottomleft")
	addnortharrow(scale=0.5, pos="bottomright")
	

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
### Leslie-Gower
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

for ( i in 1:nspp) {
  
  l1_tmp[,,i] = raster::as.matrix(daphC_pred[[i]]) 
}

#Set the reproduction rate: 
l1_a = 5

#Constant lambda across time - no variation in time
#equivalent to converting suitability to a per-capita reproduction
lambda_i = l1_tmp *l1_a

#Remove sites that start out at (approximately) 0 suitability
s_thresh = 0.01
lambda_i [lambda_i < s_thresh] = NA

#Create the population matrix with the same times: 
Nt = array( 0, dim = c(ny,nx,nspp,tend )) 

#Initial conditions (applied to every spatial point):
Nt[,,,1] = array(0.1, dim = c(ny,nx,nspp) )

#Added Nt_base to provide a baseline. This is the density without any competition. #just intraspecific copetition 
Nt_base = Nt

### For now, assume that competition and survival are independent of space
#Competition: 
#This creates a square matrix for aij with intraspecific competition = 1 on 
#the diagonal. 
aii1 = 1 
aij1 = 0.1
aij = matrix ( aij1, nspp, nspp )
diag(aij) = matrix(aii1,nspp,1)

#Survival (should be less than 1)
s_i = c(matrix(0.5,nspp,1))

#Loop over all of time: 
for (t in 1:(tend-1) ) {
  
  #Loop over all species: 
  for (s in 1:nspp) {
    
    #Multiply each competitor by aij using sweep( )
    intra_comp = sweep( Nt[, , -s, t], 3, aij[s, (-s)], FUN = "*")
    #A handy trick with rowSums to sum over 3rd dimension
    intra_comp = rowSums(intra_comp, dim=2,na.rm=T) 
    
    Nt[,,s,t+1] = ( lambda_i[,,s] / (1+ aij[s,s] *Nt[,,s,t]+ intra_comp) + 
                      s_i[s] ) * Nt[,,s,t]
    
    #The model with no inter-competition 
    #just species competing with itself
    Nt_base [,,s,t+1] =  ( lambda_i[,,s]/(1+ aij[s,s] *Nt_base[,,s,t]) + s_i[s] ) * Nt_base[,,s,t]
    
    
  }
  
}

memory.size()
memory.limit()
memory.limit(size=56000)

#Nt_base - no competition
#Nt - competition
max(Nt_base[,,,tend],na.rm=T)
max(Nt[,,,tend],na.rm=T)

#No competition - equivalent to maxent suitability
Nt_base_pop = Nt_base[,,,tend]
dim(Nt_base_pop)

#Interspecific competition
Nt_pop = Nt[,,,tend]
dim(Nt_pop)

#change in population between no competition and competition
#delta_pop = (daphC_Nt_base - Nt_base_pop)/daphC_Nt_base
delta_pop = (Nt_base[,,,tend] - Nt[,,,tend])/Nt_base[,,,tend]
dim(delta_pop)

###Convert these to rasters with same qualities as original rasters:
daphC_delta_pop = vector("list", nspp)
daphC_Nt_base = vector("list", nspp)
daphC_Nt = vector("list", nspp)


for ( i in 1:nspp) {

  daphC_delta_pop[[i]]=daphC_pred[[i]]
  daphC_Nt_base[[i]]=daphC_pred[[i]]
  daphC_Nt[[i]]=daphC_pred[[i]]
 
  daphC_delta_pop[[i]]@data=raster(delta_pop[,,i])@data
  daphC_Nt_base[[i]]@data=raster(Nt_base_pop[,,i])@data
  daphC_Nt[[i]]@data=raster(Nt_pop[,,i])@data
  
}






#==============================================================================
### Compare and Plot population density-suitability .
#==============================================================================
###MaxEnt suitability map: 
fig.name = paste("daphnia_maxent_suit.pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){

	plot(daphC_pred[[i]], main="Predicted Suitability")
	plot(wrld_simpl, fill=FALSE, add=TRUE)
	points(daph_sppC[[i]]$lon, daph_sppC[[i]]$lat, pch="+", cex=0.2)
	addscalebar(pos="bottomright")
	addnortharrow(scale=0.5, pos="bottomleft")
	
}

Sdev.off()

###Population suitability with no interspecific competition: equivalent to maxent 
fig.name = paste("daphnia_Nt_base_suit.pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

par (mfrow = c(3,1))
for (i in 1:nspp){
  
  plot(daphC_Nt_base[[i]], main=" ")
  plot(wrld_simpl, fill=FALSE, add=TRUE)
  points(daph_sppC[[i]]$lon, daph_sppC[[i]]$lat, pch=20, cex=0.75)
  addscalebar(pos="bottomright")
  addnortharrow(scale=0.5, pos="bottomleft")
  
}
dev.off()

#Population suitability with interspecific competition: 
fig.name = paste("daphnia_Nt_suit.pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)
par (mfrow = c(3,1))
for (i in 1:nspp){
  
  plot(daphC_Nt[[i]], main=" ")
  plot(wrld_simpl, fill=FALSE, add=TRUE)
  points(daph_sppC[[i]]$lon, daph_sppC[[i]]$lat, pch=20, cex=0.75)
  addscalebar(pos="bottomright")
  addnortharrow(scale=0.5, pos="bottomleft")
}
dev.off()


#Suitability change between competition and no competition
fig.name = paste("daphnia_delta_change.pdf",sep="")
pdf(file=fig.name, height=8, width=7, onefile=TRUE, family='Helvetica', pointsize=16)

options(repr.plot.width = 2, repr.plot.height = 8)

par (mfrow = c(3,1))
 for (i in 1:nspp){
   
 
 	plot(daphC_delta_pop[[i]], main=" ")
 	plot(wrld_simpl, fill=FALSE, add=TRUE)
 	points(daph_sppC[[i]]$lon, daph_sppC[[i]]$lat, pch=20, cex=0.75)
 	addscalebar(pos="bottomright")
 	addnortharrow(scale=0.5, pos="bottomleft")
 	
 }
 
 dev.off()


#==============================================================================
###Graphs - everything plotted by itself
#==============================================================================
###Maxent habitat suitability histogram

maxent_suit = vector("list",nspp)
par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  maxent_suit[[i]]=hist(daphC_pred[[i]], breaks=seq(from = 0, to = 1, by = 0.1), main="", ylab="Frequency", xlab="Habitat suitability", col=NULL)
}

#Get area of one cell
ar_maxent = vector("list",nspp)
for (i in 1:nspp){
  ar_maxent[[i]]=maxent_suit[[i]]$counts*res(daphC_pred[[i]])[1]^2
}

for (i in 1:nspp){
  maxent_suit[[i]]$ar_maxent=maxent_suit[[i]]$counts*res(daphC_pred[[i]])[1]^2
}


###Maxent habitat suitability barplot
bar_maxent = vector("list",nspp)
for (i in 1:nspp){
  bar_maxent[[i]]=barplot(maxent_suit[[i]]$ar_maxent, maxent_suit[[i]]$mids, ylab="Area (Km^2)", xlab="Suitability", col=NULL)
}

###Maxent habitat suitability boxplot
boxplot_maxent_suit = vector("list",nspp)
for (i in 1:nspp){
  boxplot_maxent_suit[[i]]=boxplot(maxent_suit[[i]]$ar_maxent, ylab="Area", col=NULL)
}

###Change in population dynamic habitat suitability histogram
delta_pop_suit = vector("list",nspp)
par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  delta_pop_suit[[i]]=hist(daphC_delta_pop[[i]], breaks=seq(from = 0, to = 1, by = 0.1), main="", ylab="Frequency", xlab="Suitability", col=NULL)
}

#Get area of one cell
ar_delta_pop = vector("list",nspp)
for (i in 1:nspp){
  ar_delta_pop[[i]]= delta_pop_suit[[i]]$counts*res(daphC_delta_pop[[i]])[1]^2
  
}

for (i in 1:nspp){
  delta_pop_suit[[i]]$ar_delta_pop=delta_pop_suit[[i]]$counts*res(daphC_delta_pop[[i]])[1]^2
}


###Change in population dynamic habitat suitability barplot
bar_delta_pop_suit = vector("list",nspp)
for (i in 1:nspp){
  bar_delta_pop_suit[[i]]=barplot(daphC_delta_pop[[i]], daphC_delta_pop[[i]]$mids, ylab="Area (Km^2)", xlab="Suitability", col=NULL)
}

###Change in population dynamic habitat suitability boxplot
boxplot_delta_pop_suit = vector("list",nspp)
par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  boxplot_delta_pop_suit[[i]]=boxplot(daphC_delta_pop[[i]], ylab="Suitability", col=NULL)
}


###Population dynamic habitat suitability (no competition) histogram
Nt_base_suit = vector("list",nspp)
par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  
  Nt_base_suit[[i]]=hist(daphC_Nt_base[[i]], breaks=seq(from = 0, to = 10, by = 1), main="", ylab="Frequency", xlab="Suitability (population units)", col=NULL)
}

#Get area of one cell
ar_Nt_base_pop = vector("list",nspp)
for (i in 1:nspp){
  #ar_maxent[[i]]= maxent_suit[[i]]$counts*res(daphC_pred[[i]])[1]^2
  ar_Nt_base_pop[[i]]= Nt_base_suit[[i]]$counts*res(daphC_Nt_base[[i]])[1]^2
  
}

for (i in 1:nspp){
  Nt_base_suit[[i]]$ar_Nt_base_pop=Nt_base_suit[[i]]$counts*res(daphC_Nt_base[[i]])[1]^2
}

###Population dynamic habitat suitability (no competition) barplot
bar_Nt_base_suit = vector("list",nspp)
for (i in 1:nspp){
  bar_Nt_base_suit[[i]]=barplot(Nt_base_suit[[i]]$ar_Nt_base_pop, Nt_base_suit[[i]]$mids, ylab="Area (Km^2)", xlab="Suitability", col=NULL)
}

###Population dynamic habitat suitability (no competition) boxplot
boxplot_Nt_base_suit = vector("list",nspp)
par (mfrow = c(1,ceiling(nspp)))

for (i in 1:nspp){
  boxplot_Nt_base_suit[[i]]=boxplot(daphC_Nt_base[[i]], ylab="Population density (number of individuals)", col=NULL)
}


###Population dynamic habitat suitability (competition) histogram
Nt_suit = vector("list",nspp)
par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  Nt_suit[[i]]=hist(daphC_Nt[[i]], breaks=seq(from = 0, to = 10, by = 1), main="", ylab="Frequency", xlab="Suitability (population units)", col=NULL)
}

#Get area of one cell
ar_Nt_pop = vector("list",nspp)
for (i in 1:nspp){
  ar_Nt_pop[[i]]=Nt_suit[[i]]$counts*res(daphC_Nt[[i]])[1]^2
  
}

for (i in 1:nspp){
  Nt_suit[[i]]$ar_Nt_pop=Nt_suit[[i]]$counts*res(daphC_Nt[[i]])[1]^2
}


###Population dynamic habitat suitability (competition) barplot
bar_Nt_suit = vector("list",nspp)
for (i in 1:nspp){
  bar_Nt_suit[[i]]=barplot(Nt_suit[[i]]$ar_Nt_pop, Nt_suit[[i]]$mids, ylab="Area (Km^2)", xlab="Suitability", col=NULL)
}

###Population dynamic habitat suitability (competition) boxplot
boxplot_Nt_suit = vector("list",nspp)
par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  boxplot_Nt_suit[[i]]=boxplot(daphC_Nt[[i]], ylab="Population density (number of individuals)", col=NULL)
}

#==============================================================================
###Graphs - plot current habitat suitability boxplots together
#==============================================================================
#Convert no competition and competition rasters to dataframe:
#Competition
df_Nt = vector("list",nspp)
for (i in 1:nspp){
  df_Nt[[i]]=as.data.frame(daphC_Nt[[i]])
}

#No competition
df_Nt_base=vector("list", nspp)
for (i in 1:nspp){
  df_Nt_base[[i]]=as.data.frame(daphC_Nt_base[[i]])
}

col_use = c("light blue","pink")
legend_use = c("No competition", "Competition")
#boxplot of Nt_base (no competition) and Nt (competition)
boxplots = vector("list",nspp)
par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
boxplots[[i]]=boxplot(df_Nt_base[[i]]$layer, df_Nt[[i]]$layer,
                      names=c("No competition", "Competition"),
                      col=c("light blue", "pink"),
                      ylab="Population density")

}

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
 addscalebar(pos="bottomleft")
 addnortharrow(scale=0.5, pos="bottomright")

  
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
co_suit = matrix(0,ny,nx)
#x = l1_tmp[,,i] - co_suit

for ( i in 1:nspp) {
  
  l1_tmp[,,i] = raster::as.matrix(daphC_pred_future[[i]]) #maxent suitability

}

#Set the reproduction rate: 
l1_a = 5

#Remove sites that start out at (approximately) 0 suitability
s_thresh = 0.01
lambda_i [lambda_i < s_thresh] = NA

#Constant lambda across time
#no variation in time
#equivalent to converting suitability to a per-capita reproduction
lambda_i = l1_tmp *l1_a


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
aij1 = 0.1
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
    intra_comp = rowSums(intra_comp, dim=2, na.rm=T) 
    
    Nt[,,s,t+1] = ( lambda_i[,,s] / (1+ aij[s,s] *Nt[,,s,t]+ intra_comp) + 
                      s_i[s] ) * Nt[,,s,t]
    
    
    #The model with no competition
    #just species competing with itself
    Nt_base [,,s,t+1] =  ( lambda_i[,,s]/(1+ aij[s,s] *Nt_base[,,s,t]) + s_i[s] ) * Nt_base[,,s,t]

    
  }
  
}


#Nt_base - no competition
#Nt - competition
max(Nt_base[,,,tend], na.rm=T)
max(Nt[,,,tend], na.rm=T)

#change in population between no competition and competition
delta_pop_fut = (Nt_base[,,,tend] - Nt[,,,tend])/Nt_base[,,,tend]
dim(delta_pop_fut)

#No competition - equivalent to maxent suitability
Nt_base_pop_fut = Nt_base[,,,tend]
dim(Nt_base_pop_fut)

#Interspecific competition
Nt_pop_fut = Nt[,,,tend]
dim(Nt_pop_fut)


###Convert these to rasters with same qualities as original rasters:
daphC_delta_pop_fut = vector("list", nspp)
daphC_Nt_base_fut=vector("list", nspp)
daphC_Nt_fut=vector("list", nspp)


for ( i in 1:nspp) {
  
  daphC_delta_pop_fut[[i]]=daphC_pred_future[[i]]
  daphC_Nt_base_fut[[i]]=daphC_pred_future[[i]]
  daphC_Nt_fut[[i]]=daphC_pred_future[[i]]
 
  daphC_delta_pop_fut[[i]]@data=raster(delta_pop_fut[,,i])@data
  daphC_Nt_base_fut[[i]]@data=raster(Nt_base_pop_fut[,,i])@data
  daphC_Nt_fut[[i]]@data=raster(Nt_pop_fut[,,i])@data
}


#==============================================================================
# Compare and Plot population density-suitability .
#==============================================================================
###Future maxEnt suitability: 
fig.name = paste("daphnia_maxent_suit_fut.pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  
  plot(daphC_pred_future[[i]], main="Predicted Future Suitability")
  plot(wrld_simpl, fill=FALSE, add=TRUE)
  points(daph_sppC[[i]]$lon, daph_sppC[[i]]$lat, pch="+", cex=0.2)
  addscalebar(pos="bottomright")
  addnortharrow(scale=0.5, pos="bottomleft")
}


dev.off()

###Future population density with no interspecific competition: equivalent to maxent 
fig.name = paste("daphnia_Nt_base_fut.pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)
par (mfrow = c(3,1))
for (i in 1:nspp){
  
  plot(daphC_Nt_base_fut[[i]], main=" ")
  plot(wrld_simpl, fill=FALSE, add=TRUE)
  points(daph_sppC[[i]]$lon, daph_sppC[[i]]$lat, pch=20, cex=0.75)
  addscalebar(pos="bottomright")
  addnortharrow(scale=0.5, pos="bottomleft")
  
}


dev.off()

###Future population density with interspecific competition: 
fig.name = paste("daphnia_Nt_fut.pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)
par (mfrow = c(3,1))
for (i in 1:nspp){
  
  plot(daphC_Nt_fut[[i]], main=" ")
  plot(wrld_simpl, fill=FALSE, add=TRUE)
  points(daph_sppC[[i]]$lon, daph_sppC[[i]]$lat, pch=20, cex=0.75)
  addscalebar(pos="bottomright")
  addnortharrow(scale=0.5, pos="bottomleft")
  
}
dev.off()


###Final future change in suitability
fig.name = paste("daphnia_pop_driven_change_fut.pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)
options(repr.plot.width = 2, repr.plot.height = 8)

par (mfrow = c(3,1))
for (i in 1:nspp){
  
  plot(daphC_delta_pop_fut[[i]], main=" ")
  plot(wrld_simpl, fill=FALSE, add=TRUE)
  points(daph_sppC[[i]]$lon, daph_sppC[[i]]$lat, pch=20, cex=0.75)
  addscalebar(pos="bottomright")
  addnortharrow(scale=0.5, pos="bottomleft")
  
}

dev.off()

#==============================================================================
###Graphs - everything plotted by itself
#==============================================================================
###Future Maxent habitat suitability histogram

maxent_suit_fut = vector("list",nspp)
par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  maxent_suit_fut[[i]]=hist(daphC_pred_future[[i]], breaks=seq(from = 0, to = 1, by = 0.1), main="", ylab="Frequency", xlab="Habitat suitability", col=NULL)
}

#Get area of one cell
ar_maxent_fut = vector("list",nspp)
for (i in 1:nspp){
  ar_maxent_fut[[i]]=maxent_suit_fut[[i]]$counts*res(daphC_pred_future[[i]])[1]^2
}

for (i in 1:nspp){
  maxent_suit_fut[[i]]$ar_maxent_fut=maxent_suit_fut[[i]]$counts*res(daphC_pred_future[[i]])[1]^2
}


#Future Maxent habitat suitability barplot
bar_maxent_fut = vector("list",nspp)
for (i in 1:nspp){
  bar_maxent_fut[[i]]=barplot(maxent_suit_fut[[i]]$ar_maxent_fut, maxent_suit_fut[[i]]$mids, ylab="Area (Km^2)", xlab="Suitability", col=NULL)
}

#Future Maxent habitat suitability boxplot
boxplot_maxent_suit_fut = vector("list",nspp)
for (i in 1:nspp){
  boxplot_maxent_suit_fut[[i]]=boxplot(maxent_suit_fut[[i]]$ar_maxent_fut, ylab="Area", col=NULL)
}

###Future Change in population dynamic habitat suitability histogram
delta_pop_suit_fut = vector("list",nspp)
par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  delta_pop_suit_fut[[i]]=hist(daphC_delta_pop_fut[[i]], breaks=seq(from = 0, to = 1, by = 0.1), main="", ylab="", xlab="Suitability", col=NULL)
}

#Get area of one cell
ar_delta_pop_fut = vector("list",nspp)
for (i in 1:nspp){
  #ar_maxent[[i]]= maxent_suit[[i]]$counts*res(daphC_pred[[i]])[1]^2
  ar_delta_pop_fut[[i]]= delta_pop_suit_fut[[i]]$counts*res(daphC_delta_pop_fut[[i]])[1]^2
  
}

for (i in 1:nspp){
  delta_pop_suit_fut[[i]]$ar_delta_pop_fut=delta_pop_suit_fut[[i]]$counts*res(daphC_delta_pop_fut[[i]])[1]^2
}

#Future Change in population dynamic habitat suitability barplot
bar_delta_pop_suit_fut = vector("list",nspp)
for (i in 1:nspp){
  bar_delta_pop_suit_fut[[i]]=barplot(delta_pop_suit_fut[[i]]$ar_delta_pop_fut, delta_pop_suit_fut[[i]]$mids, ylab="Area (Km^2)", xlab="Suitability", col=NULL)
}

#Future Change in population dynamic habitat suitability boxplot
boxplot_delta_pop_suit_fut = vector("list",nspp)
for (i in 1:nspp){
  boxplot_delta_pop_suit_fut[[i]]=boxplot(daphC_delta_pop_fut[[i]], ylab="Habitat suitability", col=NULL)
}


###Future Population dynamic habitat suitability (no competition) histogram
Nt_base_suit_fut = vector("list",nspp)
par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  Nt_base_suit_fut[[i]]=hist(daphC_Nt_base_fut[[i]], breaks=seq(from = 0, to = 10, by = 1), main="", ylab="", xlab="Suitability (population units)", col=NULL)
}

#Get area of one cell
ar_Nt_base_pop_fut = vector("list",nspp)
for (i in 1:nspp){
  #ar_maxent[[i]]= maxent_suit[[i]]$counts*res(daphC_pred[[i]])[1]^2
  ar_Nt_base_pop_fut[[i]]= Nt_base_suit_fut[[i]]$counts*res(daphC_Nt_base_fut[[i]])[1]^2
  
}

for (i in 1:nspp){
  Nt_base_suit_fut[[i]]$ar_Nt_base_pop_fut=Nt_base_suit_fut[[i]]$counts*res(daphC_Nt_base_fut[[i]])[1]^2
}

###Future Population dynamic habitat suitability (no competition) barplot
bar_Nt_base_suit_fut = vector("list",nspp)
for (i in 1:nspp){
  bar_Nt_base_suit_fut[[i]]=barplot(Nt_base_suit_fut[[i]]$ar_Nt_base_pop_fut, Nt_base_suit_fut[[i]]$mids, ylab="Area (Km^2)", xlab="Suitability", col=NULL)
}

###Future Population dynamic habitat suitability (no competition) boxplot
boxplot_Nt_base_suit_fut = vector("list",nspp)
for (i in 1:nspp){
  boxplot_Nt_base_suit_fut[[i]]=boxplot(Nt_base_suit_fut[[i]]$ar_Nt_base_pop_fut, ylab="Area", col=NULL)
}


###Future Population dynamic habitat suitability (competition) histogram
Nt_suit_fut = vector("list",nspp)
par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  Nt_suit_fut[[i]]=hist(daphC_Nt_fut[[i]], breaks=seq(from = 0, to = 10, by = 1), main="", ylab="", xlab="Suitability (population units)", col=NULL)
}

#Get area of one cell
ar_Nt_pop_fut = vector("list",nspp)
for (i in 1:nspp){
  ar_Nt_pop_fut[[i]]=Nt_suit_fut[[i]]$counts*res(daphC_Nt_fut[[i]])[1]^2
  
}


for (i in 1:nspp){
  Nt_suit_fut[[i]]$ar_Nt_pop_fut=Nt_suit_fut[[i]]$counts*res(daphC_Nt_fut[[i]])[1]^2
}


###Population dynamic habitat suitability (competition) barplot
bar_Nt_suit_fut = vector("list",nspp)
for (i in 1:nspp){
  bar_Nt_suit_fut[[i]]=barplot(Nt_suit_fut[[i]]$ar_Nt_pop_fut, Nt_suit_fut[[i]]$mids, ylab="Area (Km^2)", xlab="Suitability", col=NULL)
}

###Population dynamic habitat suitability (competition) boxplot
boxplot_Nt_suit_fut = vector("list",nspp)
for (i in 1:nspp){
  boxplot_Nt_suit_fut[[i]]=boxplot(Nt_suit_fut[[i]]$ar_Nt_pop_fut, ylab="Area", col=NULL)
}


#==============================================================================
###Graphs - future habitat suitability boxplots together
#==============================================================================
#Convert no competition and competition rasters to dataframe:
#Competition
df_Nt_fut = vector("list",nspp)
for (i in 1:nspp){
  df_Nt_fut[[i]]=as.data.frame(daphC_Nt_fut[[i]])
}

#No competition
df_Nt_base_fut=vector("list", nspp)
for (i in 1:nspp){
  df_Nt_base_fut[[i]]=as.data.frame(daphC_Nt_base_fut[[i]])
}


col_use = c("light blue","pink")
legend_use = c("No competition", "Competition")
#boxplot of Nt_base (no competition) and Nt (competition)
boxplots_fut = vector("list",nspp)
#par (mfrow = c(3,1))
par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  boxplots_fut[[i]]=boxplot(df_Nt_base_fut[[i]]$layer, df_Nt_fut[[i]]$layer,
                                                 names=c("No competition", "Competition"),
                                                 col=c("light blue", "pink"), 
                                                 border="black",
                                                 #outline=FALSE,
                                                 ylab="Population density (number of individuals)"
                                               )
  
    legend("topright", 
           legend = legend_use, text.font=c(3),
           col = col_use, 
           pch = c(19,19), 
           bty = "n", 
           cex = 1.3, 
           text.col = "black", 
           inset = c(0.0001, 0.01))
}


#==============================================================================
###Graphs - plot change in future and current boxplots together
#==============================================================================
#Convert future and current change in habitat suitability rasters to dataframe:
#Current
df_delta_pop = vector("list",nspp)
for (i in 1:nspp){
  df_delta_pop[[i]]=as.data.frame(daphC_delta_pop[[i]])
}

#Future
df_delta_pop_fut=vector("list", nspp)
for (i in 1:nspp){
  df_delta_pop_fut[[i]]=as.data.frame(daphC_delta_pop_fut[[i]])
}


col_use = c("light blue","pink")
legend_use = c("No competition", "Competition")

#boxplot of Nt_base (no competition) and Nt (competition)
boxplots_delta_pop = vector("list",nspp)
#par (mfrow = c(3,1))
par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  boxplots_delta_pop[[i]]=boxplot(df_delta_pop[[i]]$layer, df_delta_pop_fut[[i]]$layer, 
                            names=c("Current", "Future"),
                            col=NULL, 
                            border="black",
                            #outline=FALSE,
                            ylab="Habitat suitability")
  
}

#==============================================================================
###Graphs - plot future and current boxplots together 
#==============================================================================
#boxplot comparing current and future habitat suitbaility predictions of Nt_base (no competition) and Nt (competition)
boxplots_all = vector("list",nspp)
par (mfrow = c(1,ceiling(nspp)))
for (i in 1:nspp){
  boxplots_all[[i]]=boxplot(df_Nt_base[[i]]$layer, df_Nt[[i]]$layer, df_Nt_base_fut[[i]]$layer, df_Nt_fut[[i]]$layer, 
                            names=c("                Current", "", "                 Future", ""),
                            col=c("light blue", "pink"), 
                            border="black",
                            #outline=FALSE,
                            ylab="Suitability (Population units)")
  
  for (i in 1:nspp){
    #par(cex.lab=1.2) # is for y-axis
    
    
  par(cex.axis=1.2) # is for x-axis
  }
}


###################################################################################
####END OF CODE