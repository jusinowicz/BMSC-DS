#=============================================================================
# Implementations of the Lotka-Volterra and the  Leslie-Gower models. 
# This code is to compare discrete-time and continuous time implementations
# of models for two competing species.
# 
# A continuous time population model has the form:
# 	dN/dt = g(N..)*N 
# where N is the population, g(N...) is a growth rate function, and t is time
# 
# A discrete time model has the form: 
#	N[t+1] = g(N...)* N[t]
#
# The Leslie-Gower model for two competing species is defined by 
#   g_i(N_i...) = lambda_i/ (1+ a_ii *N_i+ a_ij*N_j) + s_i
# 
# 
# The subscript i indicates the species number, 
# lambda is the intrinsic rate of reproduction (in the absence of competition)
# a_ii, a_ij  are intra and interspecific competition
# s_i	is survival
#
# The Leslie-Gower model is meant to be a discrete-time analogue of the 
# Lotka-Volterra model. When models are deterministic (i.e. no stochasticity
# or environmental variation ) they are effictively the same model. 
#=============================================================================
# load libraries
#=============================================================================
library(tidyverse)
library(deSolve)

#=============================================================================
# 1. Continuous time model using the deSolve library
#=============================================================================
#=============================================================================
# In deSolve you define a model as a function that is passed 3 variables: 
# times  the time series with start, finish, and increment sizes 
#	(defined below)
# sp 	 a matrix that stores the output from increment to increment
# parms  all of the variables that define parameters in the model
#=============================================================================
LV_ct = function(times,sp,parms){

    nspp = parms$nspp #number of species. 2 for now
    N = matrix(sp[1:nspp], nspp, 1)
			    
	###Population dynamics
				dN = N

		#I made a for loop to do both species because they have identical
		#dynamics. 

		for( i in 1:nspp){ 

			#Lotka-Volterra dynamics
			#Note: This assumes that both intra and interspecific competition
			#are in the same matrix "aij." Use R's negative subscript indexing
			#to access the off-diagonal aij for interspecific competition. 
			dN[i] = ( (lambda_i[i] - aij[i,i] *N[i] - aij[i,(-i)]*N[-i]) ) * N[i]

		}

	return( list(dN) )

}  

#=============================================================================
# Set values of the population parameters 
#=============================================================================
nspp = 2
#Intrinsic reproduction:
lambda_i = c(1.2,1.1)

#Competition: 
#This creates a square matrix for aij with intraspecific competition = 1 on 
#the diagonal. 
aij = matrix ( c(1, 0.8, 
	0.8 , 1 ), nspp, nspp )

#Survival (should be less than 1)
s_i = c(0.9,0.9)

#Make sure this has everything that was definied in the function above. 
parms = list(
			nspp = nspp, lambda_i = lambda_i, aij=aij, s_i=s_i
		 )

#=============================================================================
# Run the model with initial conditions 
#=============================================================================
tend = 100
delta1 = 0.1
times  = seq(from = 0, to = tend, by = delta1)
tl = length(times)
winit = c(matrix(0.1,nspp,1)) #Initial population sizes

Nt_ct= ode(y=winit, times=times, func=LV_ct, parms=parms)

#=============================================================================
# Plot
#=============================================================================
plot( Nt_ct [,2], t="l", ylab = "Time", xlab = "Population density")
lines( Nt_ct [,3],  col = "red")



#=============================================================================
# 2. Discrete-time model
#=============================================================================
#Create the population matrix with the same times as continuous-time model: 
Nt = matrix(0,tend,nspp)
#Initial conditions:
Nt[1,] = c(0.1,0.1)

#Use the same parameters as the continuous-time model, or re-assign values: 

#Intrinsic reproduction:
lambda_i = c(1.2,1.1)

#Competition: 
#This creates a square matrix for aij with intraspecific competition = 1 on 
#the diagonal. 
aij = matrix ( c(1, 0.8, 
	0.8 , 1 ), nspp, nspp )

#Survival (should be less than 1)
s_i = c(0.9,0.9)

#Loop over all of time: 
for (t in 1:(tend-1) ) {

	#Loop over all species: 
	for (s in 1:nspp) {

		Nt[t+1,s] = ( lambda_i[s] / (1+ aij[s,s] *Nt[t,s]+ aij[s,(-s)]*Nt[t,-s]) + 
					s_i[s] ) * Nt[t,s]
	}

}

#=============================================================================
# Plot
#=============================================================================
plot( Nt [,1], t="l", ylab = "Time", xlab = "Population density", ylim=c(0,max(Nt)))
lines( Nt [,2],  col  = "red")


#=============================================================================
# Joint Plot
#=============================================================================
plot( Nt [,1], t="l", ylab = "Time", xlab = "Population density", ylim=c(0,max(Nt)))
lines( Nt [,2],  col = "red")
lines( Nt_ct [,2], t="l", lty=3)
lines( Nt_ct [,3], lty=3, col = "red")

