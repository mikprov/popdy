# North Atlantic Cod Populations - simulate timeseries
# by: Mikaela Provost

 

# ===================================================================
# 1) load libraries
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(zoo)
library(tidyr)
library(dplyr)
library(ggrepel)
library(gridExtra)

# ===================================================================
# 2) define sim model
sim_model <- function(A,timesteps,alpha,beta,sig_r,initial_eggs) {
  age_at_mat <- min(which(!A[1,] ==0)) # set inequality to min threshold for maturity
  maxage = length(A[,1])
  ages = length(seq(1,maxage,by=1))
  N0 = c(initial_eggs, rep(0,ages-1)) #vector, length=num of ages, first age is initial eggs
  set.seed(2) #Set the seed so that every simulation uses same random sequence
  
  Nt = matrix(0,ages,timesteps) #empty matrix to store output: rows=ages, cols=time
    #Initialize vector of population sizes with extra 
    #rows for egg production (top value in age vector) & 
    #recruitment before variability (output from BH)
  Nt[,1] = N0 #Put in initial values
  Nt[,2] = A %*% Nt[,1]  # multiply initial age vector with Leslie to get 2nd age vector
  eggs = c(Nt[1,1], rep(NA,timesteps-1)) #will save egg production here, this will be the input in BH
  recruits = c(initial_eggs, rep(0, timesteps-1)) #will save recruits here (output from BH)
 
    for(t in 1:(timesteps-2)){ #step through time
      #Save egg production, this is new number of age 1 individuals
      eggs[t+1] = Nt[1,t+1] 
      #save recruits, treat egg production as spawners in BH
      recruits[t+1] = eggs[t+1]/( (1/alpha) + (eggs[t+1]/beta) ) #((alpha*eggs[t+1])/(1+beta*eggs[t+1])) 
      #replace age 1 with recruits from BH, add noise    
      Nt[1,t+1] = (recruits[t+1])*exp(sig_r*rnorm(1,mean=0,sd=1)) 
      #perform population projection for one time step   
      Nt[,t+2] = A %*% Nt[,t+1] 
        
    }
  #Nsize = colSums(Nt) #total population size, sums num at age for each timestep
  Nsize = colSums(Nt[age_at_mat:maxage,]) #sum rows that correspond to spawning adults
  N_t = Nt[1,][2:(timesteps-1)] #Nt is number of adults, aka eggs
  eggs = eggs[2:(timesteps-1)] #egg production
  recruits = recruits[2:(timesteps-1)]
  
  return(list(N_t=N_t, eggs=eggs, recruits=recruits, Nsize=Nsize))
}

# ===================================================================


