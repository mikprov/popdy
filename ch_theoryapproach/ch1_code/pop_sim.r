# Simulate Populations using spawning contributions at age
# Goal: set up function that simulations populations using 
# contributions of spawning at age for any longevity. 
# Plan:
# 1.

# ===================================================================
# 1) load libraries
#library(ggplot2)
#library(reshape2)
#library(RColorBrewer)
#library(zoo)
#library(tidyr)
#library(dplyr)
#library(ggrepel)
#library(gridExtra)

initial_eggs = beta # delete from function
timesteps = 500
sig_r
beta = 10000
alpha = 1.2

A = matrix(0,nrow=maxage,ncol=maxage)
A[1,] <- c(0,0,0,0,0,0.1,0.1,0.6,0.1,0.1)
for(i in 1:(length(A[1,])-1)){
  A[i+1,i] <- 1
}
# ===================================================================
# 2) define sim model
sim_model <- function(A,timesteps,alpha,beta,sig_r,initial_eggs) {
  
  # create Leslie matrix
  age_at_mat <- min(which(!A[1,] ==0)) # set inequality to min threshold for maturity

  maxage = 10 #length(A[,1])
  
  ages = length(seq(1,maxage,by=1))
  N0 = rep(initial_eggs,length=maxage) #vector, length=num of ages, first age is initial eggs
  set.seed(2) #Set the seed so that every simulation uses same random sequence
  
  # initialize timeseries matrix before simulation
  Nt = matrix(0,ages,timesteps) #rows=ages, cols=time
  Nt[,1] = N0 #Put in initial values
  Nt[,2] = A %*% Nt[,1]  # multiply initial age vector with Leslie to get 2nd age vector
  eggs = c(Nt[1,1],Nt[1,2], rep(NA,timesteps-2)) #save egg production here, this will be the input in BH
  recruits = rep(NA, length=timesteps-2) #will save recruits here (output from BH)
  
  for(t in 1:(timesteps-2)){ #step through time
    #Save egg production, this is new number of age 1 individuals
    
    hold = A %*% Nt[,t+1] #Leslie matrix * age vector
    eggs[t+1] = hold[1,] #store egg production, this will be used in BH
    #save recruits
    recruits[t] = eggs[t]/( (1/alpha) + (eggs[t]/beta) ) #((alpha*eggs[t+1])/(1+beta*eggs[t+1])) 
    #replace age 1 with recruits from BH, add noise    
    hold[1,] = (recruits[t])*exp(sig_r*rnorm(1,mean=0,sd=1)) 
    #put new age vector, with recruits, into time series matrix
    Nt[,t+2] <- hold
    
  }
  plot(eggs,type="l")
  plot(colSums(Nt),type="l")
  plot(recruits,type="l")
  #Nsize = colSums(Nt) #total population size, sums num at age for each timestep
  Nsize = colSums(Nt[age_at_mat:maxage,]) #sum rows that correspond to spawning adults
  N_t = Nt[1,][2:(timesteps-1)] #Nt is number of spawning adults, aka eggs, after environmental noise
  eggs = eggs[2:(timesteps-1)] #egg production
  recruits = recruits[2:(timesteps-1)] #recruits produced via BH, before influence of environmental noise
  
  return(list(N_t=N_t, eggs=eggs, recruits=recruits, Nsize=Nsize))
}

# ===================================================================


