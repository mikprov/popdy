# Simulate populations with the goal of calculating probability of extinction
# Beverton Holt recruitment 

sim_model_pE <- function(A,timesteps,alpha,beta,sig_r,initial_eggs,noise) {
  age_at_mat <- min(which(!A[1,] ==0)) # set inequality to min threshold for maturity
  maxage = length(A[,1])
  ages = length(seq(1,maxage,by=1))
  #N0 = c(initial_eggs, rep(0,ages-1)) #vector, length=num of ages, first age is initial eggs
  N0 = rep(initial_eggs/maxage,length=ages)
  
  Nt = matrix(0,ages,timesteps) 
  
  
  Nt[,1] = N0 #Put in initial values
  Nt[,2] = A %*% Nt[,1]  # multiply initial age vector with Jacobian to get 2nd age vector
  eggs = c(initial_eggs, rep(NA,timesteps-2)) #will save egg production here, this will be the input in BH
  eggs[2] = Nt[1,2]
  recruits1 = eggs[1]/( (1/alpha) + (eggs[1]/beta) )
  recruits = c(0,rep(NA, timesteps-2)) #will save recruits here (output from BH)
  recruits[2] = recruits1 #we are ignoring recorinding recruits at time 0
  
  for(t in 1:(timesteps-2)){ #step through time
    #Save egg production, this is new number of age 1 individuals
    eggs[t+1] = Nt[1,t+1] 
    #save recruits, treat egg production as spawning biomass in BH
    recruits[t+1] = eggs[t]/( (1/alpha) + (eggs[t]/beta) ) #((alpha*eggs[t+1])/(1+beta*eggs[t+1])) 
    #replace age 1 in numbers-at-age matrix with recruits from BH, add noise    
    Nt[1,t+1] = (recruits[t+1])*exp(sig_r*noise[t+1]) 
    #Nt[1,t+1] = (recruits[t+1])*exp(sig_r*rnorm(1,mean=0,sd=1)) 
    #perform population projection for one time step   
    Nt[,t+2] = A %*% Nt[,t+1] 
  }
  Nsize = colSums(Nt[age_at_mat:maxage,]) #sum rows that correspond to spawning adults
  eggs = eggs[1:(timesteps-1)] #egg production
  recruits = recruits[1:(timesteps-1)] #recruits produced via BH, before influence of environmental noise
  return(list(eggs=eggs, recruits=recruits, Nsize=Nsize))
}

# ===================================================================
