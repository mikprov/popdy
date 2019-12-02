# Check eigenvalues and dynamics of Jacobian when lambda1 < and > 1
# last updated: 2019 Aug 5

# ===================================================================
# FUNCTIONS

# 1) sim_model, simulates the Jacaobian without density dependence
sim_model <- function(A,timesteps,initial_eggs,maxage) {
  age_at_mat <- min(which(!A[1,] ==0)) # set inequality to min threshold for maturity, for cod it is always age 1
  ages = length(seq(1,maxage,by=1))
  N0 = c(initial_eggs, rep(0,ages-1)) #vector, length=num of ages, first age is initial eggs
  set.seed(2) #Set the seed so that every simulation uses same random sequence
  
  Nt = matrix(0,ages,timesteps) #empty matrix to store output: rows=ages, cols=time
  
  Nt[,1] = N0 #Put in initial values
  Nt[,2] = A %*% Nt[,1]  # multiply initial age vector with Jacobian to get 2nd age vector
  eggs = c(initial_eggs, rep(NA,timesteps-2)) #will save egg production here, this will be the input in BH
  eggs[2] = Nt[1,2]
  
  for(t in 1:(timesteps-2)){ #step through time
    #Save egg production, this is new number of age 1 individuals
    eggs[t+1] = Nt[1,t+1] 
    
    #perform population projection for one time step   
    Nt[,t+2] = A %*% Nt[,t+1] 
  }
  Nsize = colSums(Nt[age_at_mat:maxage,]) #sum rows that correspond to spawning adults
  eggs = eggs[1:(timesteps-1)] #egg production
  return(list(eggs=eggs, Nsize=Nsize))
}

# 2) sim_model_dd, simulates the Jacobian with density dependence and noise
sim_model_dd <- function(A,timesteps,alpha,beta,sig_r,initial_eggs) {
  age_at_mat <- min(which(!A[1,] ==0)) # set inequality to min threshold for maturity, for cod it is always age 1
  
  ages = length(seq(1,maxage,by=1))
  N0 = c(initial_eggs, rep(0,ages-1)) #vector, length=num of ages, first age is initial eggs
  set.seed(2) #Set the seed so that every simulation uses same random sequence
  
  Nt = matrix(0,ages,timesteps) #empty matrix to store output: rows=ages, cols=time
  
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
    Nt[1,t+1] = (recruits[t+1])*exp(sig_r*rnorm(1,mean=0,sd=1)) 
    #perform population projection for one time step   
    Nt[,t+2] = A %*% Nt[,t+1] 
  }
  Nsize = colSums(Nt[age_at_mat:maxage,]) #sum rows that correspond to spawning adults
  eggs = eggs[1:(timesteps-1)] #egg production
  recruits = recruits[1:(timesteps-1)] #recruits produced via BH, before influence of environmental noise
  return(list(eggs=eggs, recruits=recruits, Nsize=Nsize))
}

# 3) laods the function: extract_first_eigen_value()
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/2_cod_functions.r")


# ===================================================================
#source(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_pops/','Northsea', '.r', sep=''))

# North Sea example
maxage = 10
K = 0.217
L_inf = 126
F.halfmax=0
B0 = -6.42
B1 = 1.72
tknot=0
MG = 0.427
# Generate Leslie matrix (A) and life table (NEAR)
Leslieout = assemble_Leslie(maxage=maxage, K=K, L_inf=L_inf, 
                            F.halfmax=0, B0=B0, B1=B1, tknot=0)
NEAR = Leslieout$NEAR #Life table
NEAR
LEP = sum(NEAR$egg_production) #sum egg production (not adjusted)
LEP #Lifetime egg production
1/LEP #slope of 1/LEP

A = Leslieout$A #Leslie matrix 
A

# ===================================================================
# Simulate Leslie and Jacobian *without* density dependence

# First, we want to standardize LEP across populations; to do this we adjust fecundities 
# to make LEP=1.1. Fecundity-at-age should be divided by LEP*0.9: f-at-age/(LEP*0.9)
conLEP = 1.1 #constant LEP
adjFec = round(1/conLEP,digits=1) #1/LEP is 1.1, which means alpha must be greater than this
A[1,] <- A[1,]/(LEP*adjFec) #adjusted Leslie matrix
extract_first_eigen_value(A) #lambda1 of Leslie matrix w/adjusted fecundity
# Simulate this Leslie matrix:
Leslie_output_withoutDD_lambda_small <- sim_model(A=A,timesteps=1000,initial_eggs=100,maxage=10)
plot(Leslie_output_withoutDD_lambda_small$Nsize,xlab="year",ylab="Pop size",main="Leslie simulation without\ndensity dependence (lambda1=1.02)")
# We see that population size increases exponentially because lambda1 > 1


# Second, convert Leslie to Jacobian: multiply top row of Leslie by 1/(alpha*LEP^2)
# The slope at equilibrium is 1/(alpha*LEP^2) (see ch 4 in Loo's book or my appendix for derivation)
# Since we adjusted LEP to 1.1, we multiply by: 1/(alpha*1.1^2) (we call this 'k')
alpha = 5.51 #with this value of alpha, slope at equilibrium is 0.15
slope_at_equilibrium = 1/(alpha*conLEP^2) #we call the slope 'k'
slope_at_equilibrium #approximately 0.15

A[1,] <- A[1,]*(1/(alpha*conLEP^2)) #convert to Jacobian
# We can check lambda1 of the adjusted Jacobian. 
lambda1 <- extract_first_eigen_value(A)
lambda1 #we see that lambda1=0.77


# When lambda1 is smaller (lambda1 = 0.77):
Jacobian_output_withoutDD_lambda_small <- sim_model(A=A,timesteps=1000, initial_eggs=100,maxage=10)
plot(Jacobian_output_withoutDD_lambda_small$Nsize,xlab="year",ylab="Pop size",main="Jacobian simulation without\ndensity dependence (lambda1=0.76)")
# We see in the plot above that the population is perturbed initially, but quickly goes back to zero. The population goes to zero because lambda1 < 1. 


# Simulate Jacobian when lambda1 is closer to 1:
# Let's choose a different alpha so that lambda1 will be closer to 1.
alpha0.97 = 0.97
A = Leslieout$A #reset Leslie matrix to original form without adjustments
A[1,] <- A[1,]/(LEP*adjFec) #adjust LEP in Leslie matrix to make LEP=1.1
slope_at_equilibrium = 1/(alpha0.97*conLEP^2) #we call the slope at equilibrium 'k'
slope_at_equilibrium
A[1,] <- A[1,]*slope_at_equilibrium #convert Leslie to the Jacobian by multiplying top row 1/(alpha*LEP^2)
# Let's check the leading eigenvalue of this Jacobian:
extract_first_eigen_value(A) # Here we see that lambda1 is closer to 1 (0.99)
# Simulate Jacobian:
Jacobian_output_withoutDD_lambda_big <- sim_model(A=A,timesteps=1000, initial_eggs=100,maxage=10)
plot(Jacobian_output_withoutDD_lambda_big$Nsize,xlab="year",ylab="Pop size",main="Jacobian simulation without\ndensity dependence (lambda1=0.99)")
# Here, we see that the population goes to zero, but it takes much longer to reach zero because lambda1 is bigger. 



# ===================================================================
# Simulate Leslie and Jacobian *with* density dependence and noise

# First, set up the fecundity-adjusted Leslie matrix:
A = Leslieout$A #reset Leslie matrix to original form without adjustments
A[1,] <- A[1,]/(LEP*adjFec) #adjust LEP in Leslie matrix to make LEP=1.1
extract_first_eigen_value(A) #lambda1 > 1
# Simulate adjusted Leslie with density dependence:
alpha5.51 <- 5.51 # define alpha for BH in DD
Leslie_output_withDD_lambda_small <- sim_model_dd(A=A,timesteps=1000,alpha=5.51,beta=1000,sig_r=0.3,initial_eggs=100)
plot(Leslie_output_withDD_lambda_small$recruits,xlab="year",ylab="recruits (before noise)",main="Leslie simulation with\ndensity dependence (lambda1=1.05)",type="l")

# Second, set up the Jacobian:
slope_at_equilibrium = 1/(alpha5.51*conLEP^2) #we call the slope at equilibrium 'k'
slope_at_equilibrium # k=0.15
A[1,] <- A[1,]*slope_at_equilibrium #convert Leslie to the Jacobian by multiplying top row 1/(alpha*LEP^2)
# Let's check the leading eigenvalue of this Jacobian:
extract_first_eigen_value(A) # Here we see that lambda1=0.77
# Simulate Jacobian w/dd:
Jacobian_output_withDD_lambda1_small <- sim_model_dd(A=A,timesteps=1000,alpha=5.51,beta=1000,sig_r=0.3,initial_eggs=100)
plot(Jacobian_output_withDD_lambda1_small$recruits,xlab="year",ylab="recruits (before noise)",main="Jacobian simulation with\ndensity dependence (lambda1=0.77)",type="l")


# Next, simulate Leslie and Jacobian for a different alpha value
alpha0.97 <- 0.97
A = Leslieout$A #reset Leslie matrix to original form without adjustments
A[1,] <- A[1,]/(LEP*adjFec) #adjust LEP in Leslie matrix to make LEP=1.1
extract_first_eigen_value(A) #lambda1 > 1
# Simulate adjusted Leslie with density dependence:
Leslie_output_withDD_lambda_big <- sim_model_dd(A=A,timesteps=1000,alpha=0.97,beta=1000,sig_r=0.3,initial_eggs=100)
plot(Leslie_output_withDD_lambda_big$recruits,xlab="year",ylab="recruits (before noise)",main="Leslie simulation with\ndensity dependence (lambda1=1.05)",type="l")

# Next, set up the Jacobian:
slope_at_equilibrium = 1/(alpha0.97*conLEP^2) #we call the slope at equilibrium 'k'
slope_at_equilibrium # k=0.85
A[1,] <- A[1,]*slope_at_equilibrium #convert Leslie to the Jacobian 
# Let's check the leading eigenvalue of this Jacobian:
extract_first_eigen_value(A) # Here we see that lambda1=0.99
# Simulate model w/dd:
Jacobian_output_withDD_lambda1_big <- sim_model_dd(A=A,timesteps=1000,alpha=0.97,beta=1000,sig_r=0.3,initial_eggs=100)
plot(Jacobian_output_withDD_lambda1_big$recruits,xlab="year",ylab="recruits (before noise)",main="Jacobian simulation with\ndensity dependence (lambda1=0.99)",type="l")




