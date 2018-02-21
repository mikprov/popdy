# North Atlantic Cod Populations
# by: Mikaela Provost

# Plan

#library(foreach)
rm(list=ls())
# ===================================================================
# 1) specify parms 
maxage = 40
ages = 1:maxage
# If I add populations, I need to specify their names in codPop
# And I need to make sure a script is saved in [path] and is labeled w/the name in codPop
codPop <- c("Celtic_Sea","3NO","North_Sea","West_Scotland",
            "West_Baltic","Kattegat","Coastal_Norway",
            "Faroes","Iceland","NE_Arctic","Georges_Bank",
            "Gulf_of_Maine","N_Gulf_St_Lawrence")
    #Irish Sea removed because Beta coefficients not significnat

# ===================================================================
# 3) set up emtpy objects for function


# ===================================================================
# 2) define function

simulate <- function(M,H,) {
  # fill in abundance at age for time=1
  n[,1] = R_0*exp(-(M+H)*(0:(maxage-1)))
  
  # for each time step store results from 
  # computing Leslie matrix
  for(j in 2:time){
 
     leslie_eggs <- assemble_Leslie(parms...)
 
     # extract Leslie matrix for time j
     A <- leslie_eggs$A
     # store egg production for time j
     eggs[j] <- leslie_eggs$E
  
     # perform population projection of time step
     n[,j] <- A %*% n[,j-1]
  
     # make sure abundances at age <1 become zero
     for(i in 1:maxage) {if(n[i,j]<1) n[i,j]=0}
     # if whole population goes below 1, stop sim
     if(sum(n[,j])<1) print("population crashed!!")
  }
  
  # Return list including numbers at age, egg production
  # for the whole timeseries, with the first 100 time
  # steps trimmed off (to make sure we are around the 
  # equilibrium)
  # n = a matrix with maxage rows, and time-100 cols
  # eggs = a vector, egg production each yr
  return(list(n=n,eggs=eggs))
     
}




# ===================================================================
# 3) apply functions to different cod populations
H_levels = seq(0,2,0.01) #specify a range of harvest values
FLSBlist = as.list(rep(NA,length(codPop)))
names(FLSBlist) <- codPop
LSBlist = as.list(rep(NA,length(codPop)))
names(LSBlist) = codPop



# ===================================================================
# 4) plot fraction of lifetime spawning biomass for varying F rates
# make sure to take the natural log of the fraction of unfished lifetime spawning biomass, ln(FLSB)






# ===================================================================
# 5) plot spawning biomass for each population, and fished at 3 levels (H=0,0.5,1)

# 5a) create function that calculates spawning biomass at age vector
SBagevector <- function(H) {
  w = wt_at_age(L_inf=L_inf, K=K, ages=ages) 
    #function to get weight @ age (using vonB eq then convert L->W)
  A = maturity_ogive(B0=B0,B1=B1,ages=ages) 
    #input ogive coefficients to calculate vector of maturity probabilities
  s = survival(H=H,M=M,theta0=theta0,theta1=theta1,ages=ages) 
    #function to calc survival @ age, has fishing mortality (H), M, and selectivity (S)
  SBatage = w * A * s
  return(SBatage)
}


# 5b) loop through populations to create spawning biomass at age vectors using the 
# function above and create a spawning biomass distribution plot for each population



# ===================================================================
# 6) spectral frequency plots (one for each population)
#