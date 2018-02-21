# North Atlantic Cod Populations
# by: Mikaela Provost

# Plan
# 1) specify parms
# 2) define functions
# 3) apply functions to different cod populations
# 4) plot fraction of lifetime spawning biomass for varying F rates

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
# 2) define functions

# function to get weight at age
wt_at_age <- function(L_inf,K,ages) {
  # first calculate length at age using von Bert eq
  L = L_inf*(1-exp(-K*(ages))) # 'ages' is a vector 1-40; L units = cm
  w = 0.00001*(L^3) # w units = kg
  return(w)
}
# test out the Celtic Sea
# here we see that as age increases (x axis), weight increases and then plateaus at 9.7
# if we plot the lengths from this function (L), then the asympotic length is 99
#plot(wt_at_age(L_inf=99, K= 0.39, ages=ages))


# function to get probability of maturity at age
maturity_ogive <- function(B0,B1,ages) {
  # A is the probability of an individual of being mature at a determined age. 
  # B0 (intcept) and B1 (slope) are parameter estimates
  #A = exp((B0+B1*ages))/(1+exp((B0+B1*ages))) # original attempt
  A = 1/(1+exp(-B0-B1*ages))
  return(A)
}

# function to get survival at age
survival <- function(H,M,theta0,theta1,ages) {
  #si = exp((theta0+theta1*ages))/(1+exp((theta0+theta1*ages))) 
  si = 1/(1+exp(-theta0-theta1*ages))
    #original selectivity curve --> not used 
  #si = exp((B0+B1*ages))/(1+exp((B0+B1*ages))) # original attempt
    #use maturation ogives to represent selectivity curve
  #si = 1/(1+exp(-B0-B1*ages))
  #use maturation ogives to represent selectivity curve
  s = exp(-(M+(H*si))*ages) 
    #probability of surviving to each age based on fishery 
    #selectivity and natural mort
  return(s)
}


LSB <- function(H) {
  w = wt_at_age(L_inf=L_inf, K=K, ages=ages) #function to get weight @ age (using vonB eq then convert L->W)
  A = maturity_ogive(B0=B0,B1=B1,ages=ages) #input ogive coefficients to calculate vector of maturity probabilities
  s = survival(H=H,M=M,theta0=theta0,theta1=theta1,ages=ages) #function to calc survival @ age, has fishing mortality (H), M, and selectivity (S)
  
  lifetimeSB = sum(w * A * s)
  return(lifetimeSB)
  
}

# ===================================================================
# 3) apply functions to different cod populations
H_levels = seq(0,2,0.01) #specify a range of harvest values
FLSBlist = as.list(rep(NA,length(codPop)))
names(FLSBlist) <- codPop
LSBlist = as.list(rep(NA,length(codPop)))
names(LSBlist) = codPop


for (i in seq_along(codPop)) { ## step through different parm sets (associated w/different cod populations)
  source(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_pops/',codPop[i], '.r', sep=''))
       LSBvec = sapply(H_levels,LSB) ## for a range of fishing rates, calculate LSB
       LSBlist[[i]] = LSBvec ## store LSB for different F rates in a list (might want later)
       FLSBvec = LSBvec/LSB(H=0) ## divide each LSB value by LSB when F=0
       FLSBlist[[i]] = FLSBvec ## store each vector (the fraction of LSB for varying F) in a list
       
}

FLSBdf <- data.frame(FLSBlist) #convert list to dataframe, easier for plotting

# ===================================================================
# 4) plot fraction of lifetime spawning biomass for varying F rates
# make sure to take the natural log of the fraction of unfished lifetime spawning biomass, ln(FLSB)
matplot(x=H_levels,
        y=log(FLSBdf),
        type="l",
        lty=1:length(codPop),
        lwd=1.2,
        col=1:length(codPop),
        ylim=c(-2.5,0),
        ylab="ln(FLSB)",
        xlab=expression(F(y^{-1}))
)
abline(h=log(0.35))
legend("topright",legend=codPop,col=1:length(codPop),lty=1:length(codPop),lwd=1.2)


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
H_3levels = c(0,0.5,1)

# set plots to print in a multi-panel plot
par(mfrow=c(4,4))

for (i in seq_along(codPop)) {
  source(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_pops/',codPop[i], '.r', sep=''))
    SBagematrix = sapply(H_3levels, SBagevector)
    sdmatrix = matrix(apply(SBagematrix,2,mean),nrow=maxage,
                      ncol=length(H_3levels),byrow=TRUE)
    SBagematrix_norm = SBagematrix / sdmatrix
    matplot(x=1:maxage, y=SBagematrix, type="l", 
            main=codPop[i], ylim=c(0,1.2),xlab="age")
    
}

par(mfrow=c(1,1))


# 5c) try calculating CV of spawning age vector (CV = sd / mean)
sd(SBagevector(H=0.2)) / mean(SBagevector(H=0.2))
matplot(x=1:maxage,y=sapply(H_levels,SBagevector),type="l")


# ===================================================================
# 6) spectral frequency plots (one for each population)
#