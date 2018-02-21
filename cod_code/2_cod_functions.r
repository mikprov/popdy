# Cod functions
# by: Mikaela Provost

library(faraway)
# --------------------
# 1. assemble_Leslie() 
# --------------------
# Definition of inputs:
# a) data = dataset for one population, this is subsetted from the cod data from Hui-Yu
#
# b) codPopname = name of the cod population, important because Leslie function will
#    loop over multiple cod populations
#
# c) littlek = the value multiple across the top row of the Leslie matrix, 
#    is a measure of fishing pressure (Lauren's analysis)

assemble_Leslie <- function(data,codPopname,littlek) {
  
  data = subset(data,Yearclass>1959)  #yearclasses vary among stocks
  data = subset(data,Yearclass<1990)
  
  # load parms for cod pop
  source(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_pops/',codPopname, '.r', sep=''))
  # this should load parms: L_inf, K, TEMP
  
  # calculate maturity, weight parms
  (mod.mat = glm(MATPROP~AGE,family=binomial,data=data)) #gives maturity beta coefficients
  
  L=L_inf*(1-exp(-K*Age))
  #MG=exp(0.55-1.61*log(L)+1.44*log(L_inf)+log(K)) #Gislason model II
  MG3 = exp(15.11-1.59*log(L)+0.82*log(L_inf)-3891/(273.15+6.75))  #Gisllason model III
  MP = 10^(-0.0066-0.279*log10(L_inf)+0.6543*log10(K)+0.4634*log10(10.56))  #Pauly model
  #Vul1 = data$CANUM/data$STNUM
  growth = 0.00001*(L_inf*(1-exp(-K*(Age-to))))^3 # weight at age
  mat1 = ilogit(mod.mat$coef[1]+mod.mat$coef[2]*Age)
  
  # -- assemble NEAR df: use NEAR to calculate Leslie matrix
  NEAR = data.frame(cbind(Age,mat1,growth))  #life table
  
  NEAR$Vul1 = mat1 #use maturity ogive for selectivity
  NEAR$M_G= 0.19+0.058*TEMP   #regression fit of Fig. 5c, loaded in parms set
  
  A = matrix(0,length(Age),length(Age))
  
  # -- for each age, get F and survival
  for(j in 1:length(Age)){ # step through ages
    NEAR$F[j] = NEAR$Vul1[j]*F.halfmax[1]	# Vul1 should be selectivity, but it's mat
    
    NEAR$SURV = exp(-(NEAR$F+NEAR$M_G)) #SURV is the fraction surviving at each age
    
    #NEAR$Survship = 0 # set up column for survivorship (amount or fraction present at age)
    #NEAR$Survship[1] = 1
    
    #for(k in 1:(nrow(NEAR)-1)){ # step through ages to calc survivorship
    #  NEAR$Survship[k+1] = NEAR$Survship[k]*NEAR$SURV[k] #amount present at age
    #}
  }
  
  
  # --- if I want to multiply fecundities by some k term, use this line --- #
  #A[1,] = NEAR[,2]*NEAR[,3] # insert fecundity (maturity * weight) in top row
  A[1,] = NEAR$mat1*NEAR$growth*littlek # insert fecundity (maturity * weight) in top row
  # ----------------------------------------------------------------------- #
  
  for(u in 2:length(Age)-1){ # insert survival into A on subdiagonal
    A[u+1,u]=NEAR$SURV[u]
  }
  return(A=A)
  return(NEAR=NEAR) # returns Leslie matrix
} # closes assemble_Leslie matrix function




# --------------------
# 2. extract_first_eigen_value()
# --------------------
# Definitions of inputs:
# Lesliematrix = leslie matrix, output from assemble_leslie()
extract_first_eigen_value <- function(Lesliematrix){
  # get leading eigenvalue - check to make sure it's positive
  ev = eigen(Lesliematrix)
  # squaring the imaginary part is not necessary b/c the first 
  # eigen value has no imaginary part. but including it for 
  # consistency 
  a.sq = (Re(ev$values[1]))^2 # square the real part of 1st eigenvalue
  b.sq = (Im(ev$values[1]))^2 # square the imaginary part of 1st eigenvalue
  firstval = sqrt(a.sq + b.sq) # magnitude of 2nd eigenvalue
  return(firstval)
  rm(a.sq) # remove from workpace
  rm(b.sq) # remove from workpace
}


extract_second_eigen_value <- function(Lesliematrix){
  # get magnitude of second eigenvalue
  # think of real value on x-axis, imaginary on y-axis
  # magnitude is the vector between them
  ev = eigen(Lesliematrix)
  a.sq = (Re(ev$values[2]))^2 # square the real part of 2nd eigenvalue
  b.sq = (Im(ev$values[2]))^2 # square the imaginary part of 2nd eigenvalue
  secondval = sqrt(a.sq + b.sq) # magnitude of 2nd eigenvalue
  return(secondval)
  rm(a.sq) # remove from workpace
  rm(b.sq) # remove from workpace
}


# --------------------
# 3. calc_LSB_at_age_by_F()
# --------------------
# Definition of inputs:
# data = dataset for one population, this is subsetted from the cod data from Hui-Yu
# codPopname = name of the cod population, important because Leslie function will
#   loop over multiple cod populations
# littlek = the value multiple across the top row of the Leslie matrix, 
#   is a measure of fishing pressure (Lauren's analysis)
calculate_LSB_at_age_by_F <- function(data,codPopname,littlek){
  
  data = subset(data,Yearclass>1959)  #yearclasses vary among stocks
  data = subset(data,Yearclass<1990)
  
  # -- load parms for cod pop:
  source(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_pops/',codPopname, '.r', sep=''))
  # this should load parms: L_inf, K, TEMP
  
  # -- calculate maturity, weight parms
  (mod.mat = glm(MATPROP~AGE,family=binomial,data=data)) #gives maturity beta coefficients
  L=L_inf*(1-exp(-K*Age))
  #MG=exp(0.55-1.61*log(L)+1.44*log(L_inf)+log(K)) #Gislason model II
  MG3 = exp(15.11-1.59*log(L)+0.82*log(L_inf)-3891/(273.15+6.75))  #Gisllason model III
  MP = 10^(-0.0066-0.279*log10(L_inf)+0.6543*log10(K)+0.4634*log10(10.56))  #Pauly model
  #Vul1 = data$CANUM/data$STNUM
  growth = 0.00001*(L_inf*(1-exp(-K*(Age-to))))^3
  mat1 = ilogit(mod.mat$coef[1]+mod.mat$coef[2]*Age)
  
  # -- assemble NEAR df: use NEAR to calculate LSB
  NEAR = data.frame(cbind(Age,mat1,growth))  #cols: age, maturity, growth (size at age)
  NEAR$Vul1 = mat1 #use maturity ogive for selectivity
  NEAR$M_G= 0.19+0.058*TEMP   #mortality at age, regression fit of Fig. 5c, 'TEMP' is loaded in parms set
  
  # -- 
  LEPdf = matrix(0,length(Age),length(F.halfmax))
  for(g in 1:length(F.halfmax)){ # step through each fishing level
    for(j in 1:length(Age)){ # step through ages
      NEAR$F[j] = NEAR$Vul1[j]*F.halfmax[g]	# Vul1 should be selectivity, but it's mat
      # calculate F = selectivity * F rate
      NEAR$SURV = exp(-(NEAR$F+NEAR$M_G)) #SURV is the fraction surviving at each age
      
      NEAR$Survship = 0 # set up column for survivorship (amount or fraction present at age)
      NEAR$Survship[1] = 1
      
      for(h in 1:(nrow(NEAR)-1)){ # step through ages to calc survivorship
        NEAR$Survship[h+1] = NEAR$Survship[h]*NEAR$SURV[h] #amount present at age
      } # closes survivorship loop
    } # closes age loop
    
    LEPdf[,g] = NEAR[,2]*NEAR[,3]*littlek*NEAR[,8]
    
  } # closes fishing level loop
  return(LEPdf)
} # closes function to calculate LSB


