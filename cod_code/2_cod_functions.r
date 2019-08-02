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

tknot = 0
assemble_Leslie <- function(B0,B1,maxage,K,L_inf,TEMP,F.halfmax,tknot) {
                            
  Age=1:maxage
  
  #data = subset(data,Yearclass>1959)  #yearclasses vary among stocks
  #data = subset(data,Yearclass<1990)
  
  # -- calculate maturity, weight parms
  #mod.mat = glm(MATPROP~AGE,family=binomial,data=data) #gives maturity beta coefficients
  #L=L_inf*(1-exp(-K*Age)) #length at age
  #MG=exp(0.55-1.61*log(L)+1.44*log(L_inf)+log(K)) #Gislason model II
  #MG3 = exp(15.11-1.59*log(L)+0.82*log(L_inf)-3891/(273.15+6.75))  #Gisllason model II (mortality)
  #MP = 10^(-0.0066-0.279*log10(L_inf)+0.6543*log10(K)+0.4634*log10(10.56))  #Pauly model (mortality)
  #mat1 = ilogit(mod.mat$coef[1]+mod.mat$coef[2]*Age) # proportion mature at age
  
  growth = 0.00001*(L_inf*(1-exp(-K*(Age-tknot))))^3 # weight at age, uses vonB
  # NOTE: allometric length-wt scale is 0.00001 b/c Wang used it.
  # looked up on Fishbase, cod stocks range 0.005-0.0117 (cm-grams), geometric average=0.0071
  # units: 'growth' at age is weight at age in kilograms
  mat1 = ilogit(B0 + B1*Age) #B0 and B1 are published in Wang et al
  
  # -- assemble NEAR df: use NEAR to calculate Leslie matrix
  NEAR = data.frame(cbind(Age,mat1,growth))  #life table: maturity & size at age
  NEAR$Vul1 = mat1 #insert selectivity at age, use maturity ogive for selectivity
  #NEAR$M_G= 0.19+0.058*TEMP #insert mortality at age, from regression fit of Fig. 5c, 
 
  # ** CHOOSE EITHER PAULY OR GIASLON ** #
  NEAR$M = rep(MG,length=length(NEAR[,1]))
  #NEAR$M = rep(Mp,length=length(NEAR[,1]))

  A = matrix(0,length(Age),length(Age)) #empty Leslie matrix - Giaslon
  
  # -- for each age, get F and survival
  for(j in 1:length(Age)){ # step through ages
    NEAR$FISH[j] = NEAR$Vul1[j]*F.halfmax	# Vul1 should be selectivity, but we are using mat (maturity)
    NEAR$SURV = exp(-(NEAR$FISH+NEAR$M)) #SURV is the fraction surviving at each age
    
    NEAR$Survship = 0 # set up column for survivorship (amount or fraction present at age)
    NEAR$Survship[1] = 1
    
    for(k in 1:(nrow(NEAR)-1)){ # step through ages to calc survivorship using MG
      NEAR$Survship[k+1] = NEAR$Survship[k]*NEAR$SURV[k] #amount present at age
    }
    
  }
  NEAR$egg_production <- NEAR$growth*NEAR$mat1*NEAR$Survship
  
  A[1,] = NEAR$mat1*NEAR$growth# insert fecundity (maturity * weight) in top row
  # ----------------------------------------------------------------------- #
  for(u in 2:length(Age)-1){ # insert survival on subdiagonal
    A[u+1,u]=NEAR$SURV[u]  
  }
  return(list(A=A,NEAR=NEAR))
  
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
  # eigen value has no imaginary part. 
  a.sq = (Re(ev$values[1]))^2 # square the real part of 1st eigenvalue
  #b.sq = (Im(ev$values[1]))^2 # square the imaginary part of 1st eigenvalue
  firstval = sqrt(a.sq)# + b.sq) # magnitude of 1st eigenvalue
  return(firstval)
  rm(a.sq) # remove from workpace
  #rm(b.sq) # remove from workpace
}


extract_second_eigen_value <- function(Lesliematrix){
  # get magnitude of second eigenvalue
  # think of real value on x-axis, imaginary value on y-axis
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
# littlek = the value multiplied across the top row of the Leslie matrix, 
#   is a measure of fishing pressure (Lauren's analysis)

calculate_LSB_at_age_by_F <- function(B0,B1,maxage,L_inf,K,TEMP,F.halfmax){
  
  Age = 1:maxage
  
  #data = subset(data,Yearclass>1959)  #yearclasses vary among stocks, 
  #data = subset(data,Yearclass<1990)  #based on Hu-Yui code
  
  # -- calculate maturity, weight parms
  #mod.mat = glm(MATPROP~AGE,family=binomial,data=data) #gives maturity beta coefficients
  #L=L_inf*(1-exp(-K*Age)) #calculate length at age, vonB
  #Hu-Yui calculated 3 equations for mortality, don't need them for now
  #MG=exp(0.55-1.61*log(L)+1.44*log(L_inf)+log(K)) #Gislason model II
  #MG3 = exp(15.11-1.59*log(L)+0.82*log(L_inf)-3891/(273.15+6.75))  #Gisllason model III
  #MP = 10^(-0.0066-0.279*log10(L_inf)+0.6543*log10(K)+0.4634*log10(10.56))  #Pauly model
  #Vul1 = data$CANUM/data$STNUM
  #mat1 = ilogit(mod.mat$coef[1]+mod.mat$coef[2]*Age) 
    #above, I was using coefficients I calculated, but Wang et al.
    #published them, so just using those (B0 and B1, below)
  growth = 0.00001*(L_inf*(1-exp(-K*(Age-tknot))))^3 #weight at age (kilograms)
  mat1 = ilogit(B0 + B1*Age) #B0 and B1 are published in Wang et al
  
  # -- assemble NEAR df: use NEAR to calculate LSB
  NEAR = data.frame(cbind(Age,mat1,growth))  #cols: age, prop maturity at age, growth (wt at age)
  NEAR$Vul1 = mat1 #use maturity ogive for selectivity ogive -- see Wang et al for justification
  #NEAR$M_G= 0.19+0.058*TEMP   #mortality at age, regression fit of Fig. 5c, 'TEMP' is loaded in parms set
  NEAR$MG = rep(MG,length=length(NEAR[,1]))
  NEAR$MP = rep(Mp,length=length(NEAR[,1]))
  
  # -- LEPdf: rows ~ age, col ~ fishing levels (F.halfmax), observations is egg production at age
  LEPdf_MG = matrix(0,nrow=length(Age),ncol=length(F.halfmax))
  LEPdf_MP = matrix(0,nrow=length(Age),ncol=length(F.halfmax))
  
  for(g in 1:length(F.halfmax)){ # step through each fishing level, for the first column in LEPdf
   fishing_at_age = NEAR$Vul1*F.halfmax[g]	# Vul1 should be selectivity, but it's mat. F=selectivity*F rate
   NEAR$SURV_MG = exp(-(fishing_at_age+NEAR$MG)) #SURV is the fraction surviving at each age
   NEAR$SURV_MP = exp(-(fishing_at_age+NEAR$MP)) #SURV is the fraction surviving at each age
   #l_suba <- rep(NA,length=length(Age)) #create empty vector to store l_suba
   
   #for(j in 1:length(Age)){ #for each age
   #    l_suba[j] = surv_at_age[j]^(j-1) } #l_suba = survival^a starting at 0
   NEAR$SurvshipMG = 0 # set up column for survivorship (amount or fraction present at age)
   NEAR$SurvshipMG[1] = 1
   for(k in 1:(nrow(NEAR)-1)){ # step through ages to calc survivorship
     NEAR$SurvshipMG[k+1] = NEAR$SurvshipMG[k]*NEAR$SURV_MG[k] #amount present at age
   }
   
   NEAR$SurvshipMP = 0 # set up column for survivorship (amount or fraction present at age)
   NEAR$SurvshipMP[1] = 1
   for(k in 1:(nrow(NEAR)-1)){ # step through ages to calc survivorship
     NEAR$SurvshipMP[k+1] = NEAR$SurvshipMP[k]*NEAR$SURV_MP[k] #amount present at age
   }
      
   LEPdf_MG[,g] = NEAR$mat1*NEAR$growth*NEAR$SurvshipMG #prop mat at age * wt at age * survival from age 0 to age a
   LEPdf_MP[,g] = NEAR$mat1*NEAR$growth*NEAR$SurvshipMP #prop mat at age * wt at age * survival from age 0 to age a
    
  } # closes fishing level loop
  
  return(list(LEP_MG=LEPdf_MG, LEP_MP=LEPdf_MP))
} # closes function to calculate LSB



