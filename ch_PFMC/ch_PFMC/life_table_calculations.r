# Run the simulation model
# by: mikaela provost
# last edited: Feb 19, 2019
# ===================================================================

library(tidyr)
library(dplyr)
library(gridExtra)
# ---
# load the simulation model
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/simulation_model_cod_v3.r")
# load functions
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/2_cod_functions.r")
# load parms df
parms = read.csv(".csv",
                      header=TRUE,stringsAsFactors = FALSE)
parms = as.data.frame(eigentable)

# *************************************** #
# Plan
# (1) convert age to length (vonB & Schnute) 
# (2) convert length(cm) to weight (kg) at age 
# (3) calculate eggs at age 
# (4) calculate proportion mature at age
# (5) calculate vulnerability-at-age to fishing 
# *************************************** #
# This script calculates life table infomation needed
# to assemble Leslie matrices for different F values:
# L_a, B_a, propmat_a, eggs_a, vul_a
# maxage & M are stored in parms df
# assemble_Leslie(maxage,L_a,B_a,propmat_a,eggs_a,vul_a,M)


# *************************************** #
# (1) convert age to length (vonB & Schnute) 
La_list <- as.list(rep(NA,length=length(parms$spp)))
names(La_list) <- parms$spp #list of cm-at-age vectors

for(i in 1:length(pfmcsp)){ #for each spp
  if(parms[parms$spp == spp[i],]$growthFUN == "vonb"){ #if vonb equation
    #define the vonb parameters maxage,t0,Linf,Kvonb
    maxage = parms[parms$spp==parms$spp[i],]$maxage
    t0 = parms[parms$spp==parms$spp[i],]$tknot
    Linf = parms[parms$spp==parms$spp[i],]$Linf
    Kvonb = parms[parms$spp==parms$spp[i],]$Kvonb
    #run the vonb model to get cm per age
    ages <- seq(from=1,to=maxage,by=1)
    L_a <- vonbertgrowth(maxage=maxage,t0=t0,Linf=Linf,Kvonb=Kvonb,ages=ages)
    
  } else {
    #define schnute parameters maxage,t1,t2,L1,L2,Ksch
    maxage = parms[parms$spp==spp[i],]$maxage
    t1 = parms[parms$spp==parms$spp[i],]$t1
    t2 = parms[parms$spp==parms$spp[i],]$t2
    L1 = parms[parms$spp==parms$spp[i],]$L1
    L2 = parms[parms$spp==parms$spp[i],]$L2
    Ksch = parms[parms$spp==parms$spp[i],]$Ksch
    #run the schnute model to get cm per age
    ages <- seq(from=1,to=maxage,by=1)
    L_a <- schnutegrowth(maxage=maxage,t1=t1,t2=t2,L1=L1,L2=L2,Ksch=Ksch,ages=ages)}
  La_list[[i]] <- L_a
  rm(maxage,ages,t1,t2,L1,L2,Ksch)}
rm(i)


# *************************************** #
# (2) convert length(cm) to weight (kg) at age 
Ba_list <- as.list(rep(NA,length=length(parms$spp)))
names(Ba_list) <- parms$spp #list of kg-at-age vectors 
for(i in 1:length(parms$spp)){
  Ba_list[[i]] <- calc_wt_at_age(
    L_a=La_list[[i]],
    cmkga=parms[parms$spp == parms$spp[i],]$cmkga,
    cmkgb=parms[parms$spp == parms$spp[i],]$cmkgb) }
rm(i,L_a,cmkga,cmkgb)


# *************************************** #
# (3) calculate eggs at age 
eggs_list <- as.list(rep(NA,length=length(parms$spp)))
names(eggs_list) <- parms$spp
for(i in 1:length(parms$spp)){
  if(parms[parms$spp == parms$spp[i],]$eggsFUN == "eggsFUN1") {
    #define parms for fecFUN1 B_a,eggs_intercept,eggs_slope
    B_a = Ba_list[[i]]
    intercept = parms[parms$spp==parms$spp[i],]$eggs_intercept
    slope = parms[parms$spp==parms$spp[i],]$eggs_slope
    eggs <- eggsFUN1(B_a=B_a,eggs_intercept=eggs_intercept,eggs_slope=eggs_slope)
  } else {
    #define parms for fecFUN2
    B_a = Ba_list[[i]] 
    eggs_intercept = parms[parms$spp==parms$spp[i],]$eggs_intercept
    eggs_slope = parms[parms$spp==parms$spp[i],]$eggs_slope
    eggs <- eggsFUN2(B_a=B_a,intercept=eggs_intercept,slope=eggs_slope)
  }
  eggs_list[[i]] <- eggs #store vector of eggs in list
  rm(B_a,intercept,slope,eggs)
}
rm(i)


# *************************************** #
# (4) calculate proportion mature at age
prop_list <- is.list(rep(NA,length=length(parms$spp)))
for(i in 1:length(parms$spp)){
  #define parameters for prop mature function
  maxage=parms[parms$spp==parms$spp[i],]$maxage
  slope=parms[parms$spp==parms$spp[i],]$maturityslope
  inflection=parms[parms$spp==parms$spp[i],]$maturityinflectionAGE
  #calculate proportion mature at age, store in list
  prop_list[[i]] <- calc_mat_at_age(maxage,slope,inflection)
}
rm(i,maxage,slope,inflection)


# *************************************** #
# (5) calculate vulnerability-at-age to fishing 
# for now, using proportion mature
vul_list <- prop_list
rm(i,maxage,slope,inflection)
