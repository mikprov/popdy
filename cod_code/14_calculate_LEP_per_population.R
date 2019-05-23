# Calculate Lifetime Egg Production for each population
# By: Mikaela Provost
# Last edited: Jan 17, 2019

# LEP = sum of l_suba times m_suba (see lecture 11 notes in WFC 122)
# l_suba is calculated as 1 for age 1, s for age 2, s^2 for age 3, etc

# load functions
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/2_cod_functions.r")
# load cod data
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/0_load_cod_data.r") 

# load functions, libraries
library(ggplot2)
library(gridExtra)
library(dplyr)
library(plyr)
library(ggrepel)
library(devtools)
library(broom)

F.halfmax = seq(from=0,to=3.0,by=0.01) #for now, F ranges from 0 to 1
F.halfmax = 0
tknot =0 #used in vonB eq. which is in calculate_LSB_at_age_by_F() function


LSBlist <- as.list(rep(NA,length(codNames))) #LSB at age for pops stored here
LEP_at_f_levels <- as.list(rep(NA,length(codNames)))
names(LEP_at_f_levels) <- codNames
names(LSBlist) <- codNames

for (i in 1:length(LSBlist)) { # for each population
  # load parms for cod pop i
  source(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_pops/',names(LSBlist)[i], '.r', sep=''))
  # this should load parms: L_inf, K, TEMP, maxage, B0, B1, tknot=0
  
  # for each pop calc LSB at age (note: k is not incorporated here)
  lsb = calculate_LSB_at_age_by_F(maxage=maxage,
                                  L_inf=L_inf, K=K, TEMP=TEMP, F.halfmax=F.halfmax,
                                  B0=B0,B1=B1)
  LSBlist[[i]] = lsb #LSBlist should have length of codNames
  LEP <- colSums(lsb)
  LEP_at_f_levels[[i]] = cbind(F.halfmax,LEP)
  #LSB = wt at age * prop mature at age * survival at age
  rm(B0,B1,K,L_inf,TEMP)
}
rm(i,lsb,LEP)

# plot: Egg production over age when F=0 (save for now)
par(mfrow=c(4,4))
for(i in 1:length(LSBlist)){
  plot(x=seq(from=1,to=length(LSBlist[[i]][,1])), y=LSBlist[[i]][,1],xlab="age", 
       ylab="egg production", main=codNames[i],type="l")
}
par(mfrow=c(1,1))

# Calculate 1/LEP for unfished biomass
#unfishedLEP <- rep(NA,length=length(codNames))
#LEP35 <- rep(NA,length=length(codNames))
#oneoverLEP <- rep(NA,length=length(codNames))
#oneoverLEP35 <- rep(NA,length=length(codNames))
LEPinfo <- matrix(NA,nrow=length(codNames),ncol=4)
for(i in 1:length(codNames)){
  p <- as.data.frame(LEP_at_f_levels[[i]]) #2->i
  LEPinfo[i,1] <- p$LEP[1] #unfished LEP (col1)
  LEPinfo[i,2] <- LEPinfo[i,1]*0.35 #35% of unfished LEP (col2)
  LEPinfo[i,3] <- 1/LEPinfo[i,1] #slope of unfished LEP (col3)
  LEPinfo[i,4] <- 1/LEPinfo[i,2] #slope of 35% LEP, this will become my alpha
}
# merge LEP vectors into one df
LEPinfo <- as.data.frame(LEPinfo)
names(LEPinfo) <- c("unfishedLEP","LEP35","oneoverLEP","oneoverLEP35")
LEPinfo <- cbind(codNames,LEPinfo)

alphas <- LEPinfo$oneoverLEP35
# one plot for each pop, plotting both oneoverLEP and oneoverLEP35
par(mfrow=c(4,4))
for(i in 1:length(codNames)){
plot(1, type="n", xlab="", ylab="", xlim=c(0, 8), ylim=c(0, 8),
     main=LEPinfo$codNames[i]) #empty plot
abline(a=0,b=LEPinfo$oneoverLEP[i],col="blue")
abline(a=0,b=LEPinfo$oneoverLEP35[i],col="red")
text(x=c(7,1),y=c(7,7), col=c("blue","red"),
     label=c(round(LEPinfo$oneoverLEP[i],2),
             round(LEPinfo$oneoverLEP35[i],2)) )
}
