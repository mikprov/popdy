# test this:
# are these the same - 
# a) p_spawn at age = (LSB at age)/sum(LSB across ages)
# b) p_spawn at age = (survivorship at age*maturity at age)/(sum of survivorship*maturity across ages)

cod = read.table("C:/Users/provo/Documents/GitHub/popdy/cod_code/huiyu/cod_all_2013.txt",header=T,na.strings="NA")
setwd("C:/Users/provo/Documents/GitHub/popdy/cod_code/huiyu/LSB_varyfishing")

NorthseaD = subset(cod, AREA=="NORTH_SEA") # test population 

#calculate_LSB_at_age_by_F <- function(data,codPopname,littlek){
  
  NorthseaD = subset(NorthseaD,Yearclass>1959)  #yearclasses vary among stocks
  NorthseaD = subset(NorthseaD,Yearclass<1990)
  
  # -- load parms for cod pop:
  source(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_pops/',"Northsea", '.r', sep=''))
  # this should load parms: L_inf, K, TEMP
  
  # -- calculate maturity, weight parms
  (mod.mat = glm(MATPROP~AGE,family=binomial,data=NorthseaD)) #gives maturity beta coefficients
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
    
    LEPdf[,g] = NEAR[,2]*NEAR[,3]*1*NEAR[,8]
    
  } # closes fishing level loop
  #return(list(LEPdf=LEPdf,NEAR=NEAR))
  
#} # closes function to calculate LSB

# prob of spawn test #1
p_spawnA = LEPdf[1:40,301]/sum(LEPdf[1:40,301])
plot(Age,p_spawnA)
# prob of spawn test #2
# note: since the NEAR matrix is made in the for loop through fishing levels, 
# the NEAR matrix that I'm referencing to corresponds to the highest fishing level
# which is why I'm calling column 301 in p_spawnA
p_spawnB = (NEAR$mat1*NEAR$Survship*NEAR$growth)/sum(NEAR$mat1*NEAR$Survship*NEAR$growth)
plot(Age,p_spawnB)

ns_test <- calculate_LSB_at_age_by_F(data=NorthseaD, codPopname="Northsea", littlek=1)
