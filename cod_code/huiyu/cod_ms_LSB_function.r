rm(list=ls())

#install.packages(faraway)
library("faraway")
library("dplyr")


cod = read.table("C:/Users/provo/Documents/GitHub/popdy/cod_code/huiyu/cod_all_2013.txt",header=T,na.strings="NA")
setwd("C:/Users/provo/Documents/GitHub/popdy/cod_code/huiyu/LSB_varyfishing")

NorthseaD = subset(cod, AREA=="NORTH_SEA")
CoasD = subset(cod, AREA=="COAS")
E_BalticD = subset(cod, AREA=="E_BALTIC")
W_BalticD = subset(cod, AREA=="W_BALTIC")
FaroeD = subset(cod, AREA=="FAROE")
NE_ArcticD = subset(cod, AREA=="NE_ARCTIC")
CelticD = subset(cod, AREA=="CELTIC_SEA")
IcelandD = subset(cod, AREA=="ICELAND")
IrishD = subset(cod, AREA=="IRISH_SEA")
KatD = subset(cod, AREA=="KATTEGAT")
W_ScotlandD = subset(cod, AREA=="W_SCOTLAND")

NGulfD = subset(cod, AREA =="N_GulfSL")
SGulfD = subset(cod, AREA =="S_GulgSL")
GBD = subset(cod, AREA =="GB")
GMD = subset(cod, AREA =="GM")
cod3MD = subset(cod,AREA=="cod3M")
cod3NOD = subset(cod,AREA=="cod3NO")
cod3PsD = subset(cod,AREA=="cod3Ps")
cod4XD = subset(cod,AREA=="cod4X")
cod2J3KLD = subset(cod,AREA=="cod2J3KL")


datalist <- list(NorthseaD,CoasD,W_BalticD,
                 FaroeD,NE_ArcticD,CelticD,
                 IcelandD,KatD,W_ScotlandD,
                 NGulfD,GBD,GMD,
                 cod3NOD,cod3MD,cod2J3KLD,
                 cod3PsD)

codNames <- c("Northsea","Coas","W_Baltic",
              "Faroe","NE_Arctic","Celtic",
              "Iceland","Kat","W_Scotland",
              "NGulf","GB","GM",
              "3NO","3M","2J3KL",
              "3Ps")

names(datalist) <- c("Northsea","Coas","W_Baltic",
                     "Faroe","NE_Arctic","Celtic",
                     "Iceland","Kat","W_Scotland",
                     "NGulf","GB","GM",
                     "3NO","3M","2J3KL",
                     "3Ps")

##Model FLSB for individual stocks 
 # run this code for each individual population and 
 # save as a new text file (make sure to change name)

# ----
 # code in this section is setting up Hui-Yu's original 
 # code to run in a loop so that I don't have to run each
 # individual population separately ***in progress***


Age = 1:40  #max age = 40 yrs
to = 0
F.halfmax = seq(0,3,by=0.01)


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
  NEAR = data.frame(cbind(Age,mat1,growth))  #life table
  NEAR$Vul1 = mat1 #use maturity ogive for selectivity
  NEAR$M_G= 0.19+0.058*TEMP   #regression fit of Fig. 5c, loaded in parms set
  
  # -- 
  LEPdf = matrix(0,length(Age),length(F.halfmax))
    for(g in 1:length(F.halfmax)){ # step through each fishing level
       for(j in 1:length(Age)){ # step through ages
          NEAR$F[j] = NEAR$Vul1[j]*F.halfmax[g]	# Vul1 should be selectivity, but it's mat
                                                # calculate F = selectivity * F rate
          NEAR$SURV = exp(-(NEAR$F+NEAR$M_G)) #SURV is the fraction surviving at each age
    
          NEAR$Survship = 0 # set up column for survivorship (amount or fraction present at age)
          NEAR$Survship[1] = 1
    
          for(k in 1:(nrow(NEAR)-1)){ # step through ages to calc survivorship
            NEAR$Survship[k+1] = NEAR$Survship[k]*NEAR$SURV[k] #amount present at age
            } # closes survivorship loop
          } # closes age loop
  
      LEPdf[,g] = NEAR[,2]*NEAR[,3]*littlek*NEAR[,8]
      
    } # closes fishing level loop
  return(LEPdf)
} # closes function to calculate LSB




# --- Export LSB tables when k=1 --- #
# loop over pop data in datalist (littlek = 1)
# prep for loop:
LSBlist = as.list(rep(NA,length(datalist))) # store matrix of LSB vs F
names(LSBlist) = codNames
for (i in 1:length(datalist)) { # step through each pop in datalist
  LSBlist[[i]] = calculate_LSB_at_age_by_F(data=datalist[[i]], codPopname=codNames[i], littlek=1)
    # store the output from calculate_LSB.. function in list. output is LSB by age and F
  write.table(LSBlist[[i]],file=paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLSB/k1/",codNames[i],".txt",sep=""))
    # export txt file of LSB
  
}

  
# --- Export LSB tables when k=0.2 --- #
# loop over pop data in datalist 
# prep for loop:
LSBlist = as.list(rep(NA,length(datalist))) # store matrix of LSB vs F
names(LSBlist) = codNames
for (i in 1:length(datalist)) { # step through each pop in datalist
  LSBlist[[i]] = calculate_LSB_at_age_by_F(data=datalist[[i]], codPopname=codNames[i], littlek=0.2)
  # store the output from calculate_LSB.. function in list. output is LSB by age and F
  write.table(LSBlist[[i]],file=paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLSB/k0.2/",codNames[i],".txt",sep=""))
  # export txt file of LSB
  
}



# --- Export LSB tables when k=0.4 --- #
# loop over pop data in datalist 
# prep for loop:
LSBlist = as.list(rep(NA,length(datalist))) # store matrix of LSB vs F
names(LSBlist) = codNames
for (i in 1:length(datalist)) { # step through each pop in datalist
  LSBlist[[i]] = calculate_LSB_at_age_by_F(data=datalist[[i]], codPopname=codNames[i], littlek=0.4)
  # store the output from calculate_LSB.. function in list. output is LSB by age and F
  write.table(LSBlist[[i]],file=paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLSB/k0.4/",codNames[i],".txt",sep=""))
  # export txt file of LSB
  
}


# --- Export LSB tables when k=0.6 --- #
# loop over pop data in datalist 
# prep for loop:
LSBlist = as.list(rep(NA,length(datalist))) # store matrix of LSB vs F
names(LSBlist) = codNames
for (i in 1:length(datalist)) { # step through each pop in datalist
  LSBlist[[i]] = calculate_LSB_at_age_by_F(data=datalist[[i]], codPopname=codNames[i], littlek=0.6)
  # store the output from calculate_LSB.. function in list. output is LSB by age and F
  write.table(LSBlist[[i]],file=paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLSB/k0.6/",codNames[i],".txt",sep=""))
  # export txt file of LSB
  
}



# --- Export LSB tables when k=0.8 --- #
# loop over pop data in datalist 
# prep for loop:
LSBlist = as.list(rep(NA,length(datalist))) # store matrix of LSB vs F
names(LSBlist) = codNames
for (i in 1:length(datalist)) { # step through each pop in datalist
  LSBlist[[i]] = calculate_LSB_at_age_by_F(data=datalist[[i]], codPopname=codNames[i], littlek=0.8)
  # store the output from calculate_LSB.. function in list. output is LSB by age and F
  write.table(LSBlist[[i]],file=paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLSB/k0.8/",codNames[i],".txt",sep=""))
  # export txt file of LSB
  
}



# --- Export LSB tables when k=0 --- #
# loop over pop data in datalist 
# prep for loop:
LSBlist = as.list(rep(NA,length(datalist))) # store matrix of LSB vs F
names(LSBlist) = codNames
for (i in 1:length(datalist)) { # step through each pop in datalist
  LSBlist[[i]] = calculate_LSB_at_age_by_F(data=datalist[[i]], codPopname=codNames[i], littlek=0)
  # store the output from calculate_LSB.. function in list. output is LSB by age and F
  write.table(LSBlist[[i]],file=paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLSB/k0/",codNames[i],".txt",sep=""))
  # export txt file of LSB
  
}












# --- plot my LSB vs F with Hui-Yu's plot:
#Figure4_new
# each df belongs to one population
# rows = ages
# cols = fishing rates
# each element = the amount of spawning biomass
# produced at that age with that fishing rate
# to get lifetime spawning biomass, get sum for each column

# plot my calculations vs Hui-Yu LSB
HYfilenames <- c("Northsea","Coas","W_Baltic",
              "Faroe","NE_Arctic","Celtic",
              "Iceland","Irish","Kat","W_Scotland",
              "NGulf","SGulf","GB","GM","3NO")
datalisttest <- datalist[HYfilenames]
names(datalisttest)


par(mfrow=c(4,3))
for(i in 1:length(HYfilenames)){

  dataHY = read.table(paste('Slo_',HYfilenames[i],'1', '.txt', sep=''),header=T)
  dataMP = LSBlist[[i]]
  test = cbind(colSums(dataHY),colSums(dataMP))
  datasource= c("Hui-Yu","Mikaela")
  matplot(x=F.halfmax, y=test, type="l",ylab="LSB",
          main=paste(HYfilenames[i]),xlab="F")
  legend("topright",legend=datasource,col=1:length(datasource),lty=1:length(datasource),lwd=1.2)

}
par(mfrow=c(1,1))

# read in Hui-Yu's calculations for plotting
data1 = read.table("Slo_Celtic1.txt",header=T)
data1.1 = LSBlist[["Celtic"]]
test = cbind(colSums(data1),colSums(data1.1))
datasource= c("from Hui-Yu","Mikaela created")
matplot(x=F.halfmax, y=test, type="l",ylab="LSB")
legend("topright",legend=datasource,col=1:length(datasource),lty=1:length(datasource),lwd=1.2)

data2 = read.table("Slo_Irish1.txt",header=T)
data3 = read.table("Slo_W_Scotland1.txt",header=T)
data4 = read.table("Slo_Northsea1.txt",header=T)
data5 = read.table("Slo_GB1.txt",header=T)
data6 = read.table("Slo_GM1.txt",header=T)
data7 = read.table("Slo_Faroe1.txt",header=T)
data8 = read.table("Slo_W_Baltic1.txt",header=T)
data9 = read.table("Slo_Coas1.txt",header=T)
data10 = read.table("Slo_Kat1.txt",header=T)
data11 = read.table("Slo_Iceland1.txt",header=T)
data12 = read.table("Slo_NE_Arctic1.txt",header=T)
data13 = read.table("Slo_3M1.txt",header=T)
data14 = read.table("Slo_3Ps1.txt",header=T)
data15 = read.table("Slo_3NO1.txt",header=T)
data16 = read.table("Slo_NGulf1.txt",header=T)
data17 = read.table("Slo_2J3KL1.txt",header=T)


FLSB = matrix(0,301,18)
# sum each column, add take the ln(), then divide each 
# column sum by the sum of the first column. The first
# column corresponds to when F=0. Every other column 
# after that corresponds to increasing F rates.
FLSB[,1]=log(apply(data1,2,sum)/apply(data1,2,sum)[1])
FLSB[,2]=log(apply(data2,2,sum)/apply(data2,2,sum)[1])
FLSB[,3]=log(apply(data3,2,sum)/apply(data3,2,sum)[1])
FLSB[,4]=log(apply(data4,2,sum)/apply(data4,2,sum)[1])
FLSB[,5]=log(apply(data5,2,sum)/apply(data5,2,sum)[1])
FLSB[,6]=log(apply(data6,2,sum)/apply(data6,2,sum)[1])
FLSB[,7]=log(apply(data7,2,sum)/apply(data7,2,sum)[1])
FLSB[,8]=log(apply(data8,2,sum)/apply(data8,2,sum)[1])
FLSB[,9]=log(apply(data9,2,sum)/apply(data9,2,sum)[1])
FLSB[,10]=log(apply(data10,2,sum)/apply(data10,2,sum)[1])
FLSB[,11]=log(apply(data11,2,sum)/apply(data11,2,sum)[1])
FLSB[,12]=log(apply(data12,2,sum)/apply(data12,2,sum)[1])
FLSB[,13]=log(apply(data13,2,sum)/apply(data13,2,sum)[1])
FLSB[,14]=log(apply(data14,2,sum)/apply(data14,2,sum)[1])
FLSB[,15]=log(apply(data15,2,sum)/apply(data15,2,sum)[1])
FLSB[,16]=log(apply(data16,2,sum)/apply(data16,2,sum)[1])
FLSB[,17]=log(apply(data17,2,sum)/apply(data17,2,sum)[1])
FLSB[,18]=seq(0,3,by=0.01)

FLSB1 = FLSB[c(1:201),]
FLSB = FLSB1
plot(FLSB[,1]~FLSB[,18], type="l",lwd=2,col=2,ylim=c(-2.6,0),xlab="F (year-1)",ylab="ln(FLSB)")
for(i in 1:16){
  lines(FLSB[,i+1]~FLSB[,18],lwd=2,lty=i+1,col=i+1)
}



for(i in 1:3){
  lines(FLSB[,1+i]~FLSB[,18],lwd=2,lty=i+1,col=2)
}
for(i in 4:5){
  lines(FLSB[,1+i]~FLSB[,18],col="dark grey",lwd=2,lty=i-3)
}
for(i in 6:8){
  lines(FLSB[,1+i]~FLSB[,18],lwd=2,lty=i-5,col=1)
}
for(i in 9:13){
  lines(FLSB[,1+i]~FLSB[,18],lwd=2,lty=i-8,col=8)
}
for(i in 14:16){
  lines(FLSB[,1+i]~FLSB[,18],lwd=2,lty=i-13,col=4)
  
}
abline(h=log(0.35))

#abline(h=log(0.1),lty=2)
legend("bottomright",legend=c("Celtic","Irish","W Scotland","North Sea","GB","GM","Faroes","W Baltic","Coas N","Kattegat","Iceland","NE Arctic","3M","3Ps","3NO","N Gulf","2J3KL"),col=c(2,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),lwd=2,lty=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),cex=0.65)



# ----------------------------------------------------------
# Re-create spawning biomass distributions from LSB txt files

# prep for plotting
codpop_name <- c("Celtic_Sea","Irish_Sea","West_Scotland","North_Sea",
             "Georgoes_Bank","Gulf_of_Maine","Faroes",
             "West_Baltic","Kattegat", "Iceland","NE_Arctic",
             "3M","3NO","N_Gulf_St_Lawrence")
codpops1 <- c("Slo_Celtic1","Slo_Irish1","Slo_W_Scotland1","Slo_Northsea1",
              "Slo_GB1","Slo_GM1","Slo_Faroe1","Slo_W_Baltic1",
              "Slo_Kat1","Slo_Iceland1","Slo_NE_Arctic1","Slo_3M1",
              "Slo_3NO1","Slo_NGulf1")

Fs <- paste("F",F.halfmax, sep="") 
cod.no <- 1:length(codpops1)
SBdfs_list = as.list(rep(NA,length(codpops1)))

# multiplot panel of spawning biomass:
pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_code/huiyu/LSB_varyfishing/SBplot.pdf', width=7, height=10)
par(mfrow=c(5,4))
for (i in seq_along(codpops1)) { 
  datax <- read.table(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_code/huiyu/LSB_varyfishing/'
                 ,codpops1[i], '.txt', sep=''),header=T)
  
  datax = rbind(Fs,datax) #add F labels to df as first row
  colnames(datax) = datax[1,] #change the first row to column name
  datax = datax[-1,] #remove first row (bc now it's the header)
  dataxx = subset(datax, select=c("F0","F0.5","F1")) #only plot for these values of F
  SBdfs_list[[i]] = datax #store data frame with fishing rate as header in a list
  matplot(x=Age, y=dataxx, type="l",ylim=c(0,0.30),main=codpop_name[i],ylab="",xlab="")

}
par(mfrow=c(1,1))
dev.off()
