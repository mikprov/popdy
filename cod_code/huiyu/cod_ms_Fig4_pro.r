#install.packages(faraway)
library("faraway")
library("dplyr")


cod = read.table("C:/Users/provo/Documents/GitHub/popdy/cod_code/huiyu/cod_all_2013.txt",header=T,na.strings="NA")
setwd("C:/Users/provo/Documents/GitHub/popdy/cod_code/huiyu/LSB_varyfishing")

Northsea = subset(cod, AREA=="NORTH_SEA")
Coas = subset(cod, AREA=="COAS")
E_Baltic = subset(cod, AREA=="E_BALTIC")
W_Baltic = subset(cod, AREA=="W_BALTIC")
Faroe = subset(cod, AREA=="FAROE")
NE_Arctic = subset(cod, AREA=="NE_ARCTIC")
Celtic = subset(cod, AREA=="CELTIC_SEA")
Iceland = subset(cod, AREA=="ICELAND")
Irish = subset(cod, AREA=="IRISH_SEA")
Kat = subset(cod, AREA=="KATTEGAT")
W_Scotland = subset(cod, AREA=="W_SCOTLAND")

NGulf = subset(cod, AREA =="N_GulfSL")
SGulf = subset(cod, AREA =="S_GulgSL")
GB = subset(cod, AREA =="GB")
GM = subset(cod, AREA =="GM")
cod3M = subset(cod,AREA=="cod3M")
cod3NO = subset(cod,AREA=="cod3NO")
cod3Ps = subset(cod,AREA=="cod3Ps")
cod4X = subset(cod,AREA=="cod4X")
cod2J3KL = subset(cod,AREA=="cod2J3KL")


##Model FLSB for individual stocks 
 # run this code for each individual population and 
 # save as a new text file (make sure to change name)

# ----
 # code in this section is setting up Hui-Yu's original 
 # code to run in a loop so that I don't have to run each
 # individual population separately ***in progress***
for (i in 1:length(namesvector))
Northsea = Celtic

# ----
#xtabs(y~Yearclass+AGE,Northsea)

data = subset(Northsea,Yearclass>1959)  #yearclasses vary among stocks
data = subset(data,Yearclass<1990)
Linf = 99 #77  #Table 2 of ms
K = 0.390 #0.26
to = 0

(mod.mat = glm(MATPROP~AGE,family=binomial,data=data)) #gives maturity beta coefficients
Age = 1:40  #max age = 40 yrs

L=Linf*(1-exp(-K*Age))
#MG=exp(0.55-1.61*log(L)+1.44*log(Linf)+log(K)) #Gislason model II
MG3 = exp(15.11-1.59*log(L)+0.82*log(Linf)-3891/(273.15+6.75))  #Gisllason model III
MP = 10^(-0.0066-0.279*log10(Linf)+0.6543*log10(K)+0.4634*log10(10.56))  #Pauly model
#Vul1 = data$CANUM/data$STNUM
growth = 0.00001*(Linf*(1-exp(-K*(Age-to))))^3
mat1 = ilogit(mod.mat$coef[1]+mod.mat$coef[2]*Age)


NEAR = data.frame(cbind(Age,mat1,growth))  #life table

NEAR$Vul1 = mat1 #use maturity ogive for selectivity
NEAR$M_G= 0.19+0.058*11.9   #regression fit of Fig. 5c
F.halfmax = seq(0,3,by=0.01)


LEP_Northsea2 = matrix(0,length(Age),length(F.halfmax))
for(i in 1:length(F.halfmax)){
  for(j in 1:length(Age)){
    NEAR$F[j] = NEAR$Vul1[j]*F.halfmax[i]	#Vul1 should be selectivity, but it's mat
    NEAR$SURV = exp(-(NEAR$F+NEAR$M_G)) #SURV is the fraction surviving at each age
    
    
    NEAR$Survship = 0
    NEAR$Survship[1] = 1
    for(k in 1:(nrow(NEAR)-1)){
      NEAR$Survship[k+1] = NEAR$Survship[k]*NEAR$SURV[k] #amount present at age
    }
  }
  #LEP_Northsea2[i] = sum(NEAR[,2]*NEAR[,3]*NEAR[,8])
  LEP_Northsea2[,i] = NEAR[,2]*NEAR[,3]*NEAR[,8]
}

write.table(LEP_Northsea2,"Slo_Celtic1_new.txt")  

# NEAR has columns for survivorship, F
#  let's plot some of these to make sure they look right
par(mfrow=c(3,2))
plot(x=NEAR$Age, y=NEAR$mat1, main="prop mature @ age")
plot(x=NEAR$Age, y=NEAR$growth, main="wt @ age")
plot(x=NEAR$Age, y=NEAR$F, main="F @ age")
plot(x=NEAR$Age, y=NEAR$SURV, main="fraction surviving @ age")
plot(x=NEAR$Age, y=NEAR$Survship, main="survivorship @ age")
par(mfrow=c(1,1))
# the text files being read in -- Hui-Yu gave to me
# I can recreat them from the code above 

#Figure4_new
# each df belongs to one population
# rows = ages
# cols = fishing rates
# each element = the amount of spawning biomass
# produced at that age with that fishing rate
# to get lifetime spawning biomass, get sum for each column
data1 = read.table("Slo_Celtic1.txt",header=T)
data1.1 = read.table("Slo_Celtic1_new.txt",header=T)
test = cbind(colSums(data1),colSums(data1.1))
datasource= c("from Hui-Yu","Mikaela created")
matplot(x=F.halfmax,y=test, type="l",ylab="LSB")
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
