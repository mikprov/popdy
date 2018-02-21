# Create Leslie matrix from spawning biomass files (from Hui-Yu)
# by: Mikaela Provost
# modified on: 2017-10-26
# ---
# plan:
# 1. read in LSB files - first column is LSB at age when F=0
rm(list=ls())

library(ggplot2)
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
  return(A) # returns Leslie matrix
} # closes assemble_Leslie matrix function



# prep for looping
Age = 1:40  #max age = 40 yrs
to = 0
F.halfmax = seq(0,3,by=0.01)

# -------------
k = 1
# -------------
# loop over pop data in datalist to generate Leslie matrices
Alist = as.list(rep(NA,length(datalist))) # store Leslie matrix
names(Alist) = codNames
for (i in 1:length(datalist)) {
  Alist[[i]]=assemble_Leslie(data=datalist[[i]], codPopname = codNames[i], littlek=k)
  write.table(Alist[[i]],
              file=paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/",
                         "k1","/",codNames[i],".txt",sep=""))
}


# -------------
k = 0.2
# -------------
# loop over pop data in datalist to generate Leslie matrices
Alist = as.list(rep(NA,length(datalist))) # store Leslie matrix
names(Alist) = codNames
for (i in 1:length(datalist)) {
  Alist[[i]]=assemble_Leslie(data=datalist[[i]], codPopname = codNames[i], littlek=k)
  write.table(Alist[[i]],
              file=paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/",
                         "k0.2","/",codNames[i],".txt",sep=""))
}


# -------------
k = 0.4
# -------------
# loop over pop data in datalist to generate Leslie matrices
Alist = as.list(rep(NA,length(datalist))) # store Leslie matrix
names(Alist) = codNames
for (i in 1:length(datalist)) {
  Alist[[i]]=assemble_Leslie(data=datalist[[i]], codPopname = codNames[i], littlek=k)
  write.table(Alist[[i]],
              file=paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/",
                         "k0.4","/",codNames[i],".txt",sep=""))
}


# -------------
k = 0.6
# -------------
# loop over pop data in datalist to generate Leslie matrices
Alist = as.list(rep(NA,length(datalist))) # store Leslie matrix
names(Alist) = codNames
for (i in 1:length(datalist)) {
  Alist[[i]]=assemble_Leslie(data=datalist[[i]], codPopname = codNames[i], littlek=k)
  write.table(Alist[[i]],
              file=paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/",
                         "k0.6","/",codNames[i],".txt",sep=""))
}


# -------------
k = 0.8
# -------------
# loop over pop data in datalist to generate Leslie matrices
Alist = as.list(rep(NA,length(datalist))) # store Leslie matrix
names(Alist) = codNames
for (i in 1:length(datalist)) {
  Alist[[i]]=assemble_Leslie(data=datalist[[i]], codPopname = codNames[i], littlek=k)
  write.table(Alist[[i]],
              file=paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/",
                         "k0.8","/",codNames[i],".txt",sep=""))
}


# -------------
k = 0
# -------------
# loop over pop data in datalist to generate Leslie matrices
Alist = as.list(rep(NA,length(datalist))) # store Leslie matrix
names(Alist) = codNames
for (i in 1:length(datalist)) {
  Alist[[i]]=assemble_Leslie(data=datalist[[i]], codPopname = codNames[i], littlek=k)
  write.table(Alist[[i]],
              file=paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/",
                         "k0","/",codNames[i],".txt",sep=""))
} 


# functions to calculate first and second eigen values
extract_first_eigen_value <- function(Lesliematrix){
  # get leading eigenvalue
  ev = eigen(Lesliematrix)
  # the zero times imaginary part, take the real part
  firstval = Re(ev$values[1])
  return(firstval=firstval)
}

extract_second_eigen_value <- function(Lesliematrix){
  # get leading eigenvalue
  ev = eigen(Lesliematrix)
  # the zero times imaginary part, take the real part
  secondval = Re(ev$values[2]) # I want to take the square root of the real part, look up this function to see if eigenvalue is already squared to incorporate the imaginary part
  return(secondval=secondval)
}


# Need to loop over Leslie matrices in Alist,
# for a given k value. Output from functions is
# two vectors: first and second eigenvalues for 
# each population (for given k).
# need to have 6 loops b/c 6 k values


# get eigenvalues for k=1
# first, read in Leslie table for one pop,
# then calculate eigenvalues & store these 
# in firsts and seconds. After the loop, 
# create vector of k values (constant) and 
# combine firsts, seconds, k vector.
firsts = rep(NA,length(codNames))
seconds = rep(NA,length(codNames))
for(i in 1:length(codNames)){
  A = read.table(paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/k1/",
                       codNames[i],".txt",sep=""))
  firsts[i] = extract_first_eigen_value(Lesliematrix=A)
  seconds[i] = extract_second_eigen_value(Lesliematrix=A)}


# calculate spawning biomass distribution -- new way: treat as probability distribution
# y axis = probability of spawning
# x axis = age
Age = 1:40
cvs = rep(NA,length(codNames))
pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLSB/SBplot_nofishing_probofspawning.pdf', width=7, height=10)
par(mfrow=c(5,3))
for (i in 1:length(codNames)) { # new way
  datax <- read.table(file = 
                        paste('C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLSB/k1/'
                              ,codNames[i], '.txt', sep=''),header=T)
  p_spawn = datax[,1] / sum(datax[,1])
  mean_age = sum(p_spawn * Age)
  sd_age = sqrt( sum( p_spawn*(Age-mean_age)^2))
  cv = sd_age/mean_age
  cvs[i] = cv
  plot(x=Age, y=p_spawn,type="l",
       main=codNames[i], ylab="probability of spawning",
       ylim=c(0,0.3))
  abline(v=mean_age,col="red",lty=2)
  legend("topright",c(paste("mean age =",round(x=mean_age,digits=2)),
                      paste("sd = ",round(x=sd_age,digits=2)),
                      paste("CV = ",round(x=cv,digits=2))))
  }
par(mfrow=c(1,1))
dev.off()
#for (i in 1:length(codNames)) { 
#  datax <- read.table(file = 
#                        paste('C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLSB/k1/'
#                              ,codNames[i], '.txt', sep=''),header=T)
#  cv[i] = sd(datax[,1])/mean(datax[,1])}
k = rep("k1",length(firsts))
eigenvaluesk1 <- as.data.frame(cbind(codNames,k,firsts,seconds,cvs))

# --------------------------------- #
# --- get eigenvalues for k=0.8 --- #
firsts = rep(NA,length(codNames))
seconds = rep(NA,length(codNames))
for(i in 1:length(codNames)){
  A = read.table(paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/k0.8/",
                       codNames[i],".txt",sep=""))
  firsts[i] = extract_first_eigen_value(Lesliematrix=A)
  seconds[i] = extract_second_eigen_value(Lesliematrix=A)}
# calculate spawning biomass distribution
cv = rep(NA,length(codNames))
for (i in 1:length(codNames)) { 
  datax <- read.table(file = 
                        paste('C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLSB/k0.8/'
                              ,codNames[i], '.txt', sep=''),header=T)
  cv[i] = sd(datax[,1])/mean(datax[,1])}
k = rep("k0.8",length(firsts))
eigenvaluesk0.8 <- as.data.frame(cbind(codNames,k,firsts,seconds,cv))


# -------------------------------- #
# combine eigenvalue datasets

eigenvalue <- rbind(
                    eigenvaluesk0.8,
                    eigenvaluesk1)
eigenvalue <- eigenvaluesk1
eigenvalue$firsts <- as.numeric(levels(eigenvalue$firsts))[eigenvalue$firsts]
eigenvalue$seconds <- as.numeric(levels(eigenvalue$seconds))[eigenvalue$seconds]
eigenvalue$cvs <- as.numeric(levels(eigenvalue$cvs))[eigenvalue$cvs]
eigenvalue$dampratio <- eigenvalue$seconds / eigenvalue$firsts
eigenvalue$temp = c(8.4,7.13,7,7.4,4,11.9,5.8,6.5,9.57,1,8,8,1.75,3.5,0,2.5)



ggplot(eigenvalue, aes(x=cvs,y=dampratio)) +
  geom_point(aes(color=temp), size=3) +  
  geom_text(aes(label=codNames),hjust=-0.1,vjust=-0.2,check_overlap = F) +
  scale_color_gradientn(colors=rev(rainbow(n=10,start=0,end=0.7)))
  
eigenvalue <- eigenvalue[!(eigenvalue$codNames == "Celtic"),]
# plot for WFCB presentation
pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLSB/cv_dampratio_probofsp.pdf', width=8, height=6.5)
ggplot(eigenvalue, aes(x=cvs,y=dampratio)) +
  theme_classic() +
  theme(axis.title.y=element_text(size=rel(1.5))) +
  theme(axis.title.x=element_text(size=rel(1.5))) +
  geom_point(aes(color=temp), size=3) +
  #scale_color_gradient(low="blue",high="red") +
  scale_color_gradientn(colors=rev(rainbow(n=10,start=0,end=0.7))) +
  geom_text(aes(label=codNames),hjust=-0.1,vjust=-0.2,check_overlap = F) +
  #ylim(0,0.75) +
  #xlim(1.2,3) +
  ylab(expression(paste(lambda[2]/lambda[1]," ","(damping ratio)"))) +
  xlab("CV of spawning distribution")
  
dev.off()
  






# color palette
TEMPvec = c(8.4,7.13,7,7.4,4,11.9,5.8,6.5,9.57,1,8,8,1.75,3.50,2.5)
rbPal = colorRampPalette(c('red','blue'))
datacol = rbPal(10)[as.numeric(cut(TEMPvec,breaks=10))]

plot(x=cv,y=firsts,
     xlab="CV",
     ylab="lambda1", 
     main="CV vs leading eigen value",
     pch=16,
     xlim=c(1,3),
     ylim=c(0.5,1.5))
#legend("topright",legend=c(1:10),col =rbPal(10),pch=16)
text(x=cv,y=firsts,
     labels=codNames,pos=3)


#plot(x=levels(droplevels(eigen_cv$cv)),y=levels(droplevels(eigen_cv$firsts)),
#     xlab="CV",ylab="lambda 1", main="CV vs leading eigen value")
text(x=cv,y=seconds/firsts,
     labels=codNames,pos=3)


plot(x=cv,y=seconds/firsts,xlab="CV",
     ylab="lambda2/lambda1", 
     main="CV vs lambda2/lambda1",
     pch=16,
     xlim=c(1,3),
     ylim=c(0,1.5))
#text(x=cv,y=seconds/firsts,labels=codNames,pos=3)

# store Alists for different values of k
Alist = Alist0.5 
Alist0.2 


plot(x=cv,y=seconds/firsts,xlab="CV",
     ylab="lambda2/lambda1", 
     main="CV vs lambda2/lambda1",
     pch=16,
     xlim=c(1,3),
     ylim=c(0,1.5))