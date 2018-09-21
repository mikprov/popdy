# Plot 5 panel plot: lambda1 vs CV of spawning biomass distribution for different k values
# k = ranges 0-1, top row of Leslie (fecudities) is multiplied by k
# different k values signify harvested populations.

# Plan
# load libraries
# 1. create (new) Leslie matrices, multiply by k, calculate lambda1
# 2. calc LSB, plot spawning biomass at age for spot check, calculate CV of distributions
# 3. assemble df: codName, k, lambda1, CV, temp
# 4. plot 

library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(ggrepel)
library(devtools)
library(broom)
library(reshape2)

# 1. create Leslie matricies

# Load functions:
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/2_cod_functions.r")
# load cod data, break into separate populations
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/0_load_cod_data.r")


# ******* New Leslie Matrix ******** 

# prep output for new Leslie matricies
# Irish gives errors, remove from below
codNames <- c("Northsea","Coas","W_Baltic",
              "Faroe","NE_Arctic","Celtic",
              "Iceland","Kat","W_Scotland",
              "NGulf","GB","GM",
              "cod3NO","cod3M","cod2J3KL",
              "cod3Ps")

F.halfmax = 0 #for now, F is 0
LSBlist <- as.list(rep(NA,length(codNames)))
names(LSBlist) <- codNames
ks <- c(0.2,0.5,0.8,1) #kvalue=0 gives errors
tknot =0

# ---
# this loop calculates LSB at age at different values of k

for (i in 1:length(LSBlist)) { # for each population
  # load parms for cod pop i
  source(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_pops/',names(LSBlist)[i], '.r', sep=''))
  # this should load parms: L_inf, K, TEMP, maxage, B0, B1
  
  # for each pop calc LSB at age for different k values
  lsb.at.k = matrix(0,nrow=maxage,ncol=length(ks)) #create empty matrix to store LSB at age at k
  for (u in 1:length(ks)) { # step through k values to create LSB at age
    lsb.at.k[,u] = calculate_LSB_at_age_by_F(data=datalist[[i]], maxage=maxage,
                                             L_inf=L_inf, K=K, TEMP=TEMP, F.halfmax=F.halfmax,
                                             B0=B0,B1=B1)
    lsb.at.k[,u] = lsb.at.k[,u]*ks[u] #multiply by k
  }
  LSBlist[[i]] = lsb.at.k #LSBlist should have length of codNames
}

# ---
# using LSB values at age for different k values, generate 'leslie' matrix and get eigenvalues

e1 <- matrix(NA,nrow=length(codNames),ncol=length(ks)) #lambda1 for each pop at diff ks
e2 <- matrix(NA,nrow=length(codNames),ncol=length(ks)) #lambda2 for each pop at diff ks
e12 <- matrix(NA,nrow=length(codNames),ncol=length(ks))#dampratio for each pop at diff ks

for (t in 1:length(LSBlist)) { # for each pop
  # load parms for cod pop t
  source(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_pops/',names(LSBlist)[t], '.r', sep=''))
  
  for (y in 1:length(ks)) { #for each k value
    A = matrix(0,nrow=maxage,ncol=maxage) #create an empty 'leslie' matrix
    
    A[1,]=LSBlist[[t]][,y]/sum(LSBlist[[t]][,y]) #fill in relative LSB at age across top row
    A[1,]=A[1,]*ks[y]
    
    for(a in 2:maxage-1){ #filling in 'leslie' matrix
      A[a+1,a]=1} #fill in survival of 1 on subdiagonal
    
    # calculate eigenvalues of A and store them
    e1[t,y] <- extract_first_eigen_value(A) #e1matrix:col=ks, row=codNames
    e2[t,y] <- extract_second_eigen_value(A)
    e12[t,y] <- abs(e2[t,y])/e1[t,y]#damping ratio
    
  }
  
  
}
# adding labels and format into dfs
e1matrix <- as.data.frame(e1)  
colnames(e1matrix)<- ks
e1matrix <- cbind(codNames,e1matrix)
e1long <- melt(e1matrix, id.vars="codNames")
e1long$eigen <- rep('lambda1',length(e1long$codNames))

e2matrix <- as.data.frame(e2)  
colnames(e2matrix)<- ks
e2matrix <- cbind(codNames,e2matrix)
e2long <- melt(e2matrix, id.vars="codNames")
e2long$eigen <- rep('lambda2',length(e2long$codNames))

e12matrix <- as.data.frame(e12)  
colnames(e12matrix)<- ks
e12matrix <- cbind(codNames,e12matrix)
e12long <- melt(e12matrix, id.vars="codNames")
e12long$eigen <- rep('dampratio',length(e12long$codNames))

long <- rbind(e1long,e2long,e12long)
# ***************************************
# ***** read in eigentable *******
# ***************************************
eigentable <- read.csv(file="C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLSB/eigentable.csv",header=T)

long <- merge(long,eigentable[c("codNames","temp")],all.x=TRUE)
long <- merge(long,eigentable[c("codNames","mode_age")],all.x=TRUE)
long <- merge(long,eigentable[c("codNames","sd_mode")],all.x=TRUE)
long <- merge(long,eigentable[c("codNames","max_ages")],all.x=TRUE)
long <- merge(long,eigentable[c("codNames","cvs_mode")],all.x=TRUE)

# regressions
# TABLE 1: plot & regression: lambda1 vs CV, at k
m1 <- as.list(rep(NA,length(ks)))

for(n in 1:length(ks)){
  dd <- long[long$variable == ks[n] & long$eigen == "lambda1",]
  m1[[n]] <- lm(dd$value ~ dd$cvs_mode)
}
m1_tidy_lambda_vs_cv_k02 <- tidy(m1[[1]])
m1_tidy_lambda_vs_cv_k05 <- tidy(m1[[2]])
m1_tidy_lambda_vs_cv_k08 <- tidy(m1[[3]])
m1_tidy_lambda_vs_cv_k1 <- tidy(m1[[3]])
rm(dd,n,m1)
# combine into one table:
lambda_vs_cv_kall <- rbind(m1_tidy_lambda_vs_cv_k02,
                           m1_tidy_lambda_vs_cv_k05,
                           m1_tidy_lambda_vs_cv_k08,
                           m1_tidy_lambda_vs_cv_k1)
write.csv(lambda_vs_cv_kall,
          "C:/Users/provo/Documents/GitHub/popdy/cod_figures/regressions/m1_lambda1_vs_cv_kall.csv")


m2 <- as.list(rep(NA,length(ks)))

for(n in 1:length(ks)){
  dd <- long[long$variable == ks[n] & long$eigen == "lambda2",]
  m2[[n]] <- lm(dd$value ~ dd$cvs_mode)
}
m2_tidy_lambda2_vs_cv_k02 <- tidy(m2[[1]])
m2_tidy_lambda2_vs_cv_k05 <- tidy(m2[[2]])
m2_tidy_lambda2_vs_cv_k08 <- tidy(m2[[3]])
m2_tidy_lambda2_vs_cv_k1 <- tidy(m2[[3]])
rm(dd,n,m2)
# combine into one table:
lambda2_vs_cv_kall <- rbind(m2_tidy_lambda2_vs_cv_k02,
                           m2_tidy_lambda2_vs_cv_k05,
                           m2_tidy_lambda2_vs_cv_k08,
                           m2_tidy_lambda2_vs_cv_k1)
write.csv(lambda_vs_cv_kall,
          "C:/Users/provo/Documents/GitHub/popdy/cod_figures/regressions/m2_lambda2_vs_cv_kall.csv")




# plot
l1 <- ggplot(long[long$eigen=="lambda1",],aes(x=cvs_mode,y=value)) +
  geom_point(aes(col=mode_age)) + 
  geom_smooth(method="lm",se=FALSE,color="black") +
  facet_grid(. ~ variable) +
  scale_y_continuous(limits=c(0,1.1)) +
  geom_text_repel(data=long[long$eigen=="lambda1",],
                  aes(label = codNames,color=mode_age),
                  segment.color = "grey",
                  size = 2,
                  na.rm = TRUE) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(angle = 0)) +
  ylab(expression(paste(lambda[1]))) 



l2 <- ggplot(long[long$eigen=="lambda2",],aes(x=cvs_mode,y=value)) +
  geom_point(aes(col=mode_age)) + 
  geom_smooth(method="lm",se=FALSE,color="black") +
  facet_grid(. ~ variable) +
  scale_y_continuous(limits=c(0,1.1)) +
  geom_text_repel(data=long[long$eigen=="lambda2",],
                  aes(label = codNames,color=mode_age),
                  segment.color = "grey",
                  size = 2,
                  na.rm = TRUE) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(angle = 0)) +
  ylab(expression(paste(lambda[2]))) 


d <- ggplot(long[long$eigen=="dampratio",],aes(x=cvs_mode,y=value)) +
  geom_point(aes(col=mode_age)) + 
  geom_smooth(method="lm",se=FALSE,color="black") +
  facet_grid(. ~ variable) +
  scale_y_continuous(limits=c(0,1.1)) +
  geom_text_repel(data=long[long$eigen=="dampratio",],
                  aes(label = codNames,color=mode_age),
                  segment.color = "grey",
                  size = 2,
                  na.rm = TRUE) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(angle = 0)) +
  ylab(expression(abs(paste(lambda[2]))/paste(lambda[1]))) 


dage <- ggplot(long[long$eigen=="dampratio",],aes(x=mode_age,y=value)) +
  geom_point(aes(col=mode_age)) + 
  geom_smooth(method="lm",se=FALSE,color="black") +
  facet_grid(. ~ variable) +
  scale_y_continuous(limits=c(0,1.1)) +
  geom_text_repel(data=long[long$eigen=="dampratio",],
                  aes(label = codNames,color=mode_age),
                  segment.color = "grey",
                  size = 2,
                  na.rm = TRUE) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(angle = 0)) +
  ylab(expression(abs(paste(lambda[2]))/paste(lambda[1]))) 


pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/k_panel_plots_updated.pdf', width=10, height=7)
l1
l2
d
dage
dev.off()










# ***************************************
# ***** traditional Leslie matrix *******
# ***************************************
# prep for looping
F.halfmax = 0 #for now, F is 0
Alist = as.list(rep(NA,length(datalist))) # store Leslie matrix
names(Alist) = codNames
k <- c(0,0.2,0.5,0.8,1)
dfexp = as.list(rep(NA,length(k)))
names(dfexp) = k
eigenvals1 = rep(NA,length(codNames))
eigenvals2 = rep(NA,length(codNames))
eigenvals1.2 = rep(NA,length(codNames))
tknot=0

for (j in 1:length(k)) { #step through k values
  
  for (i in 1:length(names(Alist))) { # step through each dataset in datalist
    # load parms for cod pop i
    source(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_pops/',names(Alist)[i], '.r', sep=''))
    # this should load parms: L_inf, K, TEMP, maxage
    
    out=assemble_Leslie(data=datalist[[i]], maxage=maxage, K=K, L_inf=L_inf, TEMP=TEMP,
                        F.halfmax=0,tknot=0)
    Alist[[i]]=out$A
    Alist[[i]][1,] <- Alist[[i]][1,]*k[j] #multiply by slope (k), Jacobian matrix
    eigenvals1[i] <- extract_first_eigen_value(Alist[[i]])
    eigenvals2[i] <- extract_second_eigen_value(Alist[[i]])
    eigenvals1.2[i] <- eigenvals2[i] / eigenvals1[i]
    # remove pop parms for next loop 
    rm(K,L_inf,maxage)
  }
  
  kval <- rep(k[j],length(codNames)) #store lambda values from kval
  dfexp[[j]] <- cbind(kval,codNames,eigenvals1,eigenvals2,eigenvals1.2)
  
}
rm(out,A50,B0,B1,H,i,j,kval,MG,Mp,name,S50,TEMP,theta0,theta1) #clean up
df1 <- do.call(rbind.data.frame,dfexp) #
rownames(df1) <- c()
head(df1)
df1[df1$codNames == "Northsea",]



# 2. calc LSB, plot spawning biomass at age, calculate CV of distributions
lsblist = as.list(rep(NA,length(k)))
names(lsblist) = k

for (i in 1:length(k)) { #walk through k values
  klist = as.list(rep(NA,length(Alist))) #create an empty list for small dfs
  names(klist) = codNames #naming it with codNames
  
  for (j in 1:length(codNames)) {
    datax = datalist[[j]] #make sure cod data is loaded
    # load parms for cod pop j
    source(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_pops/',names(datalist)[j], '.r', sep=''))
    # this should load parms: L_inf, K, TEMP, maxage
    lep <- calculate_LSB_at_age_by_F(data=datax,littlek=k[i],maxage=maxage,L_inf=L_inf,K=K,TEMP=TEMP,F.halfmax=0)
    lep <- cbind(lep,rep(k[i],length(lep)),rep(codNames[j]))
    colnames(lep) <- c("LEP","k","codNames")
    klist[[j]] <- lep
                        
  }
  kdf <- do.call(rbind.data.frame,klist)
  rownames(kdf) <- c()
  lsblist[[i]] <- kdf
}
rm(kdf,lep,klist,i,j,A50,B0,B1,K,L_inf,maxage,MG,Mp,name,S50,TEMP,theta0,theta1) #clean up
LEPdff <- do.call(rbind.data.frame,lsblist)
rownames(LEPdff) <- c()
head(LEPdff)
LEPdff$LEP <- as.numeric(levels(LEPdff$LEP))[LEPdff$LEP] #silly factors
LEPdff[LEPdff$codNames == "Irish" ,] #just looking

# 3. this set of loops calc spawning biomass distribution for different k values
# loop over k values wasn't working! doing it the long way!
cvinfo <- as.list(rep(NA,length(k)))
names(cvinfo) <- k

# 1st value of k
y = 1 #y is the index in the vector of k values
cvs_mode = rep(NA, length(codNames)) #empty vectors 
mode_age = rep(NA, length(codNames))
sd_mode = rep(NA, length(codNames))

for (i in 1:length(codNames)) { # step through each cod population
  d <- LEPdff[LEPdff$k == k[y] & LEPdff$codNames == codNames[i],]$LEP
  p_spawn = d / sum(d) # probability of spawning at age = LSB at age/total LSB
  Age = seq(from=1,to=length(d),by=1)
  p_table = data.frame(cbind(Age,p_spawn))
  #keep <- p_table[which(p_table$p_spawn > 0.01),] # remove probabilities less than 0.01
  keep <- p_table #use this if I don't specify min limit of prob of spawning
  # using mode, not mean
  mode_age[i] = keep$Age[which.max(keep$p_spawn)] # what is the age with highest probability?
  sd_mode[i] = sqrt( sum(keep$p_spawn*(keep$Age-mode_age[i])^2) ) # stdev
  cvs_mode[i] = sd_mode[i]/mode_age[i] # coefficient of variation 
  
}
kval <- rep(k[y],length(codNames))
out <- cbind(codNames,kval,mode_age,sd_mode,cvs_mode)
cvinfo[[y]] <- out #store for later
rm(i,y,kval,out,keep,d,p_spawn,Age,mode_age,sd_mode,cvs_mode) #clean up

# 2nd value of k
y = 2
cvs_mode = rep(NA, length(codNames)) #empty vectors 
mode_age = rep(NA, length(codNames))
sd_mode = rep(NA, length(codNames))

for (i in 1:length(codNames)) { # step through each cod population
  d <- LEPdff[LEPdff$k == k[y] & LEPdff$codNames == codNames[i],]$LEP
  p_spawn = d / sum(d) # probability of spawning at age = LSB at age/total LSB
  Age = seq(from=1,to=length(d),by=1)
  p_table = data.frame(cbind(Age,p_spawn))
  #keep <- p_table[which(p_table$p_spawn > 0.01),] # remove probabilities less than 0.01
  keep <- p_table #use this if I don't specify min limit of prob of spawning
  # using mode
  mode_age[i] = keep$Age[which.max(keep$p_spawn)] # what is the age with highest probability?
  sd_mode[i] = sqrt( sum(keep$p_spawn*(keep$Age-mode_age[i])^2) ) # stdev
  cvs_mode[i] = sd_mode[i]/mode_age[i] # coefficient of variation 
  
}
kval <- rep(k[y],length(codNames))
out <- cbind(codNames,kval,mode_age,sd_mode,cvs_mode)
cvinfo[[y]] <- out
rm(i,y,kval,out,keep,d,p_spawn,Age,mode_age,sd_mode,cvs_mode)

# 3rd value of k
y = 3
cvs_mode = rep(NA, length(codNames)) #empty vectors 
mode_age = rep(NA, length(codNames))
sd_mode = rep(NA, length(codNames))

for (i in 1:length(codNames)) { # step through each cod population
  d <- LEPdff[LEPdff$k == k[y] & LEPdff$codNames == codNames[i],]$LEP
  p_spawn = d / sum(d) # probability of spawning at age = LSB at age/total LSB
  Age = seq(from=1,to=length(d),by=1)
  p_table = data.frame(cbind(Age,p_spawn))
  #keep <- p_table[which(p_table$p_spawn > 0.01),] # remove probabilities less than 0.01
  keep <- p_table #use this if I don't specify min limit of prob of spawning
  # using mode
  mode_age[i] = keep$Age[which.max(keep$p_spawn)] # what is the age with highest probability?
  sd_mode[i] = sqrt( sum(keep$p_spawn*(keep$Age-mode_age[i])^2) ) # stdev
  cvs_mode[i] = sd_mode[i]/mode_age[i] # coefficient of variation 
  
}
kval <- rep(k[y],length(codNames))
out <- cbind(codNames,kval,mode_age,sd_mode,cvs_mode)
cvinfo[[y]] <- out
rm(i,y,kval,out,keep,d,p_spawn,Age,mode_age,sd_mode,cvs_mode)

# 4th value of k
y = 4
cvs_mode = rep(NA, length(codNames)) #empty vectors 
mode_age = rep(NA, length(codNames))
sd_mode = rep(NA, length(codNames))

for (i in 1:length(codNames)) { # step through each cod population
  d <- LEPdff[LEPdff$k == k[y] & LEPdff$codNames == codNames[i],]$LEP
  p_spawn = d / sum(d) # probability of spawning at age = LSB at age/total LSB
  Age = seq(from=1,to=length(d),by=1)
  p_table = data.frame(cbind(Age,p_spawn))
  #keep <- p_table[which(p_table$p_spawn > 0.01),] # remove probabilities less than 0.01
  keep <- p_table #use this if I don't specify min limit of prob of spawning
  # using mode
  mode_age[i] = keep$Age[which.max(keep$p_spawn)] # what is the age with highest probability?
  sd_mode[i] = sqrt( sum(keep$p_spawn*(keep$Age-mode_age[i])^2) ) # stdev
  cvs_mode[i] = sd_mode[i]/mode_age[i] # coefficient of variation 
  
}
kval <- rep(k[y],length(codNames))
out <- cbind(codNames,kval,mode_age,sd_mode,cvs_mode)
cvinfo[[y]] <- out
rm(i,y,kval,out,keep,d,p_spawn,Age,mode_age,sd_mode,cvs_mode)

# 5th value of k
y = 5
cvs_mode = rep(NA, length(codNames)) #empty vectors 
mode_age = rep(NA, length(codNames))
sd_mode = rep(NA, length(codNames))

for (i in 1:length(codNames)) { # step through each cod population
  d <- LEPdff[LEPdff$k == k[y] & LEPdff$codNames == codNames[i],]$LEP
  p_spawn = d / sum(d) # probability of spawning at age = LSB at age/total LSB
  Age = seq(from=1,to=length(d),by=1)
  p_table = data.frame(cbind(Age,p_spawn))
  #keep <- p_table[which(p_table$p_spawn > 0.01),] # remove probabilities less than 0.01
  keep <- p_table #use this if I don't specify min limit of prob of spawning
  # using mode
  mode_age[i] = keep$Age[which.max(keep$p_spawn)] # what is the age with highest probability?
  sd_mode[i] = sqrt( sum(keep$p_spawn*(keep$Age-mode_age[i])^2) ) # stdev
  cvs_mode[i] = sd_mode[i]/mode_age[i] # coefficient of variation 
  
}
kval <- rep(k[y],length(codNames))
out <- cbind(codNames,kval,mode_age,sd_mode,cvs_mode)
cvinfo[[y]] <- out
rm(i,y,kval,out,keep,d,p_spawn,Age,mode_age,sd_mode,cvs_mode)

# collapse list of dfs with cvs into a single df
cvinfodf <- do.call(rbind.data.frame,cvinfo)
rownames(cvinfodf) <- c()
head(cvinfodf)
cvinfodf[30:50,] #just looking

# 4. assemble df for plotting: codNames, k, lambda1, cv, mode, sd, temp
df3 <- left_join(cvinfodf,df1)
head(df3)
df3[30:55,] #just looking

# read in temp from eigentable, i should really clean this up!
eigentable = read.csv("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLSB/eigentable.csv",
                      header=TRUE,stringsAsFactors = FALSE)
eigentable = as.data.frame(eigentable)
temptable <- subset(eigentable, select=c("codNames","temp"))
df4 <- left_join(df3,temptable)
head(df4)
df4[df4$codNames %in% c("Coas","Northsea"),]
df4$mode_age <- as.numeric(levels(df4$mode_age))[df4$mode_age] #silly factors
df4$sd_mode <- as.numeric(levels(df4$sd_mode))[df4$sd_mode] #silly factors
df4$cvs_mode <- as.numeric(levels(df4$cvs_mode))[df4$cvs_mode] #silly factors
df4$eigenvals1 <- as.numeric(levels(df4$eigenvals1))[df4$eigenvals1] #silly factors
df4$eigenvals2 <- as.numeric(levels(df4$eigenvals2))[df4$eigenvals2] #silly factors
df4$eigenvals1.2 <- as.numeric(levels(df4$eigenvals1.2))[df4$eigenvals1.2] #silly factors
df4[30:55,]
df4[is.na(df4)] <- 0
rm(temptable,df3,cvinfodf,cvinfo) #clean up

# 5. plot
# our df for plotting is df4
p_temp <- ggplot(df4,aes(x=cvs_mode,y=eigenvals1.2)) +
  geom_point(aes(col=temp)) + 
  facet_grid(. ~ kval) +
  scale_color_gradientn(colors=rev(rainbow(n=16,start=0,end=0.7))) +
  geom_text_repel(data=df4,
                aes(label = codNames,color=temp),
                segment.color = "grey",
                size = 2,
                na.rm = TRUE) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ylab(expression(paste("|",lambda[2],"|","/","|",lambda[1],"|"))) 


p_mode <- ggplot(df4,aes(x=cvs_mode,y=eigenvals1.2)) +
  geom_point(aes(col=mode_age)) + 
  geom_smooth(method="lm",se=FALSE) +
  facet_grid(. ~ kval) +
  geom_text_repel(data=df4,
                  aes(label = codNames,color=mode_age),
                  segment.color = "grey",
                  size = 2,
                  na.rm = TRUE) +
  theme_bw() + 
  theme(axis.title.y = element_text(angle = 0),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  ylab(expression(frac(paste("|",lambda[2],"|"),paste("|",lambda[1],"|"))))

pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/k_panel_plots_lambda2over1.pdf', width=10, height=7)
p_temp
p_mode
dev.off()


# --------------
# regressions
# --------------

# ---
# TABLE 1: plot & regression: lambda1 vs CV, at k
m1 <- as.list(rep(NA,length(ks)))
for(n in 1:length(k)){
  dd <- df4[df4$kval==ks[n],]
  m1[[n]] <- lm(dd$eigenvals1 ~ dd$cvs_mode)
}
m1_tidy_lambda_vs_cv_k0 <- tidy(m1[[1]])
m1_tidy_lambda_vs_cv_k02 <- tidy(m1[[2]])
m1_tidy_lambda_vs_cv_k05 <- tidy(m1[[3]])
m1_tidy_lambda_vs_cv_k08 <- tidy(m1[[4]])
rm(dd,n,m1)
# combine into one table:
lambda_vs_cv_kall <- rbind(m1_tidy_lambda_vs_cv_k0,
                           m1_tidy_lambda_vs_cv_k02,
                           m1_tidy_lambda_vs_cv_k05,
                           m1_tidy_lambda_vs_cv_k08,
                           m1_tidy_lambda_vs_cv_k1)
write.csv(lambda_vs_cv_kall,
          "C:/Users/provo/Documents/GitHub/popdy/cod_figures/regressions/m1_lambda1_vs_cv_kall.csv")

#confint(m1[[2]])
# PLOT 1: plot lambda1 vs cv, at each k
p_lambda_vs_cv <- ggplot(df4,aes(x=cvs_mode,y=eigenvals1)) +
  geom_point(aes(col=mode_age)) +
  #scale_color_gradientn(colors=rev(rainbow(n=17,start=0,end=0.7))) +
  geom_smooth(method="lm",se=FALSE,color="black") +
  facet_grid(. ~ kval) +
  geom_text_repel(data=df4,
                  aes(label = codNames,color=mode_age),
                  segment.color = "grey",
                  size = 2,
                  na.rm = TRUE) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(angle = 0)) +
  ylab(expression(paste(lambda[1]))) 
pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/k_panel_plots_labmda1_vs_cv.pdf', width=10, height=5)
p_lambda_vs_cv
dev.off()
# ---
# TABLE 2 & 3: mean lambda1 vs k
# does stdev of lambda1 narrow at larger k values?
# does mean of lambda1 change with k values?
st <- rep(NA,length(k)-1) #subtract 1 bc I'm removing k=0 from regression
me <- rep(NA,length(k)-1) #subtract 1 bc I'm removing k=0 from regression
for (n in 1:length(k)-1) {
  dd <- df4[df4$kval==k[n+1],]
  st[n] <- sd(dd$eigenvals1)
  me[n] <- mean(dd$eigenvals1)} 
k.nok0 <- k[2:5]
ddst <- as.data.frame(cbind(k.nok0,st))
ddme <- as.data.frame(cbind(k.nok0,me))
colnames(ddst) <- c("k","st")
colnames(ddme) <- c("k","me")
m2st <- lm(ddst$st ~ ddst$k)
m2me <- lm(ddme$me ~ ddme$k)
m2_tidy_stlambda_vs_k <- tidy(m2st)
m2_tidy_melambda_vs_k <- tidy(m2me)

# PLOT 2 & 3: plot stdev vs k and mean vs k
ggplot(ddme,aes(x=k,y=me)) +
  geom_point() + ylab(expression(paste("mean ",lambda[1]))) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

ggplot(ddme,aes(x=k,y=st)) +
  geom_point() + ylab(expression(paste("stdev ",lambda[1]))) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))
rm(ddst,ddme,n,m2me,m2st,st,me,dd,k.nok0)

# ---
# TABLE 4: regressions for damping ratio vs cv at each k
m3dr <- as.list(rep(NA,length(k)))
for(n in 1:length(k)){
  dd <- df4[df4$kval==k[n],]
  m3dr[[n]] <- lm(dd$eigenvals1.2~dd$cvs_mode)}
m3_tidy_dr_vs_cv_k0 <- tidy(m3dr[[1]])
m3_tidy_dr_vs_cv_k02 <- tidy(m3dr[[2]])
m3_tidy_dr_vs_cv_k05 <- tidy(m3dr[[3]])
m3_tidy_dr_vs_cv_k08 <- tidy(m3dr[[4]])
m3_tidy_dr_vs_cv_k1 <- tidy(m3dr[[5]])
rm(dd,n,m3dr)

# PLOT 4:
p_mode <- ggplot(df4,aes(x=cvs_mode,y=eigenvals1.2)) +
  geom_point(aes(col=mode_age)) + 
  geom_smooth(method="lm",se=FALSE,color="black") +
  facet_grid(. ~ kval) +
  geom_text_repel(data=df4,
                  aes(label = codNames,color=mode_age),
                  segment.color = "grey",
                  size = 2,
                  na.rm = TRUE) +
  theme_bw() + 
  theme(axis.title.y = element_text(angle = 0),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  ylab(expression(frac(paste("|",lambda[2],"|"),paste(lambda[1]))))

pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/k_panel_plots_labmda2over1_vs_mode.pdf', 
    width=10, height=5)
p_mode
dev.off()
# PLOT 4 check: mean 1/damping ratio (mean lambda2/lambda1) vs k
st <- rep(NA,length(k)-1) #subtract 1 bc I'm removing k=0 from regression
me <- rep(NA,length(k)-1) #subtract 1 bc I'm removing k=0 from regression
for (n in 1:length(k)-1) {
  dd <- df4[df4$kval==k[n+1],]
  st[n] <- sd(dd$eigenvals1.2)
  me[n] <- mean(dd$eigenvals1.2)} 
k.nok0 <- k[2:5]
ddst <- as.data.frame(cbind(k.nok0,st))
ddme <- as.data.frame(cbind(k.nok0,me))
colnames(ddst) <- c("k","st")
colnames(ddme) <- c("k","me")
m3_dr_vs_k <- tidy(lm(ddst$st~ddst$k))
rho_vs_k <- ggplot(ddst,aes(x=k,y=me)) +
  geom_point()+
  ylab(expression("mean ",frac(paste("|",lambda[2],"|"),paste(lambda[1])))) 
pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/k_panel_plots_meanlabmda2over1_vs_k.pdf', 
    width=10, height=5)
rho_vs_k
dev.off()
# ---
# TABLE 5: cv vs. mode, cv vs. stdev
# subset data,df4
df4.sub <- df4[df4$kval==1,]
df4.sub[] <- lapply(df4.sub, function(x) if(is.factor(x)) factor(x) else x)

cv.stdev <- tidy(lm(df4.sub$cvs_mode ~ df4.sub$sd_mode))
cv.mode <- tidy(lm(df4.sub$cvs_mode ~ df4.sub$mode_age))

# PLOT 5: CV vs stdev
cv_vs_mode <- ggplot(df4.sub,aes(x=mode_age,y=cvs_mode)) +
  geom_point(aes(col=mode_age)) + 
  #scale_color_gradientn(colors=rev(rainbow(n=16,start=0,end=0.7))) +
  geom_text_repel(data=df4.sub,
                  aes(label = codNames,color=mode_age),
                  segment.color = "grey",
                  size = 4,
                  na.rm = TRUE) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("CV about the mode") + xlab("peak spawning age")

pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/cv_vs_mode.pdf', 
    width=7, height=6)
cv_vs_mode
dev.off()

cv_vs_stdev <- ggplot(df4.sub,aes(x=sd_mode,y=cvs_mode)) +
  geom_point(aes(col=mode_age)) + 
  #scale_color_gradientn(colors=rev(rainbow(n=16,start=0,end=0.7))) +
  geom_text_repel(data=df4.sub,
                  aes(label = codNames,color=mode_age),
                  segment.color = "grey",
                  size = 4,
                  na.rm = TRUE) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("CV about the mode") + xlab("stdev of spawning biomass distribution")
pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/cv_vs_stdev.pdf', 
    width=7, height=6)
cv_vs_stdev
dev.off()


# export tables
write.csv(m1_tidy_lambda_vs_cv_k0,
          "C:/Users/provo/Documents/GitHub/popdy/cod_figures/regressions/m1_lambda1vsmode_k0.csv")
write.csv(m1_tidy_lambda_vs_cv_k02,
          "C:/Users/provo/Documents/GitHub/popdy/cod_figures/regressions/m1_lambda1vsmode_k02.csv")
write.csv(m1_tidy_lambda_vs_cv_k05,
          "C:/Users/provo/Documents/GitHub/popdy/cod_figures/regressions/m1_lambda1vsmode_k05.csv")
write.csv(m1_tidy_lambda_vs_cv_k08,
          "C:/Users/provo/Documents/GitHub/popdy/cod_figures/regressions/m1_lambda1vsmode_k08.csv")
write.csv(m1_tidy_lambda_vs_cv_k1,
          "C:/Users/provo/Documents/GitHub/popdy/cod_figures/regressions/m1_lambda1vsmode_k1.csv")
# combine into one table:
lambda_vs_cv_kall <- rbind(m1_tidy_lambda_vs_cv_k0,
                             m1_tidy_lambda_vs_cv_k02,
                             m1_tidy_lambda_vs_cv_k05,
                             m1_tidy_lambda_vs_cv_k08,
                             m1_tidy_lambda_vs_cv_k1)
write.csv(lambda_vs_cv_kall,
          "C:/Users/provo/Documents/GitHub/popdy/cod_figures/regressions/m1_lambda1_vs_cv_kall.csv")

write.csv(m2_tidy_melambda_vs_k,
          "C:/Users/provo/Documents/GitHub/popdy/cod_figures/regressions/m2_melambda1_vs_k.csv")
write.csv(m2_tidy_stlambda_vs_k,
          "C:/Users/provo/Documents/GitHub/popdy/cod_figures/regressions/m2_stlambda1_vs_k.csv")

write.csv(m3_tidy_dr_vs_cv_k0,
          "C:/Users/provo/Documents/GitHub/popdy/cod_figures/regressions/m3_dampratio_vs_cv_k0.csv")
write.csv(m3_tidy_dr_vs_cv_k02,
          "C:/Users/provo/Documents/GitHub/popdy/cod_figures/regressions/m3_dampratio_vs_cv_k02.csv")
write.csv(m3_tidy_dr_vs_cv_k05,
          "C:/Users/provo/Documents/GitHub/popdy/cod_figures/regressions/m3_dampratio_vs_cv_k05.csv")
write.csv(m3_tidy_dr_vs_cv_k08,
          "C:/Users/provo/Documents/GitHub/popdy/cod_figures/regressions/m3_dampratio_vs_cv_k08.csv")
write.csv(m3_tidy_dr_vs_cv_k1,
          "C:/Users/provo/Documents/GitHub/popdy/cod_figures/regressions/m3_dampratio_vs_cv_k1.csv")

m3_tidy_dr_vs_cv_kall <- rbind(m3_tidy_dr_vs_cv_k0,
                               m3_tidy_dr_vs_cv_k02,
                               m3_tidy_dr_vs_cv_k05,
                               m3_tidy_dr_vs_cv_k08,
                               m3_tidy_dr_vs_cv_k1)
write.csv(m3_tidy_dr_vs_cv_kall,
          "C:/Users/provo/Documents/GitHub/popdy/cod_figures/regressions/m3_dampratio_vs_cv_kall.csv")

write.csv(m3_dr_vs_k,
          "C:/Users/provo/Documents/GitHub/popdy/cod_figures/regressions/m3_meandr_vs_k.csv")

write.csv(cv.stdev,
          "C:/Users/provo/Documents/GitHub/popdy/cod_figures/regressions/m4_tidy_cv_vs_sd.csv")
write.csv(cv.mode,
          "C:/Users/provo/Documents/GitHub/popdy/cod_figures/regressions/m4_tidy_cv_vs_mode.csv")
