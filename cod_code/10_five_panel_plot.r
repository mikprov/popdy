# Plot 5 panel plot: lambda1 vs CV of spawning biomass distribution for different k values
# k = ranges 0-1, top row of Leslie (fecudities) is multiplied by k
# different k values signify harvested populations.

# Plan
# load libraries
# 1. create Leslie matrices, multiply by k, calculate lambda1
# 2. calc LSB, plot spawning biomass at age for spot check, calculate CV of distributions
# 3. assemble df: codName, k, lambda1, CV, temp
# 4. plot 

library(ggplot2)
library(gridExtra)
library(dplyr)
library(ggrepel)

# 1. create Leslie matricies

# Load functions:
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/2_cod_functions.r")
# load cod data, break into separate populations
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/0_load_cod_data.r")

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

# for loop: for each population, create a Leslie, multiply by k, calculate & store lambda1

for (j in 1:length(k)) { #step through k values

  for (i in 1:length(names(Alist))) { # step through each dataset in datalist
    # load parms for cod pop i
    source(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_pops/',names(Alist)[i], '.r', sep=''))
    # this should load parms: L_inf, K, TEMP, maxage
    out=assemble_Leslie(data=datalist[[i]], maxage=maxage, K=K, L_inf=L_inf, TEMP=TEMP,
                             F.halfmax=0,tknot=0)
    Alist[[i]]=out$A
    Alist[[i]][1,] <- Alist[[i]][1,]*k[j]
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

# this set of loops calc spawning biomass distribution for different k values
# loop over k values wasn't working! doing it the long way!
cvinfo <- as.list(rep(NA,length(k)))
names(cvinfo) <- k

# 1st value of k
y = 1
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

# 3. assemble df for plotting: codNames, k, lambda1, cv, mode, sd, temp
df3 <- left_join(cvinfodf,df1)
head(df3)
df3[30:55,] #just looking

# read in temp from eigentable, i should really clean this up!
eigentable = read.csv("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLSB/eigentable5.csv",
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
p_temp <- ggplot(df4,aes(x=cvs_mode,y=eigenvals)) +
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
  ylab(expression(paste(lambda[2],"/","|",lambda[1],"|"))) 


p_mode <- ggplot(df4,aes(x=cvs_mode,y=eigenvals)) +
  geom_point(aes(col=mode_age)) + 
  facet_grid(. ~ kval) +
  geom_text_repel(data=df4,
                  aes(label = codNames,color=mode_age),
                  segment.color = "grey",
                  size = 2,
                  na.rm = TRUE) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ylab(expression(paste("|",lambda[2],"|","/",lambda[1])))

pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/k_panel_plots_lambda2over1.pdf', width=10, height=7)
p_temp
p_mode
dev.off()


# plot lambda1 vs mode, at each k
ggplot(df4,aes(x=mode_age,y=eigenvals1)) +
  geom_point(aes(col=temp)) +
  scale_color_gradientn(colors=rev(rainbow(n=17,start=0,end=0.7))) +
  facet_grid(. ~ kval) +
  geom_text_repel(data=df4,
                  aes(label = codNames),
                  segment.color = "grey",
                  size = 2,
                  na.rm = TRUE) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  ylab(expression(paste(lambda[1]))) 
# regression: lambda1 vs mode, at k
m1 <- as.list(rep(NA,length(k)))
for(n in 1:length(k)){
  dd <- df4[df4$kval==k[n],]
  m1[[n]] <- lm(dd$eigenvals1 ~ dd$mode_age)
}
rm(dd,n)
summary(m1[[2]])
#confint(m1)

# does stdev of lambda1 narrow at larger k values?
st <- rep(NA,length(k))
for (n in 1:length(k)) {
  dd <- df4[df4$kval==k[n],]
  st[n] <- sd(dd$eigenvals1)} 
rm(dd,n)
dd <- as.data.frame(cbind(k,st))
m2st <- lm(dd$st~dd$k)
summary(m2st)
# does mean of lambda1 change with k values?
me <- rep(NA,length(k))
for (n in 1:length(k)) {
  dd <- df4[df4$kval==k[n],]
  me[n] <- mean(dd$eigenvals1)} 
dd <- as.data.frame(cbind(k,me))
m2me <- lm(dd$me~dd$k)
summary(m2me)
rm(dd,n)

# regressions for damping ratio vs cv
m3dr <- as.list(rep(NA,length(k)))
for(n in 1:length(k)){
  dd <- df4[df4$kval==k[n],]
  m3dr[[n]] <- lm(dd$eigenvals1.2~dd$cvs_mode)}
rm(dd,n)

# plot & regression for cv vs. mode, cv vs. stdev
dfsm <- df4[df4$kval==0.2,]
ggplot(dfsm,aes(x=mode_age,y=cvs_mode)) +
  geom_point() 
ggplot(dfsm,aes(x=sd_mode,y=cvs_mode)) +
  geom_point() 
m4cv.mode <- lm(dfsm$cvs_mode ~ dfsm$mode_age)
summary(m4cv.mode)
m4cv.sd <- lm(dfsm$cvs_mode ~ dfsm$sd_mode)
summary(m4cv.sd)
