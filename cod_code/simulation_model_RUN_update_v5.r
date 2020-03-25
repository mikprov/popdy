# Run the simulation model
# by: mikaela provost
# last edited: Feb 19, 2019
# ===================================================================

library(tidyr)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(directlabels)
library(grid)
library(lattice)
library(tidyverse)
#install.packages("stargazer")
library(stargazer)
library(RColorBrewer)
library(ggrepel)


# ---
# load the simulation model
source("C:/Users/Mikaela/Documents/GitHub/popdy/cod_code/simulation_model_cod_v3.r")

# load functions --> ***** CHECK for MG or MP *****
source("C:/Users/Mikaela/Documents/GitHub/popdy/cod_code/2_cod_functions.r")

# load peak spawning age info
eigentable = read.csv("C:/Users/Mikaela/Documents/GitHub/popdy/cod_code/mikaelaLSB/eigentable_MM.csv",
                      header=TRUE,stringsAsFactors = FALSE)
eigentable = as.data.frame(eigentable)

# *******************
# depending on which M parm, run the right code
# *******************
# MG
eigentable$codNames <- factor(eigentable$codNames,levels=eigentable[order(eigentable$mode_ageMG,eigentable$temp),]$codNames)
eigentable$codNames_plot <- factor(eigentable$codNames_plot,levels=eigentable[order(eigentable$mode_ageMG,eigentable$temp),]$codNames_plot)
eigentable <- eigentable %>% select(-mode_ageMP,-sd_modeMP,-cvs_modeMP) %>% rename(mode_age=mode_ageMG,sd_mode=sd_modeMG,cvs_mode=cvs_modeMG)
# MP
#eigentable$codNames <- factor(eigentable$codNames,levels=eigentable[order(eigentable$mode_ageMP,eigentable$temp),]$codNames)
#eigentable$codNames_plot <- factor(eigentable$codNames_plot,levels=eigentable[order(eigentable$mode_ageMP,eigentable$temp),]$codNames_plot)
#eigentable <- eigentable %>% select(-mode_ageMG,-sd_modeMG,-cvs_modeMG) %>% rename(mode_age=mode_ageMP,sd_mode=sd_modeMP,cvs_mode=cvs_modeMP)
# *******************


# make all LEPs equal
conLEP = 1.1
# adjust fecundities by this much
adjFec = round(1/conLEP,digits=1)

#selectedkvals <- round(1/selectedalphas,digits=2)
selectedkvals <- c(0.15,0.4,0.6,0.85)
selectedalphas <- round(1/(selectedkvals*conLEP^2),digits=2)
#kvals = round(seq(from=0.1,to=1,by=0.01),digits=2)

alphas <- c(8.26,7.51,6.36,5.51,4.59,3.31,2.07,1.38,0.97)
kvals = round(1/(alphas*conLEP^2),digits=2)

# setting up for later
codNames <- eigentable$codNames
codNames_ordered_by_peak <- levels(eigentable$codNames)
codNames_ordered_by_peak_plot <- levels(eigentable$codNames_plot)

# *************************************** #
# (Fig S1: spawning biomass distributions)
# see script -- 3_plot_spawning_biomass_distributions.r
# *************************************** #

# *************************************** #
# (Fig 3:eigenvalue) Generate Leslie matricies for diff F values (create Leslie arrays)
# *************************************** #
Aarray = as.list(rep(NA,length(codNames))) #Leslie matrix storage for each k value for pop i
names(Aarray) <- codNames
LeslieoutNEARL <- as.list(rep(NA,length=length(codNames)))
names(LeslieoutNEARL) <- codNames
# store eigenvalues here
eigenvals1 = matrix(NA,nrow=length(kvals),ncol=length(codNames)) 
eigenvals2 = matrix(NA,nrow=length(kvals),ncol=length(codNames))
eigenvals12 = matrix(NA,nrow=length(kvals),ncol=length(codNames))
# store max fecundities
maxfecunds = matrix(NA,nrow=length(kvals),ncol=length(codNames))
# store LEP values
LEPvals = matrix(NA,nrow=length(kvals),ncol=length(codNames))
# store top row sums
toprowsums = matrix(NA,nrow=length(kvals),ncol=length(codNames))
# store new LEP with adjusted fecundities (divide by half of original LEP)
newLEPs = matrix(NA,nrow=length(kvals),ncol=length(codNames))


for (i in 1:length(codNames)){ #for each pop i
  # load parms for cod pop i: L_inf, K (for vonB), TEMP, maxage,B0,B1 (matur)
  source(file = paste('C:/Users/Mikaela/Documents/GitHub/popdy/cod_pops/',codNames[i], '.r', sep=''))
  
  Lesliearray <- array(NA,c(maxage,maxage,length(kvals))) #store Leslie matricies for simulations
  #Leslie_eigens <- array(NA,c(maxage,maxage,length(kvals))) #store Leslie matricies for eigen analysis
  e1 = rep(NA,length=length(kvals)) #store lambda1
  e2 = rep(NA,length=length(kvals)) #store lambda2
  e12 = rep(NA,length=length(kvals)) #store inverse damping ratio
  maxfec = rep(NA,length=length(kvals)) #store max fec values
  LEPs = rep(NA,length=length(kvals))
  toprowsum = rep(NA,length=length(kvals))
  newLEP = rep(NA,length=length(kvals)) #recalculte LEP with adjusted fs
  
  for (k in 1:length(kvals)){ #step through k values 
    # create Leslie matrix:
    Leslieout = assemble_Leslie(maxage=maxage, K=K, L_inf=L_inf, 
                                F.halfmax=0, B0=B0, B1=B1, tknot=0)
    Jacobian.sim <- Leslieout$A
    Jacobian.eig <- Leslieout$A
    
    # max fecundity values
    maxfec[k] <- max(Leslieout$A[1,])
    
    # LEP for this popualtion with the two mortality estimates
    LEP = sum(Leslieout$NEAR[["egg_production"]])
    LEPs[k] = LEP
    
    # re-calculate LEP with new adjusted fecundities before k is incorporated
    newfs <- Leslieout$A[1,]/(LEP*adjFec)
    term <- rep(NA,length=maxage)
    for(a in 1:maxage){term[a] <- newfs[a]*Leslieout$NEAR$Survship[a]}
    newLEP[k] <- sum(term) #if you plot 'term', it's the spawning distribution over age
    
    #[OLD WAY] Matrix for simulations analysis: adjust fecundities by multiplying 15/(current LEP) 
    #Leslieout$A[1,] <- Leslieout$A[1,]*(15/LEP)*kvals[k]
    #Lesliearray[,,k] = Leslieout$A 
    
    # Matrix for simulations analysis: adjust fecundities  
    Jacobian.sim[1,] <- (Jacobian.sim[1,]/(LEP*adjFec))
    Lesliearray[,,k] = Jacobian.sim
    
    # Matrix for eigenvalue analysis: relative fecundities at age
    Leslie_eigens <- Jacobian.eig
    #Leslie_eigens[1,] <- Leslie_eigens[1,]*kvals[k] # f*k
    #Leslie_eigens[1,] <- (Leslie_eigens[1,]/sum(Leslie_eigens[1,]))*kvals[k] # rel-f*k
    #Leslie_eigens[1,] <- Leslie_eigens[1,]*(1/LEP)*kvals[k] # f*(1/LEP)*k
    #Leslie_eigens[1,] <- (Leslie_eigens[1,])*(15/LEP)*kvals[k] # f*(15/LEP)*k
    #Leslie_eigens[1,] <- (Leslie_eigens[1,])*(30/LEP)*kvals[k] # f*(15/LEP)*k
    
    ## move survivals to top row, set sub-diag=1
    #Leslie_eigens[1,] <- Leslie_eigens[1,]*Leslie_eigens[2,1]
    #Leslie_eigens[Leslie_eigens == Leslie_eigens[2,1]] <- 1 # suvivals on subdiagonal = 1 (only works if no fishing)
    #Leslie_eigens[1,] <- Leslie_eigens[1,]/sum(Leslie_eigens[1,])*kvals[k]
    #Lesliearray[,,k] = Leslie_eigens
    Leslie_eigens[1,] <- (Leslie_eigens[1,]/(LEP*adjFec))*kvals[k]
    #Leslie_eigens[1,] <- (Leslie_eigens[1,])*kvals[k]
    
    toprowsum[k] <- sum(Leslie_eigens[1,])
    
    # Eigenvalues of matrix with real eigenvalues
    e1[k] = extract_first_eigen_value(Leslie_eigens)
    e2[k] = extract_second_eigen_value(Leslie_eigens)
    e12[k] = e2[k] / e1[k]
    }
  # store max fec values
  maxfecunds[,i] <- maxfec
  LEPvals[,i] <- LEPs
  toprowsums[,i] <- toprowsum
  newLEPs[,i] <- newLEP # newLEPs should be a matrix of all 2s
  # store Leslie matrices
  Aarray[[i]]= Lesliearray #store Leslie 3D array in list of all pops
  NEARsave <- Leslieout$NEAR #save NEAR df
  NEARsave$codNames <- rep(codNames[i],length=length(NEARsave[,1]))
  LeslieoutNEARL[[i]] <- NEARsave
  # store eigenvalues 
  eigenvals1[,i] = e1 #store lambda1 
  eigenvals2[,i] = e2 #store lambda2
  eigenvals12[,i] = e12 #store inverse of damping ratio
  rm(L_inf,K,MG,Mp,theta0,theta1,TEMP,A50,S50,
     maxage,Leslieout,Leslie_eigens,name,NEARsave,
     a,term,maxfec,LEPs,newLEP)
  print(i)
}
rm(e1,e2,e12,i,k,H) #clean up

#format list of NEAR dfs in LeslieoutNEAR
NEAR <- bind_rows(LeslieoutNEARL,id=NULL)

#format dfs, add kvals
eigenvals1 <- as.data.frame(eigenvals1)
eigenvals2 <- as.data.frame(eigenvals2)
eigenvals12 <- as.data.frame(eigenvals12)
names(eigenvals1) <- codNames
names(eigenvals2) <- codNames
names(eigenvals12) <- codNames
e1long <- eigenvals1 %>% mutate(kvals=kvals) %>% gather("codNames","value",1:16) %>% mutate(eigen="e1") 
e2long <- eigenvals2 %>% mutate(kvals=kvals) %>% gather("codNames","value",1:16) %>% mutate(eigen="e2")
e12long <- eigenvals12 %>% mutate(kvals=kvals) %>% gather("codNames","value",1:16) %>% mutate(eigen="e12")
eigendata <- rbind(e1long,e2long,e12long)

# fill in peak, max, sd, cvs in eigendata
eigendata$maxage <- eigentable[match(eigendata$codNames,eigentable$codNames),"max_ages"]
eigendata$peak <- eigentable[match(eigendata$codNames,eigentable$codNames),"mode_age"]
eigendata$cvs <- eigentable[match(eigendata$codNames,eigentable$codNames),"cvs_mode"]
eigendata$sd <- eigentable[match(eigendata$codNames,eigentable$codNames),"sd_mode"]
head(eigendata)

# -----
# use the code below to add 
# -----
# Need to adjust the code below needs to edited for MG and MP
# (1) test significance of relationship between peak age & e1 
# (1a) within each k value: 
# m.e1 <- as.list(rep(NA,length=length(selectedkvals)))
# for(n in 1:length(selectedkvals)){
#   dd <- eigendata[eigendata$kvals == selectedkvals[n] & eigendata$eigen == "e1",]
#   m.e1[[n]] <- lm(dd$value ~ dd$peak)
# }
# rm(n,dd)
# 
# e1_vs_peak_withinkval <- stargazer(m.e1[[1]],m.e1[[2]],m.e1[[3]],m.e1[[4]],m.e1[[5]],dep.var.labels = c("lambda1"),
#                                title="Regression Results for e1~peak within kvals",
#                                covariate.labels=c("peak age"),type="text", report='vc*p',
#                                column.labels = c("k=0.1","k=0.3","k=0.5","k=0.71","k=0.91"),
#                                out="C:/Users/provo/Documents/GitHub/popdy/cod_figures/model_e1vspeak_withink.txt")
# # (1b) across k values --> average within k, for each pop
# # (1b1) average within k
# meanatk_e1 <- rep(NA,length=length(kvals))
# for(n in 1:length(kvals)){
#   meanatk_e1[n] <- mean(eigendata[eigendata$kvals == kvals[n] & eigendata$eigen == "e1",]$value)}
# e1_vs_peak_acrossk <- summary(lm(meanatk_e1~kvals))
# 
# # (2) test significance of relationship between peak age & e12
# m.e12 <- as.list(rep(NA,length=length(selectedkvals)))
# for(n in 1:length(selectedkvals)){
#   dd <- eigendata[eigendata$kvals == selectedkvals[n] & eigendata$eigen == "e12",]
#   m.e12[[n]] <- lm(dd$value ~ dd$peak)
# }
# e12_vs_peak_withinkval <- stargazer(m.e12[[1]],m.e12[[2]],m.e12[[3]],m.e12[[4]],m.e12[[5]],
#                                     dep.var.labels = c("damping ratio"),
#                                    title="Regression Results for e12~peak within kvals",
#                                    covariate.labels=c("peak age"),type="text", report='vc*p',
#                                    column.labels = c("k=0.1","k=0.4","k=0.6","k=0.9"),
#                                    out="C:/Users/provo/Documents/GitHub/popdy/cod_figures/model_e1vspeak_withink.txt")
# meanatk_e12 <- rep(NA,length=length(kvals))
# for(n in 1:length(kvals)){
#   meanatk_e12[n] <- mean(eigendata[eigendata$kvals == kvals[n] & eigendata$eigen == "e12",]$value)}
# e12_vs_peak_acrossk <- summary(lm(meanatk_e12~kvals))



# *************************************** #
# (Fig 3ab) Eigenvalue 4 panel plot AND Gantt plot
# *************************************** #

lambda1 <- ggplot(eigendata[eigendata$eigen=="e1" & eigendata$kvals %in% selectedkvals,],
                  aes(x=maxage,y=value)) +
  geom_point() + 
  geom_smooth(method="lm",se=FALSE,color="black") +
  facet_grid(. ~ kvals) +
  scale_y_continuous(limits=c(0.6,1)) +
  geom_text_repel(data=eigendata[eigendata$eigen=="e1"& eigendata$kvals %in% selectedkvals,],
                  aes(label = codNames),
                  segment.color = "grey",
                  size = 2,
                  na.rm = TRUE) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 40, hjust = 1)) +
  ylab("a") +
  xlab("Maximum age") +
  ylab(expression(paste(lambda[1]))) 

lambda12 <- ggplot(eigendata[eigendata$eigen=="e12" & eigendata$kvals %in% selectedkvals,],
                   aes(x=cvs,y=value)) +
  geom_point() + 
  geom_smooth(method="lm",se=FALSE,color="black") +
  facet_grid(. ~ kvals) +
  scale_y_continuous(limits=c(0,1)) +
  geom_text_repel(data=eigendata[eigendata$eigen=="e12"& eigendata$kvals %in% selectedkvals,],
                  aes(label = codNames),
                  segment.color = "grey",
                  size = 2,
                  na.rm = TRUE) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 40, hjust = 1)) +
  ylab("b") +
  xlab("Spawning biomass distribution Stdev") +
  ylab(expression(paste(frac(abs(lambda[2]),lambda[1])))) 


tiff(file='C:/Users/Mikaela/Documents/GitHub/popdy/cod_figures/manuscript3/fig3_ab_adjustLEP_Max_&_CV.tiff', units="in", width=6, height=7, res=300)
grid.newpage()
grid.draw(rbind(ggplotGrob(lambda1), ggplotGrob(lambda12), size = "last"))
dev.off()
rm(p,lambda1,lambda12)


# *************************************** #
# (Fig 1:schematic) Choose alpha values, plot BH curves
# *************************************** #
# In my analysis I want to evaluate population spectra when
# populations sit on different slopes of the BH curve. 
# To change the position on the curve I change beta while keeping
# alpha constant. As beta increases, the slope at equilibrium
# increases (the population moves down the curve). 
# Since LEP is 1 for all populations, equilibrium occurs where
# recruits = egg production (where the 1:1 line crosses the BH.

# For each alpha value, plot the BH curve and the 1:1 line
BH <- function(alpha,beta,E){ E/((1/alpha)+(E/beta))}

a = seq(from=0.97,to=10,by=0.01) #alpha values
a = c(5.51,2.07,1.38,0.97)
b = 1000 # beta
ee <- seq(from=0,to=b,by=b/10000) #range of egg values

#Calculate recruit values for each alpha
Rlist <- as.list(rep(NA,length=length(a))) 
names(Rlist) <- a
for(i in 1:length(a)){ Rlist[[i]] <- BH(alpha=a[i],beta=b,E=ee) }
RvE <- as.data.frame(do.call(cbind,Rlist))
RvE$eggs <- ee
RvElong <- RvE %>% gather(alpha,value,1:length(a))
RvElong$alpha <- as.numeric(RvElong$alpha)

# For each alpha, calculate the intersection with the 1:1 line. Find when R=E
eq <- rep(NA,length=length(a))
interpt <- rep(NA,length=length(a))
slopes <- rep(NA,length=length(a))
for(i in 1:length(a)){
  ee <- seq(from=0,to=b,by=b/10000) #all possible egg values
  rr <- BH(alpha=a[i],beta=b,E=ee) #calc recruit values using BH
  rr.lep <- 0+round(1/conLEP,digits=2)*ee #calc recruit values using 1/LEP
  dat <- as.data.frame(cbind(ee,rr,rr.lep))
  names(dat) <- c('eggs','recruits','recruits.lep')
  dat$diff <- abs(dat$recruits - dat$recruits.lep) #diff between rr and rr.lep
  dat <- dat[-1,] #remove first row at origin
  row.names(dat) <- NULL
  #head(dat)
  #plot(dat$diff[1:1000],type="l")
  eq[i] <- which.min(dat$diff )
  interpt[i] <- dat[eq[i],]$eggs
  #dat[(eqpoint[b]-3):(eqpoint[b]+3),]
  # calculate slope around equal point
  x1 <- dat[(eq[i]-1),]$eggs
  y1 <- dat[(eq[i]-1),]$recruits
  x2 <- dat[(eq[i]+1),]$eggs
  y2 <- dat[(eq[i]+1),]$recruits
  slopes[i] <- round((y2-y1)/(x2-x1),digits=2)
  print(i)
}
datalpha <- as.data.frame(cbind(a,interpt,slopes))
plot(x=datalpha$a,y=datalpha$slopesR)

#plot_these_alpha <- datalpha %>% filter(a %in% c(1.1,"1.4",2.0,"3.3",10.0)) #weird work around :(
plot_these_alpha <- datalpha %>% filter(a %in% c("5.51","2.07","1.38","0.97")) #weird work around :(

RvElong_forplotting <- RvElong[RvElong$alpha %in% c("5.51","2.07","1.38","0.97"),]
RvElong_forplotting$interpt <- plot_these_alpha[match(RvElong_forplotting$alpha,plot_these_alpha$a),"interpt"]
RvElong_forplotting$alpha <- as.character(RvElong_forplotting$alpha)
str(RvElong_forplotting)
# add column with intercept point for each alpha
LEPlineslope = round(1/conLEP,digits=2)
top <- ggplot(data=RvElong_forplotting,aes(x=eggs,y=value,linetype=alpha)) + geom_line() +
  geom_dl(aes(label=alpha),method="last.points") +
  xlim(c(0,1100)) +
  geom_abline(intercept=0, slope=LEPlineslope, color="black",size=1) + theme_classic() + 
  ylab("Recruits") + xlab("Egg production") +
  geom_point(aes(x=interpt,y=LEPlineslope*interpt,group=alpha),size=3) +
  theme(legend.position = "none") 
# # FIG 1(b) - Schematic showing slope at intersection for multiple alpha values
# bottom <- ggplot(data=datalpha,aes(x=a,y=slopes)) + geom_line() +
#   geom_point(data=plot_these_alpha,aes(x=a,y=slopes),size=3) +
#   geom_text(data=plot_these_alpha,aes(label=a),nudge_x=0.3,nudge_y=0.04) +
#   xlab(expression(alpha)) +
#   ylab("Slope of the Beverton-Holt\n curve at the intersection with 1/LEP") +
#   theme_classic()
# export figure

tiff(file='C:/Users/Mikaela/Documents/GitHub/popdy/cod_figures/manuscript3/fig1_schematic.tiff', units="in", width=4, height=4, res=300)
#do.call(grid.arrange,c(plist,ncol=1))
top
dev.off()
rm(top)

# # -------
# # Alternative Fig for schematic - showing one BH curve and multiple 1/LEP
# # -------
# b=1000
# a=2
# ee <- seq(from=0,to=b*10,by=b/10000) #range of egg values
# rr <- BH(alpha=a,beta=b,E=ee)
# dd <- as.data.frame(cbind(ee,rr))
# # calculate x and y end points for segments on plot (1/LEP line)
# y1 = rep(0,length=length(selectedalphas))
# y2 = rep(1000,length=length(selectedalphas))
# x1 = rep(0,length=length(selectedalphas))
# slopes_to_plot <- datalpha %>% filter(a %in% c(1.1,"1.4",2.0,"3.3",10.0)) %>% select(slopes) 
# x2 = c(y2/ slopes_to_plot )
# x2 <- x2[['slopes']]
# df <- data.frame(cbind(y1,y2,x1,x2))
# 
# # this plot needs fixing: put black circle at intersection point
# ggplot(dd,aes(x=ee,y=rr)) + geom_line() +
#   ylim(0,1500) +
#   geom_segment(data=df,aes(x=x1,y=y1,xend=x2,yend=y2)) +
#   theme_classic() +
#   xlab("Egg production" ) +
#   ylab("Recruits") +
#   geom_point(data=plot_these_alpha,aes(x=a,y=interpt),size=3)
# 
# rm(y1,y2,x1,x2,slopes_to_plot,df,dd,b,a,ee,rr)


# *************************************** #
# (3) Simulate pops. Loop over Aarray list to simulate using different Leslie matrices 
# *************************************** #
# set params for simulation:
# first, calculate alpha for each value of k
# k = 1/(alpha*LEP^2) --> we adjusted LEP above to be 2 for all pops
#alphas = round(1/(kvals*conLEP^2),digits=2)
#alphas = unique(alphas)
timesteps = 1000 #need this now to create
rm_first_timesteps = 200
betas = 1000
#alphas <- seq(from=1.1, to=10, by=0.01) #these correspond to 
sig_r = 0.3


output.3d.list <- as.list(rep(NA,length=length(codNames))) #store timeseries here
names(output.3d.list) <- codNames

for (i in 1:length(Aarray)) { #step through each pop
  Leslie3d = Aarray[[i]] #select the 3d array of Leslie matricies
  # array dims: row=ts length, col=4 is number of ts (eggs,recruits,Nt,Nsize), depth=F vals
  output.matrix <- array(NA,c(timesteps-1,3,length(alphas))) 
  
  #Run this loop if you vary beta and keep alpha constant
  #for (b in 1:length(betas)) { #step through each Leslie matrix (for each F value)
  #  output = sim_model(A=Leslie3d[,,1], timesteps=timesteps, 
  #                     alpha=alphas[i], beta=betas[b], 
  #                     sig_r=sig_r, initial_eggs=betas[b])
  #  
  #  length(output$Nsize) <- length(output$N_t) #trim Nsize ts vector, -2 elements
  #  output.matrix[,,b] <- do.call(cbind,output) #fill in array for pop i
    #colnames(output.matrix) <- names(output)
  #}
  
  #Run this loop if you vary alpha and keep beta constant
  for (a in 1:length(alphas)) { #step through each Leslie matrix (for when k=1)
    output = sim_model(A=Leslie3d[,,1], timesteps=timesteps, 
                       alpha=alphas[a], beta=betas, 
                       sig_r=sig_r, initial_eggs=betas)
    
    output$eggs <- output$eggs[1:(timesteps-1)] #trim ts vector, -1 
    output$Nsize <- output$Nsize[1:(timesteps-1)]
    output.matrix[,,a] <- do.call(cbind,output) #fill in array for pop i
    #colnames(output.matrix) <- names(output)
  }
  
  output.3d.list[[i]] <- output.matrix
  print(i)
}
rm(i,a,Leslie3d,output.matrix,output) #clean up

# At this point I have one important object:
# 1. [output.3d.list] a list of 3d arrays. Each array is timeseries output
#    from simulations at different alpha levels. 


# *************************************** #
# (3) Format output ts for plotting simulations using output.3d.list
# *************************************** #
variable_type <- c("eggs","recruits","Nsize")

# --- reorganize egg timeseries data --- #
#var.number <- 2 # eggs
#df.list <- as.list(rep(NA,length=length(codNames)))
#names(df.list) <- codNames
#for (i in 1:length(output.3d.list)) {
#  # first, reformat data to work with ggplot
#  aa <- as.data.frame(output.3d.list[[i]][,var.number,])
#  aa$year <- seq(from=1, to=length(aa[,1]))
#  colnames(aa) <- c(alphas,"year")
#  aa1 <- aa %>% gather(alphavalue,value,1:length(alphas))
#  aa1$variable <- rep(variable_type[var.number],length=length(aa1[,1]))
#  aa1$codNames <- rep(codNames[i],length=length(aa[,1]))
#  df.list[[i]] <- aa1
#  rm(aa1,aa)}
#eggs.ts <- do.call(rbind,df.list)
#eggs.ts$alphavalue <- as.numeric(as.character(eggs.ts$alphavalue))
#rm(df.list,i) #clean up

# --- reorganize recruit timeseries data --- #
var.number <- which(variable_type == "recruits") # recruits
df.list <- as.list(rep(NA,length=length(codNames)))
names(df.list) <- codNames
for (i in 1:length(output.3d.list)) {
  # first, reformat data to work with ggplot
  aa <- as.data.frame(output.3d.list[[i]][,var.number,])
  aa$year <- seq(from=1, to=length(aa[,1]),by=1)
  colnames(aa) <- c(alphas,"year")
  aa1 <- aa %>% gather(alphavalue,value,1:length(alphas))
  aa1$variable <- rep(variable_type[var.number],length=length(aa1[,1]))
  aa1$codNames <- rep(codNames[i],length=length(aa[,1]))
  aa1$alphavalue <- as.numeric(as.character(aa1$alphavalue))
  df.list[[i]] <- aa1
  print(i)
  }
recruits.ts <- bind_rows(df.list,id=NULL)
#recruits.ts$alphavalue <- as.numeric(as.character(recruits.ts$alphavalue))
str(recruits.ts)
rm(df.list,i,aa1,aa,var.number) #clean up

# --- reorganize Nsize timeseries data --- #
#var.number <- 4 # Nsize
#df.list <- as.list(rep(NA,length=length(codNames)))
#names(df.list) <- codNames
#for (i in 1:length(output.3d.list)) {
#  # first, reformat data to work with ggplot
#  aa <- as.data.frame(output.3d.list[[i]][,var.number,])
#  aa$year <- seq(from=1, to=length(aa[,1]))
#  colnames(aa) <- c(alphas,"year")
#  aa1 <- aa %>% gather(alphavalue,value,1:length(alphas))
#  aa1$variable <- rep(variable_type[var.number],length=length(aa1[,1]))
#  aa1$codNames <- rep(codNames[i],length=length(aa[,1]))
#  df.list[[i]] <- aa1
#  }
#nsize.ts <- do.call(rbind,df.list)
#rm(df.list,i,aa1,aa) #clean up

#ts.data <- rbind(eggs.ts,recruits.ts) #combine data -- too big!
ts.data <- recruits.ts
rownames(ts.data) <- NULL
ts.data$peak <- eigentable[match(ts.data$codNames,eigentable$codNames),"mode_age"]
ts.data$cvs <- eigentable[match(ts.data$codNames,eigentable$codNames),"cvs_mode"]
head(ts.data)
rm(recruits.ts)

# Appendix - time series CV vs alpha
# For each cod population, calculate CV of time series (for recruits)
# and plot it against the corresponding beta value
pp <- as.list(rep(NA,length=length(codNames)))
for(i in 1:length(codNames)){
  dat <- ts.data[ts.data$codNames == codNames[i],]
  cvs <- rep(NA,length=length(alphas)) #store cv for each alpha here
  for(b in 1:length(alphas)){ #step through alpha vals
    vals <- dat[dat$alphavalue == alphas[b],]$value[rm_first_timesteps:(timesteps-2)]#rm initial tsteps
    cvs[b] <- sd(vals)/mean(vals) #cv=sd/mean
  }
  cvdat <- as.data.frame(cbind(cvs,alphas))
  pp[[i]] <- ggplot(cvdat,aes(x=alphas,y=cvs)) + geom_line() + 
    ggtitle(paste(codNames[i])) +
    ylab("") + xlab("") + ylim(c(0,0.25)) + xlim(c(0,6))
  print(i)
}

tiff(file='C:/Users/Mikaela/Documents/GitHub/popdy/cod_figures/manuscript3/SI/figS2_cv_vs_alpha.tiff', 
     units="in", width=7, height=8, res=300) 
do.call(grid.arrange,c(pp,ncol=4,left="time series CV",bottom="alpha value"))
dev.off()
rm(pp,i,cvdat,dat,cvs)

# *************************************** #
# (4) Now that timeseries data is formated, let's plot! --- This section is option, skip to section 5 to calculate spectra.
# *************************************** #

# Fig S3:
#plot recruitment - one plot per pop, similar to egg plots
#time series for 4 k values 
selectedalphasP <- c(2.07,1.38,0.97)
str(ts.data)
p.forappendix <- as.list(rep(NA,length=length(codNames_ordered_by_peak)))
p.formanuscript <- as.list(rep(NA,length=length(codNames_ordered_by_peak)))
names(p.formanuscript) <- codNames_ordered_by_peak

for (i in 1:length(codNames_ordered_by_peak)){
  dd <- ts.data[ts.data$variable == "recruits" & 
                  ts.data$codNames == codNames_ordered_by_peak[i] &
                  ts.data$year %in% seq(from=rm_first_timesteps,to=(timesteps-500),by=1) &
                  ts.data$alphavalue %in% selectedalphasP,]
  dd$kval <- factor(round(1/(dd$alphavalue*conLEP^2),digits=2))
  
  #subtract mean from ts, and adjust so all values are positive (makes it easier to interpet)
  dd$valueminusmean <- (dd$value-mean(dd$value))+500 
  
  p.forappendix[[i]] <- ggplot(dd,aes(x=year,y=valueminusmean,color=kval)) +
    xlab("") + ylab("") +
    geom_line() + theme_classic() +
    #scale_color_brewer(palette = "Reds") +
    scale_y_log10(limits=c(200,1000)) +
    ggtitle(paste(codNames_ordered_by_peak_plot[i]," (peak=",dd$peak[1]," cv=",round(dd$cvs[1],digits=2),")",sep="")) 
  
  p.formanuscript[[i]] <- ggplot(dd,aes(x=year,y=valueminusmean,color=kval)) +
    xlab("") + ylab("") +
    geom_line() + theme_classic() +
    scale_color_brewer(palette = "Reds") +
    theme(legend.position = "none") + ggtitle("") +
    scale_y_log10(limits=c(200,1000))
  
  print(i)
}
tiff(file='C:/Users/Mikaela/Documents/GitHub/popdy/cod_figures/manuscript3/SI/figS3_timeseries_recruits_sigmaR0.3.tiff', units="in", width=7, height=13, res=300) 
do.call(grid.arrange,c(p.forappendix,ncol=2,left="Recruits (before noise)", bottom="Year"))
dev.off()
rm(p.forappendix,dd)

# --------------
# Fig 2abc: 3 spawning biomass distributions, 3 recruit ts, 1 spectrum plot for all pops 

# Fig2a: 3 spawning biomass distributions
source("C:/Users/Mikaela/Documents/GitHub/popdy/cod_code/3_plot_spawing_biomass_distributions.r") # spawning distributions
pDc <- pMG$Celtic
pDw <- pMG$W_Baltic
pDi <- pMG$Iceland

# Fig2b: 3 recruit ts
dd <- ts.data[ts.data$variable == "recruits" & 
                ts.data$codNames %in% c("Celtic","W_Baltic","Iceland") &
                ts.data$year %in% seq(from=rm_first_timesteps,to=(timesteps-500),by=1) &
                ts.data$alphavalue %in% selectedalphasP,]
dd$kval <- factor(round(1/(dd$alphavalue*conLEP^2),digits=2))
#subtract mean from ts, and adjust so all values are positive (makes it easier to interpet)
dd$valueminusmean <- (dd$value-mean(dd$value))+500 


pTc <- ggplot(data=dd[dd$codNames=="Celtic",],
              aes(x=year,y=valueminusmean,color=kval)) +
  xlab("Year") + ylab("log(recruits)") +
  geom_line() + theme_classic() +
  scale_color_manual(values = c("blue2","blue3","blue4")) +
  theme(legend.position = "none",
        axis.text.y=element_blank()) + 
  ggtitle("") +  scale_y_log10(limits=c(200,1000))

pTw <- ggplot(data=dd[dd$codNames=="W_Baltic",],
              aes(x=year,y=valueminusmean,color=kval)) +
  xlab("Year") + ylab("log(recruits)") +
  geom_line() + theme_classic() +
  scale_color_manual(values = c("seagreen2","seagreen3","seagreen4")) +
  theme(legend.position = "none",
        axis.text.y=element_blank()) +
  ggtitle("") + scale_y_log10(limits=c(200,1000))

pTi <- ggplot(data=dd[dd$codNames=="Iceland",],
              aes(x=year,y=valueminusmean,color=kval)) +
  xlab("Year") + ylab("log(recruits)") +
  geom_line() + theme_classic() +
  scale_color_manual(values = c("red2","red3","red4")) +
  theme(legend.position = "none",
        axis.text.y=element_blank()) + 
  ggtitle("") + scale_y_log10(limits=c(200,1000))



# *************************************** #
# (5) Calculate frequency content from timeseries
# *************************************** #
# Plan:
# 1. Walk through each cod pop, do spectral analysis at alpha levels
# 2. Store spec values for eggs, recruits, and Nsize

# 1. Walk through each cod pop, do spectral analysis at alpha levels
#sp.eggsL <- as.list(rep(NA,length=length(codNames))) #object for spec analysis  
sp.recruitLsm <- as.list(rep(NA,length=length(codNames))) 
span.multiplier= 1.5
# testing small number of alphas for saving time, run only the alphas associated
# with plotting the grid plots
#alphas.grid <- c(8.26,7.51,6.36,5.51,4.59,3.31,2.07,1.38,0.97)
#alphas <- alphas.grid



for (i in 1:length(codNames)){
  
  ts <- ts.data[ts.data$codNames == codNames[i],] #subset data for pop i 
  
  # setting 'span' - a vector of odd integers to specify the smoothers
  tmp <- ceiling(sqrt(length(1:(timesteps-rm_first_timesteps-1)))) #sq root of timeseries lgth, rounded
  if (tmp %% 2 == 0) {m <- tmp+1} else {m <- tmp} #make it odd, if the square root is even
  m = m * span.multiplier
  
  # --- spectral analysis on EGGS --- #
  #spsaveL <- as.list(rep(NA,length=length(alphas)))
  #names(spsaveL) <- alphas
  #for (b in 1:length(alphas)){
  #  xx = ts[ts$variable == "eggs" & ts$alphavalue == alphas[b],]$value[rm_first_timesteps:(timesteps-2)] - mean(ts[ts$variable == "eggs" & ts$alphavalue == alphas[b],]$value[rm_first_timesteps:(timesteps-2)])
  #  sp = spectrum(x=xx,spans=c(m,m),plot = FALSE)
    #save spec output for plotting, Helen Wearing says to multiply by 2
  #  spsaveL[[b]] <- sp$spec*2
  #}
  #spsave <- as.data.frame(do.call(cbind,spsaveL))
  #spsave$freq <- sp$freq
  #spsavelong <- spsave %>% gather(alphavalue, value, 1:length(alphas))
  #spsavelong$codNames <- rep(codNames[i],length=length(spsavelong[,1]))
  #sp.eggsL[[i]] <- spsavelong
  #rm(spsave)
  
  # --- spectral analysis on RECRUIT --- #
  spsaveL <- as.list(rep(NA,length=length(alphas)))
  names(spsaveL) <- alphas
  
  for (b in 1:length(alphas)){
    yy = ts %>% filter(round(alphavalue,digits=2) == round(alphas[b],digits=2)) %>% select(value) %>% slice(rm_first_timesteps:(timesteps-2)) 
    #meanyy = mean(rawyy$value)
    #yy = rawyy$value-meanyy
    
    # ---
    # original code here to calculate AUC (spec*freq)
    #sp = spec.pgram(yy$value,spans=c(m,m),taper=0.1,plot = FALSE,demean = TRUE)
    #auc = sum((sp$spec) * sp$freq) #total area under curve
    #spsaveL[[b]] = sp$spec/auc #divide each spec by auc, makes auc=1
    # ---
    # code below is new stuff, hopefully this works:
    sp = spec.pgram(yy$value,spans=c(m,m),taper=0.1,plot = FALSE,demean = TRUE)
    auc = sum(sp$spec) #total area under curve
    spsaveL[[b]] = sp$spec/auc #divide each spec by auc, makes auc=1
    # ---
    print(b)
  }
  spsave <- as.data.frame(do.call(cbind,spsaveL))
  spsave$freq <- sp$freq
  spsavelong <- spsave %>% gather(alphavalue, value, 1:length(alphas))
  spsavelong$codNames <- rep(codNames[i],length=length(spsavelong[,1]))
  sp.recruitLsm[[i]] <- spsavelong
  rm(spsave,yy,sp)
  print(i)
}
#eggs <- do.call(rbind,sp.eggsL)
#eggs$variable.type <- rep("eggs",length=length(eggs$freq))
# rec <- bind_rows(sp.recruitL,id=NULL)
# rec$variable.type <- rep("recruits",length=length(rec$freq))

recsm <- bind_rows(sp.recruitLsm,id=NULL)
recsm$variable.type <- rep("recruits",length=length(recsm$freq))

specdatalong <- recsm
head(specdatalong)
#specdatalong$alphavalue <- factor(specdatalong$alphavalue,levels=selectedalphas)
specdatalong$peak <- eigentable[match(specdatalong$codNames,eigentable$codNames),"mode_age"]
specdatalong$cvs <- eigentable[match(specdatalong$codNames,eigentable$codNames),"cvs_mode"]
specdatalong$cvs <- round(specdatalong$cvs,digits=2)
str(specdatalong)
rm(recsm)
# *************************************** #
# (6) Plot spectral analysis 
# *************************************** #
# --- egg spectra --- # 
#peggs_sp <- list()
#for (i in 1:length(codNames_ordered_by_peak)){
  
#  dataforplot <- specdatalong[specdatalong$variable.type == "eggs" 
#                              & specdatalong$codNames==codNames_ordered_by_peak[i] &
#                                specdatalong$alphavalue %in% selectedalphas,]
  # store plots in list
#  peggs_sp[[i]] <- ggplot(data=dataforplot, aes(x=freq,y=value,group=alphavalue)) + 
#    geom_line(aes(color=alphavalue)) + #ylim(1,10) +
#    geom_vline(xintercept = (1/eigentable[eigentable$codNames == codNames_ordered_by_peak[i],]$mode_age),
#             linetype="dotted") +
#    ggtitle(paste(codNames_ordered_by_peak_plot[i])) + 
#    theme(plot.title = element_text(size = 1)) + 
#    theme_classic() + ylab("") + xlab("") + theme(legend.position = "none")
#}
#rm(i,dataforplot)
#do.call(grid.arrange,c(peggs_sp,ncol=4,left="Spectra",bottom="Frequency"))
#dev.off()


# --- recruits spectra --- #
prec_sp <- list()
cs <- scales::seq_gradient_pal("#bdd7e7", "#08519c", "Lab")(seq(0,1,length.out=length(selectedalphas)))
for (i in 1:length(codNames_ordered_by_peak)){
  dataforplot <- specdatalong[specdatalong$variable.type == "recruits" 
                              & specdatalong$codNames==codNames_ordered_by_peak[i] &
                                specdatalong$alphavalue %in% selectedalphas,]
  #add col for k vals, label plots with this
  dataforplot$k_value <- round(1/(as.numeric(dataforplot$alphavalue) * conLEP^2),digits=2) 
  dataforplot$k_value <- as.character(dataforplot$k_value)
  # store plots in list
  prec_sp[[i]] <- ggplot(data=dataforplot, aes(x=freq,y=value,group=k_value)) + 
    geom_line(aes(color=k_value)) + 
    #scale_y_log10(limits = c(0.001,2500),breaks=c(1,100,1000)) +
    geom_vline(xintercept = (1/eigentable[eigentable$codNames == codNames_ordered_by_peak[i],]$mode_age),
               linetype="dotted") +
    ggtitle(paste(codNames_ordered_by_peak_plot[i],
                  " peak=",dataforplot$peak[1],
                  " cv=",dataforplot$cvs[1],sep="")) + 
    theme_classic() + ylab("") + xlab("") + #theme(legend.position = "none") +
    #scale_colour_manual(values=cs) +
    scale_color_brewer(palette = "Reds") +
    theme(plot.title = element_text(size = 10)) +
    scale_y_log10() 
    
    
}
names(prec_sp) <- codNames_ordered_by_peak
rm(i,dataforplot)

tiff(file='C:/Users/Mikaela/Documents/GitHub/popdy/cod_figures/manuscript3/SI/figS4_spectra_multipanel_freqlog.tiff', units="in", width=11, height=7, res=300) 
do.call(grid.arrange,c(prec_sp,ncol=4,left="Variance",bottom="Frequency"))
dev.off()

# ---
# Fig X. 3x3 plot (PEAK)
# dataforplot <- specdatalong[specdatalong$codNames %in% c("Iceland","Northsea","GB") &
#                               specdatalong$alphavalue %in% c(10,2.5,1.1),]
# dataforplot$skvals <- round((1/as.numeric(as.character(dataforplot$alphavalue))),digits = 2)
# dataforplot$skvals <- factor(dataforplot$skvals,levels = unique(dataforplot$skvals))
# dataforplot$codNames <- factor(dataforplot$codNames,levels = c("Iceland","Northsea","GB"))
# dataforplot$peak <- eigentable[match(dataforplot$codNames, eigentable$codNames),"mode_age"]
# dataforplot$line <- 1/(dataforplot$peak)
# dataforplot$linehalf <- dataforplot$line/2
# levels(dataforplot$codNames) <- c("Iceland (peak=10,CV=0.24)","Northsea (peak=5,CV=0.46)","GB (peak=3,CV=0.61)")
# 
# tiff(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/manuscript/Fig3_3by3_peak.tiff', units="in", width=7, height=7, res=300) 
# ggplot(data=dataforplot,aes(x=freq,y=log10(value))) +
#   geom_line() + facet_grid(skvals ~ codNames) +
#   theme_bw() + 
#   geom_vline(data=dataforplot,aes(xintercept = line),linetype="dotted") +
#   geom_vline(data=dataforplot,aes(xintercept = linehalf),linetype="dashed") +
#   ggtitle("Dashed line = 1/(2T)    Dotted line = 1/T")
# dev.off()
# rm(dataforplot)
# 
# # Fig X. 3x3 plot: pick populations with range of CV
# # pick: NE_Arctic, Northsea (or W_Baltic), W_Scotland (or GM)
# dataforplot <- specdatalong[specdatalong$codNames %in% c("NE_Arctic","W_Baltic","W_Scotland") &
#                               specdatalong$alphavalue %in% c(10,2,1.1),]
# dataforplot$skvals <- round((1/as.numeric(as.character(dataforplot$alphavalue))),digits = 2)
# dataforplot$skvals <- factor(dataforplot$skvals,levels = c(0.1,0.5,0.91))
# dataforplot$codNames <- factor(dataforplot$codNames,levels = c("NE_Arctic","W_Baltic","W_Scotland"))
# dataforplot$peak <- eigentable[match(dataforplot$codNames, eigentable$codNames),"mode_age"]
# dataforplot$line <- 1/(dataforplot$peak)
# dataforplot$linehalf <- dataforplot$line/2
# levels(dataforplot$codNames) <- c("NE_Arctic (CV=0.21,peak=9)","W_Baltic (CV=0.45,peak=5)","W_Scotland (CV=0.62,peak=3)")
# tiff(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/manuscript/Fig3_3by3_cvs.tiff', units="in", width=7, height=7, res=300) 
# ggplot(data=dataforplot,aes(x=freq,y=log10(value))) +
#   geom_line() + facet_grid(skvals ~ codNames) +
#   theme_bw() + geom_vline(data=dataforplot,aes(xintercept = line),linetype="dotted") +
#   geom_vline(data=dataforplot,aes(xintercept = linehalf),linetype="dashed") +
#   ggtitle("Dashed line = 1/(2T)    Dotted line = 1/T")
# dev.off()
# rm(dataforplot)
# 

# choose only pops with peak age at 3 or 4
#t <- textGrob("")
#lowpeakage <- list(prec_sp$Celtic,prec_sp$W_Scotland, #2 pops (peak=3)
#                   prec_sp$Faroe, prec_sp$GB, prec_sp$GM, prec_sp$Kat ) #4 pops (peak=4)
#medpeakage <- list(prec_sp$cod3M, prec_sp$Northsea, prec_sp$W_Baltic, #3 pops (peak=5)
#                   prec_sp$Coas,t,t) #1 pop (peak=6)
#highpeakage <- list(prec_sp$cod2J3KL,prec_sp$cod3NO,prec_sp$cod3Ps,prec_sp$Iceland,prec_sp$NGulf, #5 pops (peak=8)
#                    prec_sp$NE_Arctic) #1 pop (peak=9)
#export low peak age
#tiff(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/1_lowpeakage_panel.tiff', 
#     units="in", width=3, height=5, res=300) 
#do.call(grid.arrange,c(lowpeakage,ncol=2,left="Spectra",top="Low peak spawning (age = 3, 4)"))
#dev.off()
# export med peak age
#tiff(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/2_medpeakage_panel.tiff', 
#     units="in", width=3, height=5, res=300) 
#do.call(grid.arrange,c(medpeakage,ncol=2,left="Spectra",top="Medium peak spawning (age = 5, 6)"))
#dev.off()
# export high peak age
#tiff(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/3_highpeakage_panel.tiff', 
#     units="in", width=3, height=5, res=300) 
#do.call(grid.arrange,c(highpeakage,ncol=2,left="Spectra",top="High peak spawning (age = 8, 9)"))
#dev.off()


# *************************************** #
# (Fig 2c) Spectra on one plot, color coded by peak
# spawning age, max age, and CV of spawning biomass distribution
# *************************************** #
# --- recruits spectra: all on one plot ---#
dataforplot <- specdatalong[specdatalong$variable.type == "recruits" &
                              specdatalong$alphavalue %in% alphas[2],] # select one alpha value
dataforplot$peak <- eigentable[match(dataforplot$codNames, eigentable$codNames),"mode_age"]
dataforplot$maxage <- eigentable[match(dataforplot$codNames, eigentable$codNames),"max_ages"]
dataforplot$sd_mode <- eigentable[match(dataforplot$codNames, eigentable$codNames),"sd_mode"]
dataforplot$cvs_mode <- eigentable[match(dataforplot$codNames, eigentable$codNames),"cvs_mode"]
dataforplot$cvs_mode <- round(dataforplot$cvs_mode,digits = 2)
dataforplot$codNames_peak <- paste(dataforplot$codNames,"(",dataforplot$peak,")",sep="")
dataforplot$codNames_maxage <- paste(dataforplot$codNames,"(",dataforplot$maxage,")",sep="")
dataforplot$codNames_cv <- paste(dataforplot$codNames,"(",dataforplot$cvs_mode,")",sep="")
dataforplot$codNames <- factor(dataforplot$codNames, 
                               levels=c("Coas","cod3M","cod3NO","cod3Ps","Northsea",
                                 "Faroe","GB","GM","Iceland","Kat","NGulf","W_Scotland",
                                 "Celtic","NE_Arctic","cod2J3KL","W_Baltic"))
dataforplot$cvs_cat <- rep(NA,length=length(dataforplot$freq))
dataforplot[dataforplot$cvs < 0.3,]$cvs_cat <- 1 #"narrow CV (0.19-0.27)"
dataforplot[dataforplot$cvs > 0.5,]$cvs_cat <- 3 #"wide CV (0.51-0.67)"
dataforplot[is.na(dataforplot$cvs_cat),]$cvs_cat <- 2  #"medium CV (0.32-0.46)"
#dataforplot[dataforplot$maxage < 10,]$cvs_cat <- 1 #"short-lived"
#dataforplot[dataforplot$maxage > 12,]$cvs_cat <- 3 #"long-lived"
#dataforplot[is.na(dataforplot$cvs_cat),]$cvs_cat <- 2 #"medium-lived"
dataforplot$cvs_cat <- factor(dataforplot$cvs_cat)
dataforplot <- subset(dataforplot,selec=c(freq,value,codNames_cv,cvs_cat))
head(dataforplot)

colourCount = length(unique(dataforplot$codNames_cv))
getPalette = colorRampPalette(brewer.pal(9, "Paired"))

# calc average spectra line for the three CV bins
dataforplot1 <- dataforplot %>% select(freq,value,codNames_cv) %>% spread(codNames_cv,value)

narrow <- rowMeans(subset(dataforplot1,
                          select=c("NE_Arctic(0.19)","Coas(0.22)","cod3NO(0.27)","Iceland(0.24)")))
medium <- rowMeans(subset(dataforplot1,
                          select=c("W_Scotland(0.53)","Celtic(0.67)","Kat(0.63)","Faroe(0.61)","GM(0.51)","GB(0.61)","cod2J3KL(0.52)")))
wide <-   rowMeans(subset(dataforplot1,
                          select=c("cod3M(0.32)","Northsea(0.46)","NGulf(0.4)","cod3Ps(0.35)","W_Baltic(0.46)")))

d1 <- as.data.frame(cbind(dataforplot1$freq,narrow))
colnames(d1) <- c("freq","narrow")
d1 <- d1 %>% gather(codNames_cv,value,2:(length(d1[1,])))
d1$cvs_cat <- rep("4",length=length(d1[,1]))

d2 <- as.data.frame(cbind(dataforplot1$freq,medium))
colnames(d2) <- c("freq","medium")
d2 <- d2 %>% gather(codNames_cv,value,2:(length(d2[1,])))
d2$cvs_cat <- rep("5",length=length(d2[,1]))

d3 <- as.data.frame(cbind(dataforplot1$freq,wide))
colnames(d3) <- c("freq","wide")
d3 <- d3 %>% gather(codNames_cv,value,2:(length(d3[1,])))
d3$cvs_cat <- rep("6",length=length(d3[,1]))

dataforplot2 <- dataforplot %>% select("freq","codNames_cv","value","cvs_cat")
dataforplot.combo <- rbind(dataforplot2,d1,d2,d3)

sp.whole <- ggplot() + 
  geom_line(data=dataforplot.combo[dataforplot.combo$cvs_cat %in% c(1,2,3),],
            aes(x=freq,y=value,group=codNames_cv,color=cvs_cat)) +
  geom_line(data=dataforplot.combo[dataforplot.combo$cvs_cat %in% c(4),],
            aes(x=freq,y=value),color="red",size=2) +
  geom_line(data=dataforplot.combo[dataforplot.combo$cvs_cat %in% c(5),],
            aes(x=freq,y=value),color="green",size=2) +
  geom_line(data=dataforplot.combo[dataforplot.combo$cvs_cat %in% c(6),],
            aes(x=freq,y=value),color="blue",size=2) +
  scale_y_log10() +
  xlim(0,0.7) + 
  theme_classic() +
  #ylab(expression(log[10])) + 
  ylab("") + 
  xlab("frequency") + 
  theme(legend.position = "none") +
  ggtitle("") +
  geom_text_repel(data=dataforplot[dataforplot.combo$freq==0.5,],
                  aes(label=codNames_cv,x=0.5,y=value,color=cvs_cat),
                  size = 3, xlim=c(0.5,0.7), force=10,
                  na.rm = TRUE) 

pDc <- pMG$Celtic
pDw <- pMG$W_Baltic
pDi <- pMG$Iceland

# time series
pTc 
pTw 
pTi 

# ==============
# Fig 3abc
# ==============
fig2plotlist <- list(pDc,pDw,pDi,pTc,pTw,pTi,sp.whole)
lay <- rbind(c(1,4,7,7,7),
             c(2,5,7,7,7),
             c(3,6,7,7,7))

tiff(file='C:/Users/Mikaela/Documents/GitHub/popdy/cod_figures/manuscript3/fig2abc_example_pops.tiff',units="in", width=8, height=7, res=300) 
grid.arrange(grobs = fig2plotlist,
             layout_matrix = lay)
dev.off()
rm(pDc,pDw,pDi,pTc,pTw,pTi,sp.whole,fig2plotlist)



# *******
# Fig S5: Cumulative Variance from spectra plot
# *******
codNames <- as.character(codNames)
sumL <- as.list(rep(NA,length=length(codNames)))
for(c in 1:length(codNames)){
  
  dataforplot <- specdatalong[specdatalong$codNames == codNames[c] &
                              specdatalong$alphavalue == alphas[2],]
  # vector to store cum var at each freq level
  sums <- rep(NA,length=length(dataforplot$freq))
  
  #step through each level of freq
  for(i in 1:length(sums)){ 
    #first identify which frequencies you care about
    sumoverthese <- seq(from=dataforplot$freq[1],to=dataforplot$freq[i],by=dataforplot$freq[1])
    #subset on those freq, then sum the variance values
    sums[i] <- sum(dataforplot[dataforplot$freq %in% sumoverthese,]$value)
  }
   
  #store df: codNames, freq, sum of var up through that freq, and peak age
  sumdf <- as.data.frame(cbind(dataforplot$freq,
                               sums,
                               rep(codNames[c],length=length(sums)), 
                               rep(eigentable[eigentable$codNames == codNames[c],]$mode_age, length=length(sums))))
  colnames(sumdf) <- c("freq","sums","codNames","peak")
  sumdf$freq <- as.numeric(as.character(sumdf$freq)) #fixing factors 
  sumdf$peak <- as.numeric(as.character(sumdf$peak)) #fixing factors
  sumdf$sums <- as.numeric(as.character(sumdf$sums)) #fixing factors
  sumdf$freqpeak <- (sumdf$freq)*(sumdf$peak)
  sumL[[c]] <- sumdf 
  rm(sumoverthese,sums,sumdf)
}
sumsdf <- bind_rows(sumL,id=NULL)
sumsdf$sums <- as.numeric(sumsdf$sums)
sumsdf$freq <- as.numeric(as.character(sumsdf$freq))
head(sumsdf)
str(sumsdf)

poverall <- ggplot(data=sumsdf,aes(x=freqpeak,y=sums,group=codNames,color=codNames,label=codNames)) +
  geom_line() + xlab("freq*peak") +
  theme_bw() +
  scale_x_continuous(limits=c(0,5.5),breaks=seq(from=0,to=5.5,by=0.5)) +
  #ggtitle("(c) Cumulative variance (x-axis is freq*peak)") +
  ylab("Cumulative variance") + xlab("Frequency*Peak age") +
  scale_color_discrete(name="Cod Population") #+ scale_y_log10()

  
# pzoom <- ggplot(data=sumsdf,aes(x=freqpeak,y=sums,group=codNames,color=codNames,label=codNames)) +
#   geom_line() + xlab("freq*peak") +
#   theme_bw() +
#   scale_x_continuous(limits=c(0,1.5),breaks=seq(from=0,to=1.5,by=0.5)) +
#   scale_y_continuous(limits=c(15000,30000)) +
#   ggtitle("(d) Cumulative variance, zoomed in (x-axis is freq*peak)") +
#   ylab("cumulative variance")
# 
# ggplot(dataforplot, aes(x=freq,y=value,group=codNames_peak,color=codNames_peak)) + 
#   geom_line() + theme_classic() + 
#   geom_text_repel(data=dataforplot[dataforplot$freq == 0.5,],
#                   aes(label=codNames_peak,color=codNames_peak,x=0.5,y=value),
#                   size = 3, xlim=c(0,0.7), force=10,
#                   na.rm = TRUE) +
#   ggtitle("Spectra for all populations (smooth = 43y, Pauly)") +
#   ylab("log scale") + xlab("frequency") +
#   scale_color_manual(values = getPalette(colourCount)) +
#   #scale_colour_discrete(guide = 'none')  + 
#   theme(plot.margin =margin(0,0,0,0,"cm"),legend.position = "right") +
#   xlim(0,0.7) + scale_y_log10() 

# Appendix
tiff(file='C:/Users/Mikaela/Documents/GitHub/popdy/cod_figures/manuscript3/SI/FigS5_cumulative_cd_smooth43_freqpeak.tiff', units="in", width=6, height=6.5, res=300) 
poverall
#pab <- list(poverall,pzoom)
#do.call(grid.arrange,c(pab,ncol=2))
dev.off()
rm(poverall)
# jj <- ggplot(dataforplot, aes(x=freq,y=value,group=codNames)) +
#   geom_line(aes(color=maxage)) + theme_classic() +
#   ggtitle("b") +
#   ylab("") + xlab("") +
#   theme(plot.title = element_text(size = 12)) +
#   guides(fill=guide_legend(title="Max age"))
# 
# jjj <- ggplot(dataforplot, aes(x=freq,y=value,group=codNames)) +
#   geom_line(aes(color=peak)) + theme_classic() +
#   guides(fill=guide_legend(title="CV")) +
#   ggtitle("c") +
#   ylab("") + xlab("") +
#   theme(plot.title = element_text(size = 12)) +
#   scale_y_log10()

#set all pops to one color
# line.cols <- rep("grey",16) 
# CN <- c("Coas","cod3M","cod3NO","cod3Ps","Northsea",
#         "Faroe","GB","GM","Iceland","Kat","NGulf","W_Scotland",
#         "Celtic","NE_Arctic","cod2J3KL","W_Baltic")
# CN.t <- as.data.frame(cbind(CN,line.cols),stringsAsFactors = FALSE)
# names(CN.t) <- c("pops","line.cols")
# CN.t[CN.t$pops == "NE_Arctic",]$line.cols <- "tomato"
# CN.t[CN.t$pops == "Celtic",]$line.cols <- "orange"
# CN.t[CN.t$pops == "cod2J3KL",]$line.cols <- "dodgerblue"
# CN.t[CN.t$pops == "W_Baltic",]$line.cols <- "green3"
# 
# jjjj <- ggplot(dataforplot) + 
#   geom_line(aes(x=freq,y=value,color=codNames),size=1) + theme_classic() + 
#   scale_color_manual(values=CN.t$line.cols) +
#   ggtitle("b") + 
#   ylab("Variance") + xlab("") +
#   theme(plot.title = element_text(size = 14,hjust = -0.2),legend.position = "none" ) 
# 
# # create spectra plot for a=1.1 (depleted)
# dataforplotD <- specdatalong[specdatalong$variable.type == "recruits" &
#                               specdatalong$alphavalue %in% c(1.1),] # select one alpha value
# dataforplotD$peak <- eigentable[match(dataforplotD$codNames, eigentable$codNames),"mode_age"]
# dataforplotD$maxage <- eigentable[match(dataforplotD$codNames, eigentable$codNames),"max_ages"]
# dataforplotD$sd_mode <- eigentable[match(dataforplotD$codNames, eigentable$codNames),"sd_mode"]
# dataforplotD$cvs_mode <- eigentable[match(dataforplotD$codNames, eigentable$codNames),"cvs_mode"]
# dataforplotD$codNames <- factor(dataforplotD$codNames, 
#                                levels=c("Coas","cod3M","cod3NO","cod3Ps","Northsea",
#                                         "Faroe","GB","GM","Iceland","Kat","NGulf","W_Scotland",
#                                         "Celtic","NE_Arctic","cod2J3KL","W_Baltic"))
# 
# ddddcut <- ggplot(dataforplotD) + 
#   geom_line(aes(x=freq,y=value,color=codNames),size=1) + theme_classic() + 
#   scale_color_manual(values=CN.t$line.cols) +
#   #geom_line(data=dataforplot[dataforplot$codNames=="NE_Arctic",],aes(x=freq,y=value)) +
#   ggtitle("d") +
#   ylab("Variance") + xlab("") +
#   theme(plot.title = element_text(size = 14,hjust = -0.2),legend.position = "none" ) +
#   ylim(0,2000)
# ddddylog <- ggplot(dataforplotD) + 
#   geom_line(aes(x=freq,y=value,color=codNames),size=1) + theme_classic() + 
#   scale_color_manual(values=CN.t$line.cols) +
#   #geom_line(data=dataforplot[dataforplot$codNames=="NE_Arctic",],aes(x=freq,y=value)) +
#   ggtitle("e") +
#   ylab("Variance") + xlab("") +
#   theme(plot.title = element_text(size = 14,hjust = -0.2),legend.position = "none" ) +
#   scale_y_log10()
# #p <- list(jjjj,dddd)
# p <- list(ddddcut,ddddylog)
# tiff(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/manuscript/fig4bc_spectra_log_ycutoffsigR0.3.tiff', units="in", width=4, height=7, res=300) 
# do.call(grid.arrange,c(p,ncol=1,bottom="Frequency"))
# dev.off()


# *************************************** #
# (Fig S6) Spectra plots arranged by k and peak spawning age
# *************************************** #
# Each plot is one spectra curve for a pop and k value
# Pops are arranged in order of peak spawning age - horizontal axis
# For each pop, plots for different k values are arranged vertically

# add column for k values (k=1/alpha)
specdatalong$alphavalue <- as.numeric(specdatalong$alphavalue)
specdatalong$kval <- round(1/specdatalong$alphavalue,digits=3)
head(specdatalong)
# choose only a few k values to plot/show other analyses
plotdat <- as.data.frame(specdatalong[specdatalong$alphavalue %in% selectedalphas &                                        specdatalong$variable.type == "recruits",])
# plot 1: pops in order of peak age
plotdat$codNames <- factor(plotdat$codNames, levels=codNames_ordered_by_peak)
str(plotdat)
plotdat$kval <- factor(plotdat$kval,levels=rev(levels(factor(plotdat$kval))))

tiff(file='C:/Users/Mikaela/Documents/GitHub/popdy/cod_figures/manuscript3/SI/figS6a_k_peakageorder_spectra_sigmaR0.1.tiff', units="in", width=10, height=7, res=300) 
ggplot(plotdat[plotdat$codNames == codNames[1:8],], aes(x=freq,y=value)) + 
  geom_line() + facet_grid(kval~codNames) +  scale_y_log10() +
  ggtitle("Plot 1: populations ordered by peak spawning age") #ylim(0,1000)
dev.off()
tiff(file='C:/Users/Mikaela/Documents/GitHub/popdy/cod_figures/manuscript3/SI/figS6b_k_peakageorder_spectra_sigmaR0.1.tiff', units="in", width=10, height=7, res=300) 
ggplot(plotdat[plotdat$codNames == codNames[9:16],], aes(x=freq,y=value)) + 
  geom_line() + facet_grid(kval~codNames) +  scale_y_log10() +
  ggtitle("Plot 2: populations ordered by peak spawning age") 
dev.off()
 
# plot 2: remove pops with truncated distributions (I did a visual check)
#rmpops <- c("Coas","NE_Arctic","Kat","W_Scotland","NGulf","cod3NO","cod3M")
#plotdatsub <- subset(plotdat, !(codNames %in% rmpops))
#plotdat$codNames <- factor(plotdat$codNames, levels=codNames_ordered_by_peak)
#pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/k_grid_peakageorder_truncrm.pdf', width=10, height=7) 
#ggplot(plotdatsub, aes(x=freq,y=value)) + 
#  geom_line() + facet_grid(kval~codNames) + 
#  ggtitle("Populations ordered by peak spawning age, truncated distributions removed") +
#  ylim(0,60000)
#dev.off()

# plot 3: plot pops in order of max age
#plotdat$maxage <- eigentable[match(plotdat$codNames,eigentable$codNames),"max_ages"]
#plotdat$codNames <- reorder(plotdat$codNames,plotdat$maxage)
#pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/fig4_k_maxageorder_spectra.pdf', width=10, height=7) 
#ggplot(plotdat[plotdat$codNames == codNames[1:8],], aes(x=freq,y=value)) + 
#  geom_line() + facet_grid(kval~codNames) + 
#  ggtitle("Plot 1: populations ordered by max age") +
#  ylim(0,60000)
#ggplot(plotdat[plotdat$codNames == codNames[9:16],], aes(x=freq,y=value)) + 
#  geom_line() + facet_grid(kval~codNames) + 
#  ggtitle("Plot 2: populations ordered by max age") +
#  ylim(0,60000)
#dev.off()

# *************************************** #
# (11) Area under curve
# *************************************** #
head(specdatalong)
# Plan:
# 1. calculate AUC at high and low frequencies (threshold=peak spawning age, others)
# 2. plot %AUClow and %AUChigh vs peak spawning age 
# 3. plot %AUClow and %AUChigh vs max age
# 4. plot %AUClow and %AUChigh vs spawning distribution sd & CV
# 5. make these calculations for each value of k


freq <- 0.00125
#threshold for high/low frequencies is 1/(2T) 
eigentable$AUCthreshold <- 1/(eigentable$mode_age*2)
AUCthreshold_ordered_by_peak <- eigentable %>% arrange(mode_age) %>% pull(AUCthreshold)
#OR threshold could be the same for all pops:
AUCthreshold_lowest <- rep(min(eigentable$AUCthreshold),length=length(eigentable$AUCthreshold))

AUC_amount_highL <- as.list(rep(NA,length=length(alphas)))
AUC_amount_lowL <- as.list(rep(NA,length=length(alphas)))
AUC_percent_highL <- as.list(rep(NA,length=length(alphas)))
AUC_percent_lowL <- as.list(rep(NA,length=length(alphas)))
AUC_total_L <- as.list(rep(NA,length=length(alphas)))

names(AUC_amount_highL) <- alphas
names(AUC_amount_lowL) <- alphas
names(AUC_percent_highL) <- alphas
names(AUC_percent_lowL) <- alphas
names(AUC_total_L) <- alphas

for (j in 1:length(alphas)){ #for each alpha (ie kval)...
  
  AUC_amount_high <- rep(NA,length=length(codNames_ordered_by_peak))
  AUC_percent_high <- rep(NA,length=length(codNames_ordered_by_peak))
  AUC_amount_low <- rep(NA,length=length(codNames_ordered_by_peak))
  AUC_percent_low <- rep(NA,length=length(codNames_ordered_by_peak))
  AUC_total <- rep(NA,length=length(codNames_ordered_by_peak))
  
  for (i in 1:length(codNames_ordered_by_peak)){ #step through the pops
    
    # Check inequality sign for percent high or low variability
    # multiple ht of each bar under curve by the width (one unit of x, freq)
    AUC_amount_high[i] <- sum(freq*specdatalong[specdatalong$variable.type == "recruits" 
                                  & specdatalong$codNames==codNames_ordered_by_peak[i]
                                  & specdatalong$alphavalue==alphas[j]
                                  & specdatalong$freq > AUCthreshold_ordered_by_peak[i],]$value)
                               # > is for %high, =< is for %low
    AUC_amount_low[i] <- sum(freq*specdatalong[specdatalong$variable.type == "recruits" 
                                                & specdatalong$codNames==codNames_ordered_by_peak[i]
                                                & specdatalong$alphavalue==alphas[j]
                                                & specdatalong$freq <= AUCthreshold_ordered_by_peak[i],]$value)
    
    AUC_total[i] <- sum(freq*specdatalong[specdatalong$variable.type == "recruits" 
                                          & specdatalong$codNames==codNames_ordered_by_peak[i]
                                          & specdatalong$alphavalue==alphas[j],]$value)
    AUC_percent_high[i] <- AUC_amount_high[i]/AUC_total[i]
    AUC_percent_low[i] <- AUC_amount_low[i]/AUC_total[i]
    
  }
  #store percents for each k value
  AUC_amount_highL[[j]] <- AUC_amount_high 
  AUC_amount_lowL[[j]] <- AUC_amount_low 
  AUC_percent_highL[[j]] <- AUC_percent_high 
  AUC_percent_lowL[[j]] <- AUC_percent_low
  AUC_total_L[[j]] <- AUC_total
  print(j)
}
rm(i,j)
AUC_amount_highdf <- data.frame(do.call(cbind,AUC_amount_highL))
AUC_amount_lowdf <- data.frame(do.call(cbind,AUC_amount_lowL))
AUC_percent_highdf <- data.frame(do.call(cbind,AUC_percent_highL))
AUC_percent_lowdf <- data.frame(do.call(cbind,AUC_percent_lowL))
AUC_total_df <- data.frame(do.call(cbind,AUC_total_L))

AUC_amount_highdf$codNames <- codNames_ordered_by_peak
AUC_amount_lowdf$codNames <- codNames_ordered_by_peak
AUC_percent_highdf$codNames <- codNames_ordered_by_peak
AUC_percent_lowdf$codNames <- codNames_ordered_by_peak
AUC_total_df$codNames <- codNames_ordered_by_peak

# convert dfs to long format
AUC_amount_highdflong <- AUC_amount_highdf %>% 
  gather(alpha,value,1:length(alphas)) %>% 
  separate(alpha,c("addedX","alphaval"),sep="X") %>% 
  select(-addedX) %>% 
  mutate(AUCdes=rep("amt_high"))

AUC_amount_lowdflong <- AUC_amount_lowdf %>% 
  gather(alpha,value,1:length(alphas)) %>% 
  separate(alpha,c("addedX","alphaval"),sep="X") %>% 
  select(-addedX) %>% 
  mutate(AUCdes=rep("amt_low"))

AUC_percent_highdflong <- AUC_percent_highdf %>% 
  gather(alpha,value,1:length(alphas)) %>% 
  separate(alpha,c("addedX","alphaval"),sep="X") %>% 
  select(-addedX) %>% 
  mutate(AUCdes=rep("per_high"))

AUC_percent_lowdflong <- AUC_percent_lowdf %>% 
  gather(alpha,value,1:length(alphas)) %>% 
  separate(alpha,c("addedX","alphaval"),sep="X") %>% 
  select(-addedX) %>% 
  mutate(AUCdes=rep("per_low"))

AUC_total_dflong <- AUC_total_df %>% 
  gather(alpha,value,1:length(alphas)) %>% 
  separate(alpha,c("addedX","alphaval"),sep="X") %>% 
  select(-addedX) %>% 
  mutate(AUCdes=rep("total"))

AUCdat <- rbind(AUC_amount_highdflong,
                AUC_amount_lowdflong,
                AUC_percent_highdflong,
                AUC_percent_lowdflong,
                AUC_total_dflong  )

# add columns
#eigentable$peakovermax <- round(eigentable$mode_age / eigentable$max_ages,digits=2)
#AUCdat$peakovermax <- eigentable[match(AUCdat$codNames,eigentable$codNames),"peakovermax"]
AUCdat$peak <- eigentable[match(AUCdat$codNames,eigentable$codNames),"mode_age"]
AUCdat$maxage <- eigentable[match(AUCdat$codNames,eigentable$codNames),"max_ages"]
AUCdat$codNames_plot <- eigentable[match(AUCdat$codNames,eigentable$codNames),"codNames_plot"]
AUCdat$cvs <- eigentable[match(AUCdat$codNames,eigentable$codNames),"cvs_mode"]

AUCdat$cvs <- round(AUCdat$cvs,digits=2)
AUCdat$alphaval <- as.numeric(AUCdat$alphaval)
AUCdat$kval <- round(1/(AUCdat$alphaval*conLEP^2),digits = 2)
AUCdat$codNames_plot_no <- paste(AUCdat$codNames_plot,"(",AUCdat$cvs,")",sep="") #new col w/codName+cv
AUCdat$codNames_plot_no_peak <- paste(AUCdat$codNames_plot,"(",AUCdat$peak,")",sep="") #new col w/peak
AUCdat$codNames_plot_no_maxage <- paste(AUCdat$codNames_plot,"(",AUCdat$maxage,")",sep="")
#AUCdat$codNames_plot_no_peakovermax <- paste(AUCdat$codNames_plot,"(",AUCdat$peakovermax," p/m)",sep="")


# set factor levels
AUCdat$kval <- factor(AUCdat$kval,levels=unique(AUCdat$kval))

# set codNames in order of CV
codNames_plot_cvs_order <- unique(AUCdat$codNames_plot_no)
codNames_plot_cvs_order <- factor(codNames_plot_cvs_order,
                                  levels=unique(AUCdat[order(AUCdat$cvs),]$codNames_plot_no))
AUCdat$codNames_plot_no <- factor(AUCdat$codNames_plot_no,levels=levels(codNames_plot_cvs_order))

# set codNames in order of peak age
codNames_plot_peak_order <- unique(AUCdat$codNames_plot_no_peak)
codNames_plot_peak_order <- factor(codNames_plot_peak_order,
                                   levels=unique(AUCdat[order(AUCdat$peak),]$codNames_plot_no_peak))
AUCdat$codNames_plot_no_peak <- factor(AUCdat$codNames_plot_no_peak,
                                       levels=levels(codNames_plot_peak_order))
# set codNames in order of max age
codNames_plot_max_order <- unique(AUCdat$codNames_plot_no_maxage)
codNames_plot_max_order <- factor(codNames_plot_max_order,
                                  levels=unique(AUCdat[order(AUCdat$maxage),]$codNames_plot_no_maxage))
AUCdat$codNames_plot_no_maxage <- factor(AUCdat$codNames_plot_no_maxage,
                                         levels=levels(codNames_plot_max_order))
# set codNames in order of peakovermax
#codNames_plot_peakovermax_order <- unique(AUCdat$codNames_plot_no_peakovermax)
#codNames_plot_peakovermax_order <- factor(codNames_plot_peakovermax_order,
#                                  levels=unique(AUCdat[order(AUCdat$peakovermax),]$codNames_plot_no_peakovermax))
#AUCdat$codNames_plot_no_peakovermax <- factor(AUCdat$codNames_plot_no_peakovermax,
#                                         levels=levels(codNames_plot_peakovermax_order))

# *************************************** #
# Gantt plot for fraction high frequency variance
# *************************************** #
# (a) lambda2/lambda1
# start <- rep(NA,length=length(selectedkvals))
# end <- rep(NA,length=length(selectedkvals))
# for(i in 1:length(selectedkvals)){
#   start[i] <- min(eigendata[eigendata$kvals == selectedkvals[i] & eigendata$eigen == "e12",]$value)
#   end[i] <- max(eigendata[eigendata$kvals == selectedkvals[i] & eigendata$eigen == "e12",]$value)}
# df <- as.data.frame(cbind(selectedkvals,start,end))
# colnames(df) <- c("k","start","end")
# df$k <- as.factor(df$k)
# df$diff <- as.factor(round(df$end-df$start,digits=2))
# gantt1 <- ggplot(data=df,aes(x=start,xend=end,y=k,yend=k,label=diff)) + 
#   theme_bw() +
#   geom_segment(size=8,color="grey") +
#   labs(x=expression(paste(abs(lambda[2])/lambda[1])),y="k") +
#   xlim(0.86,1) +
#   geom_label(x=0.99,size=4)
# rm(start,end,i,df)
# 
# # (b) fraction high frequency variance
# start <- rep(NA,length=length(selectedkvals))
# end <- rep(NA,length=length(selectedkvals))
# for(i in 1:length(selectedkvals)){
#   start[i] <- min(AUCdat[AUCdat$kval == selectedkvals[i] & AUCdat$AUCdes == "percent",]$value)
#   end[i] <- max(AUCdat[AUCdat$kval == selectedkvals[i] & AUCdat$AUCdes == "percent",]$value)}
# df <- as.data.frame(cbind(selectedkvals,start,end))
# colnames(df) <- c("k","start","end")
# df$k <- as.factor(df$k)
# df$diff <- as.factor(round(df$end-df$start,digits=2))
# gantt2 <- ggplot(data=df,aes(x=start,xend=end,y=k,yend=k,label=diff)) + 
#   theme_bw() +
#   geom_segment(size=8,color="grey") +
#   labs(x="Fraction of high frequency variance",y="k") +
#   xlim(0,0.6) +
#   geom_label(x = 0.55, size=4)
# tiff(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/manuscript/figXX_gantt_highfreq_threholdpeak.tiff', units="in", width=4, height=3.5, res=300)
# g <- list(gantt1,gantt2)
# do.call(grid.arrange,c(g,ncol=2))
# dev.off()
# rm(start,end,i,df,gantt1,gantt2,g)


# *************************************** #
# (Fig 4a) Amount of total variance
# (Fig 4b) Fraction of high frequency variance

# *************************************** #
head(ts.data)
ts.data$alphavalue <- as.numeric(as.character(ts.data$alphavalue))
#plotalpha <- unique(ts.data$alphavalue)[c(1,10,20,30,40,50,60,70,80,90)]
plotalpha <- c(8.26,7.51,6.36,5.51,4.59,3.31,2.07,1.38,0.97)
plotalpha <- c(5.51,3.31,2.07,1.38,0.97)
alphas <- plotalpha
#plotalpha <- round(1/(kvals*1.1^2),digits=2)
round(1/(alphas*conLEP^2),digits=2)

# ---
# Calculate Total Variance in each time series:
# ---
varL = as.list(rep(NA,length=length(plotalpha))) 
for(i in 1:length(codNames)){ # for each pop i
  
  #subset ts data to pop i only
  dat <- ts.data[ts.data$codNames == codNames[i] & 
                 ts.data$variable == "recruits" & 
                 ts.data$alphavalue %in% plotalpha,]
  #store variance in ts data at different alpha levels
  variance <- rep(NA,length=length(plotalpha))
  #store mean of ts data at different alpha levels
  means <- rep(NA,length=length(plotalpha))
  
  for(b in 1:length(plotalpha)){ #step through alpha values
    
    # calculate equilibrium value (mean) for each time series (rm first 200 ts)
    means[b] <- mean(dat[dat$alphavalue == plotalpha[b],]$value[rm_first_timesteps:(timesteps-2)])
    # substract mean from time series (vals_meanrm)
    vals = dat[dat$alphavalue == plotalpha[b],]$value[rm_first_timesteps:(timesteps-2)] #- means[b]
    # sq root the variance & then divide by mean 
    variance[b] <- sqrt(var(vals))/means[b]
  
  }
  varL[[i]] <- as.data.frame(cbind(variance,plotalpha,rep(as.character(codNames[i]),length=length(plotalpha))))
  names(varL[[i]]) <- c("variance","alphavalue","codNames")
}
vardat <- do.call(rbind,varL)
rm(i,b,dat,variance,means,vals) #clean up
head(vardat)
vardat$peak <- eigentable[match(vardat$codNames,eigentable$codNames),"mode_age"]
vardat$kval <- round((1/(as.numeric(as.character(vardat$alphavalue))*conLEP^2)),digits = 2)
vardat$codNames_plot <- eigentable[match(vardat$codNames,eigentable$codNames),"codNames_plot"]
vardat$cvs <- eigentable[match(vardat$codNames,eigentable$codNames),"cvs_mode"]
vardat$maxage <- eigentable[match(vardat$codNames,eigentable$codNames),"max_ages"]
vardat$codNames_plot_no_maxage <- paste(vardat$codNames_plot,"(",vardat$maxage,")",sep="")
vardat$sd <- eigentable[match(vardat$codNames,eigendata$codNames),"sd_mode"]
#vardat$peakovermax <- round(vardat$peak/vardat$maxage,digits=2)
#vardat$codNames_plot_no_peakovermax <- paste(vardat$codNames_plot,"(",vardat$peakovermax," p/m)",sep="")

# set k slopes to factor, order by increasing slope
vardat$kval <- factor(vardat$kval,levels=unique(vardat$kval))

# make sure variance values are numeric
vardat$variance <- as.numeric(as.character(vardat$variance))

# set codNames factor levels to increase with peak age OR CV
#vardat$codNames_plot <- factor(vardat$codNames_plot,levels=codNames_ordered_by_peak_plot)
codNames_ordered_by_cvs_plot <- eigentable %>% arrange(cvs_mode,codNames) %>% pull(codNames_plot)
codNames_ordered_by_cvs <- eigentable %>% arrange(cvs_mode,codNames) %>% pull(codNames)

# set codNames in order of peak spawning age
vardat$codNames <- factor(vardat$codNames,levels=codNames_ordered_by_cvs)
vardat$codNames_plot <- factor(vardat$codNames_plot,levels=codNames_ordered_by_cvs_plot)

# create codNames col ordered by CV
vardat$cvs <- round(vardat$cvs,digits=2) #round off cv values
vardat$codNames_plot_no <- paste(vardat$codNames_plot,"(",vardat$cvs,")",sep="") #new col w/codName+cv
codNames_plot_no_order <- unique(vardat$codNames_plot_no) #creating factor levels for new names
codNames_plot_no_order <- factor(codNames_plot_no_order,levels=unique(vardat[order(vardat$cvs),]$codNames_plot_no),ordered=TRUE)
vardat$codNames_plot_no <- factor(vardat$codNames_plot_no,levels=levels(codNames_plot_no_order))
# create codNames col ordered by peak
vardat$codNames_plot_no_peak <- paste(vardat$codNames_plot,"(",vardat$peak,")",sep="") #new col w/codName+peak
codNames_plot_no_peak_order <- unique(vardat$codNames_plot_no_peak)
codNames_plot_no_peak_order <- factor(codNames_plot_no_peak_order,levels=unique(vardat[order(vardat$peak),]$codNames_plot_no_peak))
vardat$codNames_plot_no_peak <- factor(vardat$codNames_plot_no_peak,levels=levels(codNames_plot_no_peak_order))
# set codNames in order of max age
codNames_plot_max_order <- unique(vardat$codNames_plot_no_maxage)
codNames_plot_max_order <- factor(codNames_plot_max_order,
                                  levels=unique(vardat[order(vardat$maxage),]$codNames_plot_no_maxage))
vardat$codNames_plot_no_maxage <- factor(vardat$codNames_plot_no_maxage,
                                         levels=levels(codNames_plot_max_order))

# get range of total variance at least depleted (k=0.1)
max.var.k0.15 <- max(vardat[vardat$kval == 0.15,]$variance)
min.var.k0.15 <- min(vardat[vardat$kval == 0.15,]$variance)
# get range of total variance at most depleted (k=0.9)
max.var.k0.85 <- max(vardat[vardat$kval == 0.85,]$variance)
min.var.k0.85 <- min(vardat[vardat$kval == 0.85,]$variance)
head(vardat)

# # ************************
# # Average pops with same peak age in vardat & AUCdat
# # ************************
# head(vardat)
# ks <- round(1/alphas,digits=1)
# peaks <- sort(unique(eigentable$mode_age),decreasing=FALSE) #peak ages of pops
# codNames_peakavg <- c("peak_avg_3","peak_avg_4","peak_avg_5","Coas_6","peak_avg_8","NE_Arctic_9")
# var_meansL <- as.list(rep(NA,length=length(ks)))
# names(var_meansL) <- ks
# for(i in 1:length(ks)){
#   var_means <- rep(NA,length=length(peaks))
#   for(j in 1:length(peaks)){
#     var_means[j] <- mean(vardat[vardat$kval==ks[i] & vardat$peak==peaks[j],]$variance)
#   }
#   kval <- rep(ks[i],length=length(peaks))
#   store <- as.data.frame(cbind(codNames_peakavg,peaks,var_means,kval))
#   var_meansL[[i]]<-store
#   rm(kval,store,var_means)
# }
# var_meansdf <- bind_rows(var_meansL,id=NULL)
# var_meansdf$var_means <- as.numeric(var_meansdf$var_means)
# var_meansdf$codNames_peakavg <- factor(var_meansdf$codNames_peakavg,
#                                          levels=c("peak_avg_3","peak_avg_4","peak_avg_5",
#                                                   "Coas_6","peak_avg_8","NE_Arctic_9"))
# fig5a_avg <- ggplot(var_meansdf,aes(x=kval,y=codNames_peakavg)) +
#   geom_raster(aes(fill=var_means)) + 
#   xlab("Slope on S-R curve (k)") + ylab("") +
#   scale_fill_gradient(low="purple", high="orange") + 
#   theme_classic() +
#   guides(fill=guide_legend(title="Standard deviation\n normalized to the mean", reverse=FALSE)) +
#   theme(axis.text.x = element_text(angle = 70, hjust = 1),
#         legend.position = "top",
#         legend.title = element_text(size = 8), 
#         legend.text = element_text(size = 8)) +
#   ylab("ordered by peak")
# # average populations in AUCdat
# head(AUCdat)
# ks <- round(1/alphas,digits=1)
# peaks <- sort(unique(eigentable$mode_age),decreasing=FALSE) #peak ages of pops
# codNames_peakavg <- c("peak_avg_3","peak_avg_4","peak_avg_5","Coas_6","peak_avg_8","NE_Arctic_9")
# auc_meansL <- as.list(rep(NA,length=length(ks)))
# names(auc_meansL) <- ks
# for(i in 1:length(ks)){
#   auc_means <- rep(NA,length=length(peaks))
#   for(j in 1:length(peaks)){
#     auc_means[j] <- mean(AUCdat[AUCdat$kval==ks[i] & AUCdat$peak==peaks[j] & AUCdat$AUCdes=="percent",]$value)
#   }
#   kval <- rep(ks[i],length=length(peaks))
#   store <- as.data.frame(cbind(codNames_peakavg,peaks,auc_means,kval))
#   auc_meansL[[i]]<-store
#   rm(kval,store,auc_means)
# }
# auc_meansdf <- bind_rows(auc_meansL,id=NULL)
# auc_meansdf$auc_means <- as.numeric(auc_meansdf$auc_means)
# auc_meansdf$codNames_peakavg <- factor(auc_meansdf$codNames_peakavg,
#                                        levels=c("peak_avg_3","peak_avg_4","peak_avg_5",
#                                                 "Coas_6","peak_avg_8","NE_Arctic_9"))
# 
# fig5b_avg <- ggplot(auc_meansdf,aes(x=kval,y=codNames_peakavg)) +
#   geom_raster(aes(fill=auc_means)) + 
#   xlab("Slope on S-R curve (k)") + ylab("") +
#   scale_fill_gradient(low="purple", high="orange") + 
#   theme_classic() +
#   guides(fill=guide_legend(title="High frequency variance\n(fraction range 0.06-0.5)", reverse=FALSE)) +
#   theme(axis.text.x = element_text(angle = 70, hjust = 1),
#         legend.position = "top",
#         legend.title = element_text(size = 8), 
#         legend.text = element_text(size = 8)) +
#   ylab("ordered by peak")
# rm(auc_means,i,j,auc_meansL)
# 
# # bar chart for averaged pops
# auc_meansdf
# diffs <- rep(NA,length=length(codNames_peakavg))
# for(i in 1:length(unique(auc_meansdf$codNames_peakavg))){
#   diffs[i] <- auc_meansdf[auc_meansdf$codNames_peakavg == unique(auc_meansdf$codNames_peakavg)[i] &
#   auc_meansdf$kval == 0.1,]$auc_means - auc_meansdf[auc_meansdf$codNames_peakavg == unique(auc_meansdf$codNames_peakavg)[i] & auc_meansdf$kval == 0.9,]$auc_means
#   
# }
# diffsdf <- as.data.frame(cbind(codNames_peakavg,peaks,diffs))
# diffsdf$codNames_peakavg <- factor(diffsdf$codNames_peakavg,
#                                        levels=c("peak_avg_3","peak_avg_4","peak_avg_5",
#                                                 "Coas_6","peak_avg_8","NE_Arctic_9"))
# diffsdf$diffs <- as.numeric(as.character(diffsdf$diffs))
# str(diffsdf)
# 
# fig5c_avg <- ggplot(data=diffsdf, aes(x=codNames_peakavg, y=diffs)) + geom_bar(stat="identity") + 
#   theme_classic() +
#   ggtitle("")  +
#   theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
#   coord_flip() +
#   ylab("Change in the fraction of high freq\n variance as k goes from 0.1-0.9 ") + 
#   xlab("") 
# 
# tiff(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/manuscript/fig5abc_peak_avg_2T_sigR0.3_span1.5.tiff', units="in", width=11, height=6, res=300) 
# p <- list(fig5a_avg,fig5b_avg,fig5c_avg)
# do.call(grid.arrange,c(p,ncol=3))
# dev.off()
# rm(fig5a_avg,fig5b_avg,fig5c_avg)
# 
# 
# # ************************
# # Plot total var at different k vals
# # ************************


fig4bFractionHIGH <- ggplot(data=AUCdat[AUCdat$alphaval %in% alphas & AUCdat$AUCdes=="per_high",],
                            aes(x=kval,y=codNames_plot_no)) + 
  geom_raster(aes(fill=value)) + 
  xlab("Slope on egg-recruit curve at equilibrium (k)") + ylab("") + 
  scale_fill_gradient(low="purple", high="orange") + 
  scale_colour_gradient(limits = c(0, 1)) +
  theme_classic() +
  guides(fill=guide_legend(title=paste("Fraction of variance\nat high frequencies"))) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1),
        legend.position="top",
        legend.justification = c(0, 1),
        plot.title = element_text(hjust = -0.5, vjust=-0.1),
        legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8)) + theme(plot.margin=unit(c(0.1,1.5,0.1,0.1),"cm"))

# for the population with the smallest CV, what is the range of fraction high freq var?
min(AUCdat[AUCdat$cvs==min(AUCdat[AUCdat$AUCdes=="per_high",]$cvs),]$value)
max(AUCdat[AUCdat$cvs==min(AUCdat[AUCdat$AUCdes=="per_high",]$cvs),]$value)

# for the population with the largest CV,, what is the range of fraction high freq var?
min(AUCdat[AUCdat$cvs==max(AUCdat[AUCdat$AUCdes=="per_high",]$cvs),]$value)
max(AUCdat[AUCdat$cvs==max(AUCdat[AUCdat$AUCdes=="per_high",]$cvs),]$value)

fig4aTotalVar <- ggplot(data=vardat[vardat$alphaval%in% alphas,],
                   aes(x=kval,y=codNames_plot_no_maxage)) + 
  geom_raster(aes(fill=variance)) + 
  xlab("Slope on egg-recruit curve at equilibrium (k)") + ylab("") + 
  scale_fill_gradient(low="purple", high="orange") + 
  scale_colour_gradient(limits = c(0, 1)) +
  theme_classic() +
  guides(fill=guide_legend(title=paste("Total Variance"))) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1),
        legend.position="top",
        legend.justification = c(0, 1),
        legend.box="horizontal",
        plot.title = element_text(hjust = -0.8, vjust=-0.1),
        legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8)) + theme(plot.margin=unit(c(0.1,1.5,0.1,0.1),"cm"))

# when k=0.15, what is the range of total variance?
min(vardat[vardat$kval==0.15,]$variance)
max(vardat[vardat$kval==0.15,]$variance)
# when k=0.85, what is the range of total variance?
min(vardat[vardat$kval==0.85,]$variance)
max(vardat[vardat$kval==0.85,]$variance)

tiff(file="C:/Users/Mikaela/Documents/GitHub/popdy/cod_figures/manuscript3/fig4ab_totalvar_fracHigh_v1.tiff",
     units="in", width=9, height=6, res=300)
     
grid.newpage()
#grid.draw(cbind(ggplotGrob(fig5bHIGH), ggplotGrob(fig5bLOW), size = "first"))
grid.draw(cbind(ggplotGrob(fig4aTotalVar), ggplotGrob(fig4bFractionHIGH)))
#f <- list(fig5aLOW,fig5bHIGH)
#do.call(grid.arrange,c(f,ncol=2))
dev.off()

# Testing new figure: 4 panel plot
# panels = k values
# xaxis = max age
# yaxis = CV in recruitment time series
head(vardat)
newfigtest <- ggplot(vardat[vardat$kval %in% selectedkvals,],
                  aes(x=sd,y=variance)) +
  geom_point() + 
  geom_smooth(method="lm",se=FALSE,color="black") +
  facet_grid(. ~ kval) +
  #scale_y_continuous(limits=c(0.6,1)) +
  geom_text_repel(data=vardat[vardat$kval %in% selectedkvals,],
                  aes(label = codNames_plot),
                  segment.color = "grey",
                  size = 2,
                  na.rm = TRUE) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(angle = 90),
        axis.text.x = element_text(angle = 40, hjust = 1)) +
  ylab("Total variance in recruitment time series") +
  xlab("Spawning distrition Stdev") 
tiff(file="C:/Users/Mikaela/Documents/GitHub/popdy/cod_figures/manuscript3/recruitCV_vs_maxage_byKvals.tiff",
     units="in", width=9, height=6, res=300)
newfigtest
dev.off()

# **********************************
# Figure 5 -- don't include in manuscript
# Fig5a: lambda1 v total variance 
# Fig5b: lambda2/1 v fraction high freq var
# **********************************

# Plot 1
# Merge vardat and eigendata dfs based on codNames & kval
#subset eigendata for lambda1
eigendata1 <- eigendata %>% 
  filter(kvals %in% selectedkvals & eigen=="e1") %>% 
  select(kvals,codNames,value,eigen) 
colnames(eigendata1)[which(names(eigendata1) == "kvals")] <- "kval"
colnames(eigendata1)[which(names(eigendata1) == "value")] <- "value_e1"
eigendata1$kval <- as.factor(eigendata1$kval)

# subset vardat
vardat1 <- vardat %>%
  filter(kval %in% selectedkvals) %>%
  select(codNames,kval,variance,codNames_plot)
#rm(eigendata1)

# join dfs
totvar.e1.df <- right_join(vardat1,eigendata1,by=c("codNames","kval"))

fig5a <- ggplot(data=totvar.e1.df,aes(x=value_e1,y=variance)) +
  geom_point() + 
  geom_smooth(method="lm",se=FALSE,color="black") +
  facet_grid(. ~ kval) +
  geom_text_repel(data=totvar.e1.df,
                  aes(label = codNames_plot),
                  segment.color = "grey",
                  size = 2,
                  na.rm = TRUE) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) + #,axis.title.y = element_text(angle = 0)) +
  #ylab("Amount of low frequency variance") +
  #ylab(expression(paste("Rerecruitment variance ",frac(sigma,"m")))) +
  ylab(expression(paste("Rerecruitment variance ","(",sigma,"/",mu,")"))) +
  xlab(expression(paste(lambda[1]))) 


# fig5b: fraction of high frequency variance vs 1/DR
# Merge AUCdat and eigendata dfs based on codNames & kval
eigendata1 <- eigendata %>% 
  filter(kvals %in% selectedkvals & eigen=="e12") %>% 
  select(kvals,codNames,value,eigen)
colnames(eigendata1)[which(names(eigendata1) == "kvals")] <- "kval"
colnames(eigendata1)[which(names(eigendata1) == "value")] <- "value_e12"
eigendata1$kval <- as.factor(eigendata1$kval)
# AUCdat_plot_high1 <- AUCdat_plot_high %>% 
#   filter(kval %in% selectedkvals & AUCdes=="amt_high")
# AUCdat_plot_high1$kval <- as.numeric(as.character(AUCdat_plot_high1$kval))
# eigendata1$kval <- as.numeric(as.character(eigendata1$kval))
# AUC_e12 <- left_join(AUCdat_plot_high1,eigendata1,all.x=TRUE)
# rm(eigendata1,AUCdat_plot_high1)
# head(AUC_e12)

AUCdat_per_high <- AUCdat %>% 
  filter(kval %in% selectedkvals & AUCdes=="per_high")
AUCdat_per_high$kval <- as.numeric(as.character(AUCdat_per_high$kval))
eigendata1$kval <- as.numeric(as.character(eigendata1$kval))
AUC_e12 <- left_join(AUCdat_per_high,eigendata1,all.x=TRUE)
rm(eigendata1,AUCdat_per_high)
head(AUC_e12)
head(AUCdat)


fig5b <- ggplot(data=AUC_e12,aes(x=value_e12,y=value)) +
  geom_point() + 
  geom_smooth(method="lm",se=FALSE,color="black") +
  facet_grid(. ~ kval) +
  #scale_y_continuous(limits=c(0,1.1)) +
  geom_text_repel(data=AUC_e12,
                  aes(label = codNames_plot),
                  segment.color = "grey",
                  size = 2,
                  na.rm = TRUE) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(angle = 90)) +
  ylab("Fraction of high frequency variance") +
  xlab(expression(paste(frac(abs(lambda[2]),lambda[1])))) 

#e12vDR.rsq <- rep(NA,length=length(selectedkvals))
# for(n in 1:length(selectedkvals)){
#   dd <- AUC_e12[AUC_e12$kval == selectedkvals[n],]
#   e12vDR.rsq[n] <- summary(lm(value ~ value_e12,dd))$r.squared
# }
# rm(n,dd)
#ylab(expression(paste(abs(lambda[2])/lambda[1]))) +
#ylim(0.5,1.01)
tiff(file='C:/Users/Mikaela/Documents/GitHub/popdy/cod_figures/manuscript3/fig5ab_eigens_vs_simoutput_sigR0.3_span1.5_oneoverpeakhighfreq_V3.tiff', units="in", width=6, height=7, res=300) 
#p12 <- list(p1,p2)
#do.call(grid.arrange,c(p12,ncol=1))
grid.newpage()
grid.draw(rbind(ggplotGrob(fig5a), ggplotGrob(fig5b), size = "last"))
dev.off()
rm(fig5a,fig5b)



# ---------------- Extra code ------------------------ #
# ---
# Redo Fig 5a and 5b 
dataforplot <- AUCdat[AUCdat$AUCdes == "fraction",]
fig5HIGH <- ggplot(data=dataforplot,aes(x=kval,y=codNames_plot_no)) + 
  geom_raster(aes(fill=value)) + 
  xlab("Slope on egg-recruit curve at equilibrium (k)") + ylab("") + 
  scale_fill_gradient(low="purple", high="orange") + 
  scale_colour_gradient(limits = c(0, 1)) +
  theme_classic() +
  guides(fill=guide_legend(title=paste("High frequency variance\n(fraction range 0.06-0.5)"))) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1),
        legend.position="top",
        plot.title = element_text(hjust = -0.5, vjust=-0.1),
        legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8)) 


# calculate difference in fraction variance when 
# k=0.9 (alpha=1.1) vs k=0.1 (alpha=10)
pervar_k0.1 <- AUCdat[AUCdat$AUCdes == "percent" & AUCdat$alphaval == 10,]
pervar_k0.91 <- AUCdat[AUCdat$AUCdes == "percent" & AUCdat$alphaval == 1.1,]
diff <- rep(NA,length=length(codNames_ordered_by_peak))

for(d in 1:length(codNames_ordered_by_peak)){
    diff[d] <- round(AUCdat[AUCdat$AUCdes == "percent" & AUCdat$alphaval == 10 & AUCdat$codNames==codNames_ordered_by_peak[d],]$value - AUCdat[AUCdat$AUCdes == "percent" & AUCdat$alphaval == 1.1 & AUCdat$codNames==codNames_ordered_by_peak[d],]$value, digits = 2)
}
diffdata <- as.data.frame(cbind(rev(codNames_ordered_by_peak),rev(diff)))
names(diffdata) <- c("codNames","diff")
diffdata$diff <- as.numeric(as.character(diffdata$diff))

# add columns
diffdata$cvs <- round(eigentable[match(diffdata$codNames,eigentable$codNames),"cvs_mode"],digits=2)
diffdata$peak <- eigentable[match(diffdata$codNames,eigentable$codNames),"mode_age"]
diffdata$codNames_plot <- eigentable[match(diffdata$codNames,eigentable$codNames),"codNames_plot"]
diffdata$maxage <- eigentable[match(diffdata$codNames,eigentable$codNames),"max_ages"]
#diffdata$peakovermax <- round(diffdata$peak/diffdata$maxage, digits=2)
diffdata$codNames_peak <- paste(diffdata$codNames_plot,"(",diffdata$peak,")",sep="")
diffdata$codNames_cv <- paste(diffdata$codNames_plot,"(",diffdata$cv,")",sep="")  
diffdata$codNames_max <- paste(diffdata$codNames_plot,"(",diffdata$maxage,")",sep="")
#diffdata$codNames_peakovermax <- paste(diffdata$codNames_plot,"(",diffdata$peakovermax," p/m)",sep="")

# set factor levels
diffdata$codNames_peak <- factor(diffdata$codNames_peak,levels=levels(codNames_plot_no_peak_order))
diffdata$codNames_cv <- factor(diffdata$codNames_cv,levels=levels(codNames_plot_cvs_order))
diffdata$codNames_max <- factor(diffdata$codNames_max,levels=levels(codNames_plot_max_order))
#diffdata$codNames_peakovermax <- factor(diffdata$codNames_peakovermax,levels=levels(codNames_plot_peakovermax_order))

# plot diffdata barplot
fig5c <- ggplot(data=diffdata, aes(x=codNames_cv, y=diff)) + geom_bar(stat="identity") + 
  theme_classic() +
  ggtitle("")  +
  theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
  coord_flip() +
  ylab("Difference in fraction of\n high frequency variance") + 
  xlab("") 

tiff(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/manuscript2/fig7abc_highvar_CVPeak_sigR0.3_span1.5.tiff', units="in", width=9, height=6, res=300) 
f <- list(fig5a,fig5b)
do.call(grid.arrange,c(f,ncol=2))
dev.off()
rm(fig5a,fig5b,fig5c)
tiff(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/manuscript/fig5c_bars_CVorder_sigR0.3.tiff', units="in", width=4, height=6, res=300) 
fig5c
dev.off()

# **************************************
# --- Average the differences --- #
# **************************************
# try averaging the difference in % populations with the same peak age
diffdata$peak <- eigentable[match(diffdata$codNames_ordered_by_peak,eigentable$codNames),"mode_age"]
plot(diffdata %>% group_by(peak) %>% summarise(avg = mean(diff)))
diff.means <- diffdata %>% group_by(peak) %>% summarise(avg = mean(diff))
# CI for peak ages with multiple populations:
#peak_need_avg <- unique(diffdata$peak[duplicated(diffdata$peak)])
peak <- unique(diffdata$peak)
#error <- rep(NA,length=length(peak_need_avg))
upper <- rep(NA,length=length(peak_need_avg))
lower <- rep(NA,length=length(peak_need_avg))
m <- rep(NA,length=length(peak_need_avg))
for (p in 1:length(peak_need_avg)){
  n = length(diffdata[diffdata$peak == peak_need_avg[p],]$diff)
  m[p] = mean(diffdata[diffdata$peak == peak_need_avg[p],]$diff)
  sd = sd(diffdata[diffdata$peak == peak_need_avg[p],]$diff)
  #error = qt(0.975,df=n-1)*sd/sqrt(n)
  upper[p] = m[p]+(qt(0.975,df=n-1)*sd/sqrt(n))
  lower[p] = m[p]-(qt(0.975,df=n-1)*sd/sqrt(n))
}
upper[is.nan(upper)] <- 0
lower[is.nan(lower)] <- 0
diff.means <- as.data.frame(cbind(peak_need_avg,m,upper,lower))

my_breaks = seq(3,9,by=1)

fig5c.avg <- ggplot(data=diff.means, aes(x=peak, y=m)) + 
  geom_bar(stat="identity") + 
  geom_errorbar(data=diff.means, mapping=aes(x=peak, ymin=lower, ymax=upper), width=0.2) +
  scale_y_continuous(limits = c(0,0.18)) +
  scale_x_continuous(breaks=my_breaks,labels=my_breaks) +
  theme_classic() +
  ggtitle("")  +
  theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
  coord_flip() +
  ylab("Difference in percent of\n low frequency variance") + xlab("")






# ------------------------------------------------------------------------------------------------------

#Instead I could use total AUC to examine variance in time series
dataforplot <- AUCdat[AUCdat$AUCdes == "total" &
                        AUCdat$alphaval %in% plotalpha,]
ggplot(data=dataforplot,aes(x=kval,y=codNames)) +
  geom_raster(aes(fill=value)) + xlab("Slope of stock-recruit curve at equilibrium (k)") +
  ylab("") + scale_fill_gradient(low="purple", high="orange") + 
  guides(fill=guide_legend(title="Total AUC")) +
  ggtitle("Total variance (area under curve at all frequencies)") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))



# order pops by max age, facet by truncated distribution
AUCdat$codNames <- factor(AUCdat$codNames, 
                          levels=unique(AUCdat$codNames[order(AUCdat$maxage)]))
dataforplot <- AUCdat[AUCdat$AUCdes == "perhigh" &
                        AUCdat$kval %in% c(0.1,0.2,0.3,0.4,0.5,0.91),]
ht<- ggplot(data=dataforplot,aes(x=codNames,y=kval)) +
  geom_raster(aes(fill=value)) + xlab("") + 
  ylab("") + scale_fill_gradient(low="green", high="red") + 
  scale_colour_gradient(limits = c(0, 1)) +
  guides(fill=guide_legend(title="Percent at high\n frequency")) +
  ggtitle("(g) Percent of AUC at high frequencies (>0.1), percents range 0.02-0.29") +
  facet_grid(.~truncdis) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
  
dataforplot <- AUCdat[AUCdat$AUCdes == "perlow" &
                        AUCdat$kval %in% c(0.1,0.2,0.3,0.4,0.5,0.91),]
lt<- ggplot(data=dataforplot,aes(x=codNames,y=kval)) +
  geom_raster(aes(fill=value)) + xlab("") +
  ylab("") + scale_fill_gradient(low="green", high="red") + 
  scale_colour_gradient(limits = c(0, 1)) +
  guides(fill=guide_legend(title="Percent at low\n frequency")) +
  ggtitle("(h) Percent of AUC at low frequencies (<0.1), percents range 0.71-0.97") +
  facet_grid(.~truncdis) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
dataforplot <- AUCdat[AUCdat$AUCdes == "total" &
                        AUCdat$kval %in% c(0.1,0.2,0.3,0.4,0.5,0.91),]
tt<- ggplot(data=dataforplot,aes(x=codNames,y=kval)) +
  geom_raster(aes(fill=value)) + xlab("") +
  ylab("") + scale_fill_gradient(low="green", high="red") + 
  #scale_colour_gradient(limits = c(0, 1)) +
  guides(fill=guide_legend(title="Total AUC")) +
  ggtitle("(i) Total variance (area under curve at all frequencies)") +
  facet_grid(.~truncdis) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
ordered_max_trunc <- list(ht,lt,tt)

tiff(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/3_gridplot_bymaxage_truncsplit.tiff', 
     units="in", width=7, height=10, res=300) 
do.call(grid.arrange,c(ordered_max_trunc,ncol=1,left="Slope (k)",bottom="Cod Populations (increasing maximum age) split by truncated and not truncated distributions"))
dev.off()


# Extra code --------------------------------------------------------------
# *************************************** #
# (x) Area under curve - remove pops with truncated distributions
# *************************************** #
keep <- c("Northsea","W_Baltic","Faroe","Celtic","Iceland","GB","GM",
              "cod2J3KL","cod3Ps")
codNames_ordered_by_peak_keep <- codNames_ordered_by_peak[codNames_ordered_by_peak %in% keep]
eigentable.trunc <- eigentable[eigentable$codNames %in% codNames_ordered_by_peak_keep,]

AUC_greater <- rep(NA,length=length(codNames_ordered_by_peak_keep))
AUC_less <- rep(NA,length=length(codNames_ordered_by_peak_keep))
AUC_total <- rep(NA,length=length(codNames_ordered_by_peak_keep))
AUCperlow <- rep(NA,length=length(codNames_ordered_by_peak_keep))
AUCperhigh <- rep(NA,length=length(codNames_ordered_by_peak_keep))

AUCthreshold <- c(0.1,0.2,1)
plottitles <- c("threshold freq=0.1","threshold freq=0.2","threshold freq=peak spawning age")

par(mfrow=c(4,3))
for (j in 1:length(AUCthreshold)){ #step through different AUC thresholds
  
  for (i in 1:length(codNames_ordered_by_peak_keep)){
    if(AUCthreshold[j] == 1) 
    {AUCthreshold[j] <- (1/eigentable[eigentable$codNames == codNames_ordered_by_peak_keep[i],]$mode_age)} else
    {AUCthreshold[j]}
    
    AUC_greater[i] <- sum(freq[1]*specdatalong[specdatalong$variable.type == "recruit" 
                                               & specdatalong$codNames==codNames_ordered_by_peak_keep[i]
                                               & specdatalong$freq > AUCthreshold[j],]$value)
    
    AUC_less[i] <- sum(freq[1]*specdatalong[specdatalong$variable.type == "recruit" 
                                            & specdatalong$codNames==codNames_ordered_by_peak_keep[i]
                                            & specdatalong$freq <= AUCthreshold[j],]$value)
    if(j==3){AUC_less[i] <- AUC_less[i]*AUCthreshold[j]} else {AUC_less[i]}
    if(j==3){AUC_greater[i] <- AUC_greater[i]*(0.5-AUCthreshold[j])} else {AUC_greater[i]}
    AUC_total[i] <- sum(freq[1]*specdatalong[specdatalong$variable.type == "recruit" 
                                             & specdatalong$codNames==codNames_ordered_by_peak_keep[i],]$value)
    AUCperlow[i] <- AUC_less[i]/AUC_total[i]
    AUCperhigh[i] <- AUC_greater[i]/AUC_total[i]
    print(AUCperlow[i]+AUCperhigh[i])
  }
  plot(x=eigentable.trunc[order(eigentable.trunc$mode_age),]$mode, 
       xlab="peak spawning age",y=AUCperlow,ylab="Percent AUC at low freq",main=plottitles[j])
  plot(x=eigentable.trunc[order(eigentable.trunc$mode_age),]$max_ages, 
       xlab="max age",y=AUCperlow,ylab="Percent AUC at low freq",main=plottitles[j])
  plot(x=eigentable.trunc[order(eigentable.trunc$mode_age),]$cvs_mode, 
       xlab="Spawning biomass distribution CV",y=AUCperlow,
       ylab="Percent AUC at low freq",main=plottitles[j])
  #plist <- as.list(p1peak,p2max,p3CV)
  #pList[[j]] <- plist
  
}
par(mfrow=c(1,1))

par(mfrow=c(2,2))
plot(x=eigentable.trunc[order(eigentable.trunc$mode_age),]$mode, 
     xlab="peak spawning age",y=AUC_total,ylab="Total AUC")
plot(x=eigentable.trunc[order(eigentable.trunc$mode_age),]$max_ages, 
     xlab="max age",y=AUC_total,ylab="Total AUC")
plot(x=eigentable.trunc[order(eigentable.trunc$mode_age),]$cvs_mode, 
     xlab="Spawning biomass distribution CV",y=AUC_total,
     ylab="Total AUC")
par(mfrow=c(1,1))

# *************************************** #
# (11) Area under curve - shaded boxes
# *************************************** #

ggplot(data, aes(variable, id)) +
  geom_raster(aes(fill = value)) + 
  scale_fill_gradient(low = "white",
                      high = "steelblue")
