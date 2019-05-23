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

# ---
# load the simulation model
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/simulation_model_cod_v3.r")

# load functions
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/2_cod_functions.r")

# load peak spawning age info
eigentable = read.csv("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLSB/eigentable.csv",
                      header=TRUE,stringsAsFactors = FALSE)
eigentable = as.data.frame(eigentable)

# load cod data, break into separate populations
#source("C:/Users/provo/Documents/GitHub/popdy/cod_code/0_load_cod_data.r")
# I don't think I need the actual data since I'm pulling maturity 
# information from Wang et al. 
#codNames <- c("Northsea","Coas","W_Baltic",
#              "Faroe","NE_Arctic","Celtic",
#              "Iceland","Kat","W_Scotland",
#              "NGulf","GB","GM",
#              "cod3NO","cod3M","cod2J3KL",
#              "cod3Ps")
codNames <- eigentable$codNames

# *************************************** #
# (Fig 1:distributions) spawning distributions for all pops
# see script -- 3_plot_spawning_biomass_distributions.r
# *************************************** #

# *************************************** #
# (Fig 4:eigenvalue) Generate Leslie matricies for diff F values (create Leslie arrays)
# *************************************** #
kvals = seq(from=0.1,to=1,by=0.1) #Check max F value (some pops can withstand high F)
Aarray = as.list(rep(NA,length(codNames))) #Leslie matrix storage for each F value for pop i
names(Aarray) <- codNames
eigenvals1 = matrix(NA,nrow=length(kvals),ncol=length(codNames)) 
eigenvals2 = matrix(NA,nrow=length(kvals),ncol=length(codNames))
eigenvals12 = matrix(NA,nrow=length(kvals),ncol=length(codNames))

for (i in 1:length(codNames)){ #for each pop i
  # load parms for cod pop i: L_inf, K (for vonB), TEMP, maxage,B0,B1 (matur)
  source(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_pops/',codNames[i], '.r', sep=''))
  
  Lesliearray <- array(NA,c(maxage,maxage,length(kvals))) #store Leslie matricies
  e1 = rep(NA,length=length(kvals)) #store lambda1
  e2 = rep(NA,length=length(kvals)) #store lambda2
  e12 = rep(NA,length=length(kvals)) #store inverse damping ratio
  
  for (k in 1:length(kvals)){ #step through F values 
    # create Leslie matrix:
    Leslieout = assemble_Leslie(maxage=maxage, K=K, L_inf=L_inf, TEMP=TEMP,
                                F.halfmax=0, B0=B0, B1=B1, tknot=0)
    # move survivals along subdiagonal to fecundities
    Leslieout$A[1,] <- Leslieout$A[1,]*Leslieout$A[2,1]
      #I can do it this way bc survival is constant w/age (ie no fishing)
    
    # set suvivals on subdiagonal = 1 (only works if no fishing)
    Leslieout$A[Leslieout$A == Leslieout$A[2,1]] <- 1
    
    # transform fecundity-at-age to probability density curve & multiply by k (slope)
    Leslieout$A[1,] <- (Leslieout$A[1,]/sum(Leslieout$A[1,]))*kvals[k]
    
    #NEAR$Survship = 0 # set up column for survivorship (amount or fraction present at age)
    #NEAR$Survship[1] = 1
    
    #for(k in 1:(nrow(NEAR)-1)){ # step through ages to calc survivorship
    #  NEAR$Survship[k+1] = NEAR$Survship[k]*NEAR$SURV_at_age[k] #amount present at age
    #}
    
    Lesliearray[,,k] = Leslieout$A #3D array of Leslie matricies
    e1[k] = extract_first_eigen_value(Leslieout$A)
    e2[k] = extract_second_eigen_value(Leslieout$A)
    e12[k] = e2[k] / e1[k]
    }
    
  Aarray[[i]]= Lesliearray #store Leslie 3D array in list of all pops
  eigenvals1[,i] = e1 #store lambda1 
  eigenvals2[,i] = e2 #store lambda2
  eigenvals12[,i] = e12 #store inverse of damping ratio
}
rm(e1,e2,e12,i,k,Lesliearray) #clean up

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
eigendata$peak <- eigentable[match(eigendata$codNames,eigentable$codNames),"mode_age"]

# (1) test significance of relationship between peak age & e1 
# (1a) within each k value: 
m.e1 <- as.list(rep(NA,length=length(kvals)))
for(n in 1:length(kvals)){
  dd <- eigendata[eigendata$kvals == kvals[n] & eigendata$eigen == "e1",]
  m.e1[[n]] <- lm(dd$value ~ dd$peak)
}
rm(n,dd)
e1_vs_peak_withinkval <- stargazer(m.e1[[1]],m.e1[[4]],m.e1[[6]],m.e1[[9]],dep.var.labels = c("lambda1"),
                               title="Regression Results for e1~peak within kvals",
                               covariate.labels=c("peak age"),type="text", report='vc*p',
                               column.labels = c("k=0.1","k=0.4","k=0.6","k=0.9"),
                               out="C:/Users/provo/Documents/GitHub/popdy/cod_figures/model_e1vspeak_withink.txt")
# (1b) across k values --> average within k, for each pop
# (1b1) average within k
meanatk_e1 <- rep(NA,length=length(kvals))
for(n in 1:length(kvals)){
  meanatk_e1[n] <- mean(eigendata[eigendata$kvals == kvals[n] & eigendata$eigen == "e1",]$value)}
e1_vs_peak_acrossk <- lm(meanatk_e1~kvals)

# (2) test significance of relationship between peak age & e12
m.e12 <- as.list(rep(NA,length=length(kvals)))
for(n in 1:length(kvals)){
  dd <- eigendata[eigendata$kvals == kvals[n] & eigendata$eigen == "e12",]
  m.e12[[n]] <- summary(lm(dd$value ~ dd$peak))
}
e12_vs_peak_withinkval <- stargazer(m.e12[[1]],m.e12[[4]],m.e12[[6]],m.e12[[9]],dep.var.labels = c("damping ratio"),
                                   title="Regression Results for e12~peak within kvals",
                                   covariate.labels=c("peak age"),type="text", report='vc*p',
                                   column.labels = c("k=0.1","k=0.4","k=0.6","k=0.9"),
                                   out="C:/Users/provo/Documents/GitHub/popdy/cod_figures/model_e1vspeak_withink.txt")
meanatk_e12 <- rep(NA,length=length(kvals))
for(n in 1:length(kvals)){
  meanatk_e12[n] <- mean(eigendata[eigendata$kvals == kvals[n] & eigendata$eigen == "e12",]$value)}
e12_vs_peak_acrossk <- summary(lm(meanatk_e12~kvals))



# *************************************** #
# (Fig 4) Eigenvalue 4 panel plot
# *************************************** #

lambda1 <- ggplot(eigendata[eigendata$eigen=="e1" & eigendata$kvals %in% c(0.1,0.4,0.6,0.9),],aes(x=peak,y=value)) +
  geom_point() + 
  geom_smooth(method="lm",se=FALSE,color="black") +
  facet_grid(. ~ kvals) +
  scale_y_continuous(limits=c(0,1.1)) +
  geom_text_repel(data=eigendata[eigendata$eigen=="e1"& eigendata$kvals %in% c(0.1,0.4,0.6,0.9),],
                  aes(label = codNames),
                  segment.color = "grey",
                  size = 2,
                  na.rm = TRUE) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(angle = 0)) +
  ylab("a") +
  xlab("Peak spawning age") +
  #ylab(expression(paste(lambda[1]))) +
  ylim(0.5,1.01)

lambda12 <- ggplot(eigendata[eigendata$eigen=="e12" & eigendata$kvals %in% c(0.1,0.4,0.6,0.9),],aes(x=peak,y=value)) +
  geom_point() + 
  geom_smooth(method="lm",se=FALSE,color="black") +
  facet_grid(. ~ kvals) +
  scale_y_continuous(limits=c(0,1.1)) +
  geom_text_repel(data=eigendata[eigendata$eigen=="e12"& eigendata$kvals %in% c(0.1,0.4,0.6,0.9),],
                  aes(label = codNames),
                  segment.color = "grey",
                  size = 2,
                  na.rm = TRUE) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(angle = 0)) +
  ylab("b") +
  xlab("Peak spawning age") +
  #ylab(expression(paste(abs(lambda[2])/lambda[1]))) +
  ylim(0.5,1.01)
p <- list(lambda1,lambda12)
tiff(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/manuscript/fig4_2xpanelplot.tiff', units="in", width=5, height=7, res=300)
do.call(grid.arrange,c(p,ncol=1))
dev.off()
rm(p)
# *************************************** #
# (Fig 2:schematic) Choose alpha values, plot BH curves
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
#betasH <- seq(from=10^5,to=10^7,by=((10^7)-(10^5))/4)
#betas <- seq(from=10^3,to=10^5,by=((10^5)-(10^3))/6)

#a = c(1.2,2,5,10) #alpha values
a = seq(from=1.1, to=10, by=0.1)
b = 1000 # beta
ee <- seq(from=0,to=b,by=b/10000) #range of egg values

Rlist <- as.list(rep(NA,length=length(a))) #recruit values for each alpha
names(Rlist) <- a
for(i in 1:length(a)){ Rlist[[i]] <- BH(alpha=a[i],beta=b,E=ee) }
RvE <- as.data.frame(do.call(cbind,Rlist))
RvE$eggs <- ee
RvElong <- RvE %>% gather(alpha,value,1:length(a))

# For each alpha, calculate the intersection with the 1:1 line. Find when R=E
eq <- rep(NA,length=length(a))
interpt <- rep(NA,length=length(a))
slopes <- rep(NA,length=length(a))
for(i in 1:length(a)){
  ee <- seq(from=0,to=b,by=b/10000)
  rr <- BH(alpha=a[i],beta=b,E=ee)
  dat <- as.data.frame(cbind(ee,rr))
  names(dat) <- c('eggs','recruits')
  dat$diff <- abs(dat$eggs - dat$recruits) #diff between R and E
  dat <- dat[-1,] #remove first row at origin
  row.names(dat) <- NULL
  #head(dat)
  #plot(dat$diff,type="l")
  eq[i] <- which.min(dat$diff )
  interpt[i] <- dat[eq[i],]$eggs
  #dat[(eqpoint[b]-3):(eqpoint[b]+3),]
  # calculate slope around equal point
  x1 <- dat[(eq[i]-1),]$eggs
  y1 <- dat[(eq[i]-1),]$recruits
  x2 <- dat[(eq[i]+1),]$eggs
  y2 <- dat[(eq[i]+1),]$recruits
  slopes[i] <- (y2-y1)/(x2-x1)
}
datalpha <- as.data.frame(cbind(a,interpt,slopes))
datalpha$slopesR <- round(datalpha$slopes,digits=3)
#datalpha[1:20,]
plot(x=datalpha$a,y=datalpha$slopesR)

plot_these_alpha <- datalpha[datalpha$a %in% c(1.1,2.5,5,10),]

RvElong_forplotting <- RvElong[RvElong$alpha %in% plot_these_alpha$a,]
# add column with intercept point for each alpha
RvElong_forplotting$interpt <- plot_these_alpha[match(RvElong_forplotting$alpha,plot_these_alpha$a),"interpt"]
top <- ggplot(data=RvElong_forplotting,aes(x=eggs,y=value,linetype=alpha)) + geom_line() +
  geom_dl(aes(label=alpha),method="last.points") +
  geom_abline(intercept=0, slope=1, color="black",size=1) + theme_classic() + 
  ylab("Recruits") + xlab("Egg production") +
  geom_point(aes(x=interpt,y=interpt,group=alpha),size=3) +
  theme(legend.position = "none") 
# FIG 2(b) - Schematic showing slope at intersection for multiple alpha values
bottom <- ggplot(data=datalpha,aes(x=a,y=slopes)) + geom_line() +
  geom_point(data=plot_these_alpha,aes(x=a,y=slopes),size=3) +
  geom_text(data=plot_these_alpha,aes(label=a),nudge_x=0.3,nudge_y=0.04) +
  xlab(expression(alpha)) +
  ylab("Slope of the Beverton-Holt\n curve at the intersection with 1/LEP") +
  theme_classic()
# export figure
plist <- list(top,bottom)
tiff(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/manuscript/fig2_schematic.tiff', units="in", width=4, height=7, res=300)
#jpeg(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/Fig1_schematic.jpg') #note: file name specifies the alpha used in simluation model
do.call(grid.arrange,c(plist,ncol=1))
dev.off()

# -------
# Alternative Fig for schematic - showing one BH curve and multiple 1/LEP
# -------
b=1000
a=2
ee <- seq(from=0,to=b*10,by=b/10000) #range of egg values
rr <- BH(alpha=a,beta=b,E=ee)
dd <- as.data.frame(cbind(ee,rr))
# calculate x and y end points for segments on plot (1/LEP line)
y1 = rep(0,length=length(datalpha[,1]))
y2 = rep(1000,length=length(datalpha[,1]))
x1 = rep(0,length=length(datalpha[,1]))
x2 = c(y2/datalpha$slopes)
df <- data.frame(cbind(y1,y2,x1,x2))

ggplot(dd,aes(x=ee,y=rr)) + geom_line() +
  ylim(0,1500) +
  geom_segment(data=df,aes(x=x1,y=y1,xend=x2,yend=y2)) +
  theme_classic()




# *************************************** #
# (3) Simulate pops. Loop over Aarray list to simulate using different Leslie matrices 
# *************************************** #
# set params for simulation:
timesteps = 1000 #need this now to create
rm_first_timesteps = 200
betas = 1000
alphas <- seq(from=0.1, to=10, by=0.1) #these correspond to 
sig_r = 0.3
span.multiplier = 1 # adjusting the span in spec.prgm()
#alphas <- rep(alpha, length=length(codNames)) #alpha could be diff for pops


output.3d.list <- as.list(rep(NA,length=length(codNames))) #store timeseries here
names(output.3d.list) <- codNames

for (i in 1:length(Aarray)) { #step through each pop
  Leslie3d = Aarray[[i]] #select the 3d array of Leslie matricies
  # array dims: row=ts length, col=4 is number of ts (eggs,recruits,Nt,Nsize), depth=F vals
  output.matrix <- array(NA,c(timesteps-2,4,length(alphas))) 
  
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
    output = sim_model(A=Leslie3d[,,10], timesteps=timesteps, 
                       alpha=alphas[a], beta=betas, 
                       sig_r=sig_r, initial_eggs=betas)
    
    length(output$Nsize) <- length(output$N_t) #trim Nsize ts vector, -2 elements
    output.matrix[,,a] <- do.call(cbind,output) #fill in array for pop i
    #colnames(output.matrix) <- names(output)
  }
  
  output.3d.list[[i]] <- output.matrix
}
rm(i,a,Leslie3d,output.matrix,output) #clean up

# At this point I have one important object:
# 1. [output.3d.list] a list of 3d arrays. Each array is timeseries output
#    from simulations at different F levels. 


# *************************************** #
# (3) Format output ts for plotting simulations using output.3d.list
# *************************************** #
variable_type <- c("Nt","eggs","recruits","Nsize")

# --- reorganize egg timeseries data --- #
var.number <- 2 # eggs
df.list <- as.list(rep(NA,length=length(codNames)))
names(df.list) <- codNames
for (i in 1:length(output.3d.list)) {
  # first, reformat data to work with ggplot
  aa <- as.data.frame(output.3d.list[[i]][,var.number,])
  aa$year <- seq(from=1, to=length(aa[,1]))
  colnames(aa) <- c(alphas,"year")
  aa1 <- aa %>% gather(alphavalue,value,1:length(alphas))
  aa1$variable <- rep(variable_type[var.number],length=length(aa1[,1]))
  aa1$codNames <- rep(codNames[i],length=length(aa[,1]))
  df.list[[i]] <- aa1
  rm(aa1,aa)}
eggs.ts <- do.call(rbind,df.list)
rm(df.list,i) #clean up

# --- reorganize recruit timeseries data --- #
var.number <- 3 # recruits
df.list <- as.list(rep(NA,length=length(codNames)))
names(df.list) <- codNames
for (i in 1:length(output.3d.list)) {
  # first, reformat data to work with ggplot
  aa <- as.data.frame(output.3d.list[[i]][,var.number,])
  aa$year <- seq(from=1, to=length(aa[,1]))
  colnames(aa) <- c(alphas,"year")
  aa1 <- aa %>% gather(alphavalue,value,1:length(alphas))
  aa1$variable <- rep(variable_type[var.number],length=length(aa1[,1]))
  aa1$codNames <- rep(codNames[i],length=length(aa[,1]))
  df.list[[i]] <- aa1
  rm(aa1,aa)}
recruits.ts <- do.call(rbind,df.list)
rm(df.list) #clean up

# --- reorganize Nsize timeseries data --- #
var.number <- 4 # Nsize
df.list <- as.list(rep(NA,length=length(codNames)))
names(df.list) <- codNames
for (i in 1:length(output.3d.list)) {
  # first, reformat data to work with ggplot
  aa <- as.data.frame(output.3d.list[[i]][,var.number,])
  aa$year <- seq(from=1, to=length(aa[,1]))
  colnames(aa) <- c(alphas,"year")
  aa1 <- aa %>% gather(alphavalue,value,1:length(alphas))
  aa1$variable <- rep(variable_type[var.number],length=length(aa1[,1]))
  aa1$codNames <- rep(codNames[i],length=length(aa[,1]))
  df.list[[i]] <- aa1
  rm(aa1,aa)}
nsize.ts <- do.call(rbind,df.list)
rm(df.list) #clean up

ts.data <- rbind(eggs.ts,recruits.ts,nsize.ts) #combine data
rownames(ts.data) <- NULL
head(ts.data)

# Put in appendix:
# For each cod population, calculate CV of time series (for recruits)
# and plot it against the corresponding beta value
pp <- as.list(rep(NA,length=length(codNames)))
for(i in 1:length(codNames)){
  dat <- ts.data[ts.data$codNames == codNames[i] & ts.data$variable == "recruits",]
  cvs <- rep(NA,length=length(alphas))
  for(b in 1:length(alphas)){
    vals <- dat[dat$alphavalue == alphas[b],]$value[rm_first_timesteps:(timesteps-2)]
    cvs[b] <- sd(vals)/mean(vals)
  }
  cvdat <- as.data.frame(cbind(cvs,alphas))
  pp[[i]] <- ggplot(cvdat,aes(x=alphas,y=cvs)) + geom_line() + 
    ggtitle(paste(codNames[i],"(b=1000)")) +
    ylab("") + xlab("") + ylim(c(0,0.1))
}

tiff(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/manuscript/figS_cv_vs_alpha.tiff', 
     units="in", width=7, height=8, res=300) 
do.call(grid.arrange,c(pp,ncol=4,left="time series CV",bottom="alpha value"))
dev.off()


# *************************************** #
# (4) Now that timeseries data is formated, let's plot! --- This section is option, skip to section 5 to calculate spectra.
# *************************************** #

# Put in appendix:
#plot recruitment - one plot per pop, similar to egg plots
#time series for 4 k values 
selectedalphas <- factor(c(1.1,1.4,2,3.3,10))
prec <- list()
codNames_ordered_by_peak <- eigentable %>% arrange(mode_age) %>% pull(codNames)
codNames_ordered_by_peak_plot <- eigentable %>% arrange(mode_age) %>% pull(codNames_plot)

for (i in 1:length(codNames_ordered_by_peak)){
  dd <- ts.data[ts.data$variable == "recruits" & 
                  ts.data$codNames == codNames_ordered_by_peak[i] &
                  ts.data$year %in% seq(from=rm_first_timesteps,to=(timesteps-2),by=1) &
                  ts.data$alphavalue %in% selectedalphas,]
  dd$kval <- as.character(round(1/as.numeric(dd$alphavalue),digits=2))
  
  prec[[i]] <- ggplot(dd, #aes(x=year,y=value,color=Fval)) +
                   aes(x=year,y=value,color=kval)) +
    xlab("") + ylab("") +
    geom_line() + theme_classic() + ylim(c(75,1050)) +
    #scale_color_brewer(palette = "Reds") +
    ggtitle(paste(codNames_ordered_by_peak_plot[i]))
}
tiff(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/manuscript/figS_timeseries_recruits_sigmaR0.3.tiff', units="in", width=7, height=13, res=300) 
do.call(grid.arrange,c(prec,ncol=2,left="Recruits (before noise)", bottom="Year"))
dev.off()
rm(selectedalphas,prec)

selectedalphas <- factor(c(1.1,1.4,2,3.3,10))
prec <- list()
theseones <- c("Celtic","W_Baltic","cod2J3KL","NE_Arctic")
for (i in 1:length(theseones)){
  dd <- ts.data[ts.data$variable == "recruits" & 
                  ts.data$codNames == theseones[i] &
                  ts.data$year %in% seq(from=rm_first_timesteps,to=600,by=1) &
                  ts.data$alphavalue %in% selectedalphas,]
  dd$kval <- as.character(round(1/as.numeric(dd$alphavalue),digits=2))
  
  prec[[i]] <- ggplot(dd, #aes(x=year,y=value,color=Fval)) +
                      aes(x=year,y=value,color=kval)) +
    xlab("") + ylab("") +
    geom_line() + theme_classic() + ylim(c(75,1050)) +
    scale_color_brewer(palette = "Reds") +
    labs(subtitle = paste(theseones[i]))
    #ggtitle(paste(theseones[i]))
}
tiff(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/manuscript/fig4a_timeseries_sigR0.3.tiff', units="in", width=4, height=6, res=300) 
do.call(grid.arrange,c(prec,ncol=1,left="Recruits (before noise)", bottom="Year"))
dev.off()
rm(selectedalphas,prec)

my_palette = brewer.pal(5, "Reds")
# *************************************** #
# (5) Calculate frequency content from timeseries
# *************************************** #
# Plan:
# 1. Walk through each cod pop, do spectral analysis at alpha levels
# 2. Store spec values for eggs, recruits, and Nsize

# 1. Walk through each cod pop, do spectral analysis at alpha levels
sp.eggsL <- as.list(rep(NA,length=length(codNames))) #object for spec analysis  
sp.recruitL <- as.list(rep(NA,length=length(codNames))) 
#ts.for.spec.eg <- as.list(rep(NA,length=length(codNames)))
#ts.for.spec.re <- as.list(rep(NA,length=length(codNames)))

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
    yy = ts[ts$variable == "recruits" & ts$alphavalue == alphas[b],]$value[rm_first_timesteps:(timesteps-2)] - mean(ts[ts$variable == "recruits" & ts$alphavalue == alphas[b],]$value[rm_first_timesteps:(timesteps-2)])
    sp = spec.pgram(yy,spans=c(m,m),plot = FALSE)
    spsaveL[[b]] = 2*sp$spec # save matrix of spec values for different FLEP, index by pop i
  }
  spsave <- as.data.frame(do.call(cbind,spsaveL))
  spsave$freq <- sp$freq
  spsavelong <- spsave %>% gather(alphavalue, value, 1:length(alphas))
  spsavelong$codNames <- rep(codNames[i],length=length(spsavelong[,1]))
  sp.recruitL[[i]] <- spsavelong
  rm(spsave)
  print(i)
}
#eggs <- do.call(rbind,sp.eggsL)
#eggs$variable.type <- rep("eggs",length=length(eggs$freq))
rec <- do.call(rbind,sp.recruitL)
rec$variable.type <- rep("recruits",length=length(rec$freq))
specdatalong <- rec
head(specdatalong)

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
  # store plots in list
  prec_sp[[i]] <- ggplot(data=dataforplot, aes(x=freq,y=value,group=alphavalue)) + 
    geom_line(aes(color=alphavalue)) + ylim(0,50000) +
    geom_vline(xintercept = (1/eigentable[eigentable$codNames == codNames_ordered_by_peak[i],]$mode_age),
               linetype="dotted") +
    geom_text(x=((1/eigentable[eigentable$codNames == codNames_ordered_by_peak[i],]$mode_age)+0.06), 
              y=40000, label=eigentable[eigentable$codNames == codNames_ordered_by_peak[i],]$mode_age, size=4) +
    ggtitle(paste(codNames_ordered_by_peak_plot[i])) + 
    theme_classic() + ylab("") + xlab("") + theme(legend.position = "none") +
    scale_colour_manual(values=cs) +
    theme(plot.title = element_text(size = 10)) 
    
}
names(prec_sp) <- codNames_ordered_by_peak
rm(i,dataforplot)

tiff(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/manuscript/figS_spectra_sigmaR0.3.tiff', units="in", width=7, height=7, res=300) 
do.call(grid.arrange,c(prec_sp,ncol=4,left="Spectra",bottom="Frequency"))
dev.off()

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
# (Fig 3) Spectra on one plot, color coded by peak
# spawning age, max age, and CV of spawning biomass distribution
# *************************************** #
# --- recruits spectra: all on one plot ---#
dataforplot <- specdatalong[specdatalong$variable.type == "recruits" &
                              specdatalong$alphavalue %in% c(10),] # select one alpha value
dataforplot$peak <- eigentable[match(dataforplot$codNames, eigentable$codNames),"mode_age"]
dataforplot$maxage <- eigentable[match(dataforplot$codNames, eigentable$codNames),"max_ages"]
dataforplot$sd_mode <- eigentable[match(dataforplot$codNames, eigentable$codNames),"sd_mode"]
dataforplot$cvs_mode <- eigentable[match(dataforplot$codNames, eigentable$codNames),"cvs_mode"]
dataforplot$codNames <- factor(dataforplot$codNames, 
                               levels=c("Coas","cod3M","cod3NO","cod3Ps","Northsea",
                                 "Faroe","GB","GM","Iceland","Kat","NGulf","W_Scotland",
                                 "Celtic","NE_Arctic","cod2J3KL","W_Baltic"))
j <- ggplot(dataforplot, aes(x=freq,y=value,group=codNames)) + 
  geom_line(aes(color=peak)) + theme_classic() + 
  ggtitle("a") +
  ylab("") + xlab("") +
  theme(plot.title = element_text(size = 12)) +
  guides(fill=guide_legend(title="Peak spawning age"))
  
jj <- ggplot(dataforplot, aes(x=freq,y=value,group=codNames)) + 
  geom_line(aes(color=maxage)) + theme_classic() + 
  ggtitle("b") +
  ylab("") + xlab("") +
  theme(plot.title = element_text(size = 12)) +
  guides(fill=guide_legend(title="Max age"))

jjj <- ggplot(dataforplot, aes(x=freq,y=value,group=codNames)) + 
  geom_line(aes(color=cvs_mode)) + theme_classic() + 
  guides(fill=guide_legend(title="CV")) +
  ggtitle("c") +
  ylab("") + xlab("") +
  theme(plot.title = element_text(size = 12)) 

#set all pops to one color
line.cols <- rep("grey",16) 
CN <- c("Coas","cod3M","cod3NO","cod3Ps","Northsea",
        "Faroe","GB","GM","Iceland","Kat","NGulf","W_Scotland",
        "Celtic","NE_Arctic","cod2J3KL","W_Baltic")
CN.t <- as.data.frame(cbind(CN,line.cols),stringsAsFactors = FALSE)
names(CN.t) <- c("pops","line.cols")
CN.t[CN.t$pops == "NE_Arctic",]$line.cols <- "tomato"
CN.t[CN.t$pops == "Celtic",]$line.cols <- "orange"
CN.t[CN.t$pops == "cod2J3KL",]$line.cols <- "dodgerblue"
CN.t[CN.t$pops == "W_Baltic",]$line.cols <- "green3"

jjjj <- ggplot(dataforplot) + 
  geom_line(aes(x=freq,y=value,color=codNames),size=1) + theme_classic() + 
  scale_color_manual(values=CN.t$line.cols) +
  ggtitle("b") + 
  ylab("Variance") + xlab("") +
  theme(plot.title = element_text(size = 14,hjust = -0.2),legend.position = "none" ) 

# create spectra plot for a=1.1 (depleted)
dataforplotD <- specdatalong[specdatalong$variable.type == "recruits" &
                              specdatalong$alphavalue %in% c(1.1),] # select one alpha value
dataforplotD$peak <- eigentable[match(dataforplotD$codNames, eigentable$codNames),"mode_age"]
dataforplotD$maxage <- eigentable[match(dataforplotD$codNames, eigentable$codNames),"max_ages"]
dataforplotD$sd_mode <- eigentable[match(dataforplotD$codNames, eigentable$codNames),"sd_mode"]
dataforplotD$cvs_mode <- eigentable[match(dataforplotD$codNames, eigentable$codNames),"cvs_mode"]
dataforplotD$codNames <- factor(dataforplotD$codNames, 
                               levels=c("Coas","cod3M","cod3NO","cod3Ps","Northsea",
                                        "Faroe","GB","GM","Iceland","Kat","NGulf","W_Scotland",
                                        "Celtic","NE_Arctic","cod2J3KL","W_Baltic"))

dddd <- ggplot(dataforplotD) + 
  geom_line(aes(x=freq,y=value,color=codNames),size=1) + theme_classic() + 
  scale_color_manual(values=CN.t$line.cols) +
  #geom_line(data=dataforplot[dataforplot$codNames=="NE_Arctic",],aes(x=freq,y=value)) +
  ggtitle("c") +
  ylab("Variance") + xlab("") +
  theme(plot.title = element_text(size = 14,hjust = -0.2),legend.position = "none" ) 
p <- list(jjjj,dddd)
tiff(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/manuscript/fig4bc_spectra_sigR0.3.tiff', units="in", width=4, height=7, res=300) 
do.call(grid.arrange,c(p,ncol=1,bottom="Frequency"))
dev.off()


# *************************************** #
# (Fig 5) Spectra plots arranged by k and peak spawning age
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

tiff(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/manuscript/fig5a_k_peakageorder_spectra_sigmaR0.1.tiff', units="in", width=10, height=7, res=300) 
ggplot(plotdat[plotdat$codNames == codNames[1:8],], aes(x=freq,y=value)) + 
  geom_line() + facet_grid(kval~codNames) + 
  ggtitle("Plot 1: populations ordered by peak spawning age") #ylim(0,1000)
dev.off()
tiff(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/manuscript/fig5b_k_peakageorder_spectra_sigmaR0.1.tiff', units="in", width=10, height=7, res=300) 
ggplot(plotdat[plotdat$codNames == codNames[9:16],], aes(x=freq,y=value)) + 
  geom_line() + facet_grid(kval~codNames) + 
  ggtitle("Plot 2: populations ordered by peak spawning age") +
  ylim(0,1000)
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

kvals <- unique(specdatalong$kval) #all possible kvals
alphas <- unique(specdatalong$alphavalue)

AUC_less_L <- as.list(rep(NA,length=length(alphas)))
AUCperlow_L <- as.list(rep(NA,length=length(alphas)))
AUC_total_L <- as.list(rep(NA,length=length(alphas)))

names(AUC_less_L) <- as.character(alphas)
names(AUC_total_L) <- alphas
names(AUCperlow_L) <- alphas

for (j in 1:length(alphas)){ #for each alpha (ie kval)...
  
  AUC_less <- rep(NA,length=length(codNames_ordered_by_peak))
  AUC_total <- rep(NA,length=length(codNames_ordered_by_peak))
  AUCperlow <- rep(NA,length=length(codNames_ordered_by_peak))
  
  for (i in 1:length(codNames_ordered_by_peak)){ #step through the pops
    
    AUC_less[i] <- sum(freq*specdatalong[specdatalong$variable.type == "recruits" 
                                  & specdatalong$codNames==codNames_ordered_by_peak[i]
                                  & specdatalong$alphavalue==alphas[j]
                                  & specdatalong$freq <= AUCthreshold_ordered_by_peak[i],]$value)
    
    AUC_total[i] <- sum(freq*specdatalong[specdatalong$variable.type == "recruits" 
                                          & specdatalong$codNames==codNames_ordered_by_peak[i]
                                          & specdatalong$alphavalue==alphas[j],]$value)
    AUCperlow[i] <- AUC_less[i]/AUC_total[i]
    
  }
  #store percents for each k value
  AUC_less_L[[j]]    <- AUC_less 
  AUC_total_L[[j]]   <- AUC_total
  AUCperlow_L[[j]]   <- AUCperlow
  print(j)
}
rm(i,j)
AUC_less_df <- data.frame(do.call(cbind,AUC_less_L))
AUC_total_df <- data.frame(do.call(cbind,AUC_total_L))
AUCperlow_df <- data.frame(do.call(cbind,AUCperlow_L))

AUC_less_df$codNames <- codNames_ordered_by_peak
AUC_total_df$codNames <- codNames_ordered_by_peak
AUCperlow_df$codNames <- codNames_ordered_by_peak


# convert dfs to long format
AUC_less_dflong <- AUC_less_df %>% gather(alpha,value,1:length(alphas)) %>% separate(alpha,c("addedX","alphaval"),sep="X") %>% select(-addedX) %>% mutate(AUCdes=rep("less"))
AUC_total_dflong <- AUC_total_df %>% gather(alpha,value,1:length(alphas)) %>% separate(alpha,c("addedX","alphaval"),sep="X") %>% select(-addedX) %>% mutate(AUCdes=rep("total"))
AUCperlow_dflong <- AUCperlow_df %>% gather(alpha,value,1:length(alphas)) %>% separate(alpha,c("addedX","alphaval"),sep="X") %>% select(-addedX) %>% mutate(AUCdes=rep("perlow"))

AUCdat <- rbind(AUC_less_dflong,
                AUC_total_dflong,
                AUCperlow_dflong)

AUCdat$peak <- eigentable[match(AUCdat$codNames,eigentable$codNames),"mode_age"]
AUCdat$maxage <- eigentable[match(AUCdat$codNames,eigentable$codNames),"max_ages"]
AUCdat$alphaval <- as.numeric(AUCdat$alphaval)
AUCdat$kval <- round(1/AUCdat$alphaval,digits = 2)
AUCdat <- AUCdat[AUCdat$kval <1,]
# set k slopes to factor, order by increasing slope
AUCdat$kval <- factor(AUCdat$kval,levels=unique(AUCdat$kval))
AUCdat$kval <- factor(AUCdat$kval,levels=rev(levels(unique(AUCdat$kval))))




# *************************************** #
# (Fig 5a) Total variance for different k values
# (Fig 5b) % of low frequency variance
# (Fig 5c) Difference in low frequency bars
# *************************************** #
head(ts.data)
ts.data$alphavalue <- as.numeric(as.character(ts.data$alphavalue))
varL = as.list(rep(NA,length=length(plotalpha)))
for(i in 1:length(codNames)){
  dat <- ts.data[ts.data$codNames == codNames[i] & 
                 ts.data$variable == "recruits" & 
                 ts.data$alphavalue %in% plotalpha,]
  variance <- rep(NA,length=length(plotalpha))
  means <- rep(NA,length=length(plotalpha))
  
  for(b in 1:length(plotalpha)){ #step through alpha values
    
    # calculate equilibrium value for each time series (rm first 200 ts)
    means[b] <- mean(dat[dat$alphavalue == plotalpha[b],]$value[rm_first_timesteps:(timesteps-2)])
    # substract mean from time series (vals_meanrm)
    vals = dat[dat$alphavalue == plotalpha[b],]$value[rm_first_timesteps:(timesteps-2)] - means[b]
    # sq root the variance & then divide by mean 
    variance[b] <- sqrt(var(vals))/means[b]
  }
  varL[[i]] <- as.data.frame(cbind(variance,plotalpha,rep(codNames[i],length=length(plotalpha))))
  names(varL[[i]]) <- c("variance","alphavalue","codNames")
}
vardat <- do.call(rbind,varL)
rm(i,b,dat,variance,means,vals) #clean up
vardat$peak <- eigentable[match(vardat$codNames,eigentable$codNames),"mode_age"]
vardat$kval <- round((1/as.numeric(as.character(vardat$alphavalue))),digits = 2)
vardat$codNames_plot <- eigentable[match(vardat$codNames,eigentable$codNames),"codNames_plot"]
# set k slopes to factor, order by increasing slope
vardat$kval <- factor(vardat$kval,levels=unique(vardat$kval))
#vardat$kval <- factor(vardat$kval,levels=rev(levels(unique(vardat$kval))))
# make sure cv values are numeric
vardat$variance <- as.numeric(as.character(vardat$variance))
# set codNames factor levels to increase with peak age
vardat$codNames_plot <- factor(vardat$codNames_plot,levels=codNames_ordered_by_peak_plot)
# set codNames in order of peak spawning age
vardat$codNames <- factor(vardat$codNames,levels=codNames_ordered_by_peak)

# get range of total variance at least depleted (k=0.1)
max.var.k0.1 <- max(vardat[vardat$kval == 0.1,]$variance)
min.var.k0.1 <- min(vardat[vardat$kval == 0.1,]$variance)
# get range of total variance at most depleted (k=0.9)
max.var.k0.9 <- max(vardat[vardat$kval == 0.91,]$variance)
min.var.k0.9 <- min(vardat[vardat$kval == 0.91,]$variance)


fig5a <- ggplot(vardat,aes(x=kval,y=codNames)) +
  geom_raster(aes(fill=variance)) + 
  xlab("Slope on S-R curve (k)") + ylab("") +
  scale_fill_gradient(low="purple", high="orange") + 
  theme_classic() +
  guides(fill=guide_legend(title="Total variance", reverse=TRUE)) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1),
        legend.position = "top",
        legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8))
tiff(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/manuscript/fig5a_totalvar_sigR0.3.tiff', units="in", width=4, height=6, res=300) 
fig5a
dev.off()
# (Fig 5b) Percent at low frequencies for different k values

# order pops by peak
AUCdat$codNames <- factor(AUCdat$codNames, 
                          levels=unique(AUCdat$codNames[order(AUCdat$peak)]))
plotalpha <- c(10,5,3.3,2.5,2,1.7,1.4,1.2,1.1)
# FIG 5E: Percent of AUC at low frequencies for different k values
dataforplot <- AUCdat[AUCdat$AUCdes == "perlow" &
                        AUCdat$alphaval %in% plotalpha,]
# get range of % low freq variance when populations are not depleted:
min.per.k0.1 <- min(dataforplot[dataforplot$kval == 0.1,]$value)
max.per.k0.1 <- max(dataforplot[dataforplot$kval == 0.1,]$value)
# get range of % high freq variance when pops are not depleted:
max.per.k.01high <- 1-min.per.k0.1
min.per.k.01high <- 1-max.per.k0.1
# get range of % low freq variance when populations are most depleted:
min.per.k0.9 <- min(dataforplot[dataforplot$kval == 0.91,]$value)
max.per.k0.9 <- max(dataforplot[dataforplot$kval == 0.91,]$value)
# get range of % high freq variance when pops are most depleted:
max.per.k0.9high <- 1-min.per.k0.9
min.per.k0.9high <- 1-max.per.k0.9

fig5b <- ggplot(data=dataforplot,aes(x=kval,y=codNames)) +
  geom_raster(aes(fill=value)) + 
  xlab("Slope on S-R curve (k)") + ylab("") + 
  scale_fill_gradient(low="purple", high="orange") + 
  scale_colour_gradient(limits = c(0, 1)) +
  theme_classic() +
  guides(fill=guide_legend(title="Low frequency variance\n(% range 0.52-0.99)")) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1),
        legend.position="top",
        plot.title = element_text(hjust = -0.1, vjust=-0.1),
        legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8)) 
  
tiff(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/manuscript/fig5b_lowvar_sigR0.3.tiff', units="in", width=4, height=6, res=300) 
fig5b
dev.off()
# calculate difference in % of low frequency variance when 
# k=0.9 (alpha=1.1) vs k=0.1 (alpha=10)
pervar_k0.1 <- AUCdat[AUCdat$AUCdes == "perlow" & AUCdat$alphaval == 10,]
pervar_k0.91 <- AUCdat[AUCdat$AUCdes == "perlow" & AUCdat$alphaval == 1.1,]
diff <- rep(NA,length=length(codNames_ordered_by_peak))
for(d in 1:length(codNames_ordered_by_peak)){
  
  diff[d] <- round(AUCdat[AUCdat$AUCdes == "perlow" & AUCdat$alphaval == 1.1 & AUCdat$codNames==codNames_ordered_by_peak[d],]$value - AUCdat[AUCdat$AUCdes == "perlow" & AUCdat$alphaval == 10 & AUCdat$codNames==codNames_ordered_by_peak[d],]$value, digits = 2)
  
}
diffdata <- as.data.frame(cbind(rev(codNames_ordered_by_peak),rev(diff)))
names(diffdata) <- c("codNames_ordered_by_peak","diff")
diffdata$codNames_ordered_by_peak <- factor(diffdata$codNames_ordered_by_peak,
                                            levels=codNames_ordered_by_peak)
diffdata$diff <- as.numeric(as.character(diffdata$diff))

# try averaging the difference in % populations with the same peak age
diffdata$peak <- eigentable[match(diffdata$codNames_ordered_by_peak,eigentable$codNames),"mode_age"]
plot(diffdata %>% group_by(peak) %>% summarise(avg = mean(diff)))

fig5c <- ggplot(data=diffdata, aes(x=codNames_ordered_by_peak, y=diff)) + geom_bar(stat="identity") + 
  scale_y_continuous(limits = c(0,0.17)) +
  theme_classic() +
  ggtitle("") +
  xlab("Percent") + ylab("") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
  coord_flip() 

tiff(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/manuscript/fig5c_bars_sigR0.3.tiff', units="in", width=4, height=6, res=300) 
fig5c
dev.off()






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