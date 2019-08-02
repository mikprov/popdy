# RUN
# by: mikaela provost
# last edited: Aug 1, 2019
# ===================================================================

library(tidyr)
library(dplyr)
library(gridExtra)

# *************************************** #
# Plan:
# (1) run script to calculate life table information and load parms df
# (2) assemble Leslie matrices for different F values
# (3) simulate time series 
# (4) calculate spectra



# *************************************** #
# (1) run script to get life table information and load parms df
source("C:/.../life_table_calculations.r")
# print lists of vectors if I want to verify they were successful
La_list
Ba_list 
eggs_list
prop_list
vul_list

# *************************************** #
# (2) assemble Leslie matrices for different F values
# getting ready
Fvals = seq(from=0,to=3,by=0.1) #instantaneous F rates
A3d_list = as.list(rep(NA,length(parms$spp))) #Leslie arrays in list
LTABLE_list <- as.list(rep(NA,length(parms$spp))) #LTABLE arrays in list
names(A3d_list) <- parms$spp
names(LTABLE_list) <- parms$spp
eigan1 <- as.matrix(NA,nrow=length(Fvals),ncol=length(parms$spp))
eigan12 <- as.matrix(NA,nrow=length(Fvals),ncol=length(parms$spp))

for(i in 1:length(parms$spp)){ #for spp i
  
  # define parms for spp i
  M = parms[parms$spp == parms$spp[i],]$M
  maxage = parms[parms$spp == parms$spp[i],]$maxage
  
  # empty arrays for Leslie and LTABLE, vectors for eigens
  A3d <- array(NA,c(maxage,maxage,length(Fvals)))
  LTABLE3d <- array(NA,c(maxage,12,length(Fvals))) #12 cols
  e1 <- rep(NA,length=length(Fvals))
  e12 <- rep(NA,length=length(Fvals))
  
  # assemble Leslie for each value of Fvals
  for(f in 1:length(Fvals)){
    out <- assemble_Leslie(maxage=maxage,
                           L_a=La_list[[i]],
                           B_a=Ba_list[[i]],
                           propmat_a=prop_list[[i]],
                           eggs_a=eggs_list[[i]],
                           vul_a=vul_list[[i]],
                           M=M,
                           Fval=Fvals[f])
    A3d[,,f] <- out$A           #store Leslies in 3d array
    LTABLE3d[,,f] <- out$LTABLE #store LTABLE in 3d array
    #some notes on A generated in assemble_Leslie(): function
    #has option to adjust LEP so that it is constant among spp
    #check out it out in functions.r
    
    # extract lambda1, lambda2/lambda1
    e1[f] = extract_first_eigen_value(out$A)
    e2[f] = extract_second_eigen_value(out$A)
    e12[f] = e2[f] / e1[f]
  }
  A3d_list[[i]] <- A3d 
  LTABLE_list[[i]] <- LTABLE3d
  eigan1[,i] <- e1
  eigen12[,i] <- e12
  #clean up before next i in loop
  rm(A3d,LTABLE3d,out,M,maxage) 
}
rm(i,j) #clean up




# *************************************** #
# (1) Generate Leslie matricies for diff F values (create Leslie arrays)
# *************************************** #
Fvalues = seq(0,0,by=1) #Check max F value (some pops can withstand high F)
Aarray = as.list(rep(NA,length(codNames))) #Leslie matrix storage for each F value for pop i
names(Aarray) <- codNames
eigenvals1 = matrix(NA,nrow=length(Fvalues),ncol=length(codNames)) 
eigenvals2 = matrix(NA,nrow=length(Fvalues),ncol=length(codNames))
eigenvals12 = matrix(NA,nrow=length(Fvalues),ncol=length(codNames))


for (i in 1:length(codNames)){ #for each pop i
  # load parms for cod pop i: L_inf, K (for vonB), TEMP, maxage,B0,B1 (matur)
  source(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_pops/',codNames[i], '.r', sep=''))
  
  Lesliearray <- array(NA,c(maxage,maxage,length(Fvalues))) #store Leslie matricies
  e1 = rep(NA,length=length(Fvalues)) #store lambda1
  e2 = rep(NA,length=length(Fvalues)) #store lambda2
  e12 = rep(NA,length=length(Fvalues)) #store inverse damping ratio
  
  for (f in 1:length(Fvalues)){ #step through F values 
    # create Leslie matrix:
    Leslieout = assemble_Leslie(maxage=maxage, K=K, L_inf=L_inf, TEMP=TEMP,
                                F.halfmax=Fvalues[f], B0=B0, B1=B1, tknot=0)
    # move survivals along subdiagonal to fecundities
    Leslieout$A[1,] <- Leslieout$A[1,]*Leslieout$A[2,1]
    
    # set suvivals on subdiagonal = 1 (only works if no fishing)
    Leslieout$A[Leslieout$A == Leslieout$A[2,1]] <- 1
    
    # transform fecundity-at-age to probability density curve
    Leslieout$A[1,] <- (Leslieout$A[1,]/sum(Leslieout$A[1,]))*1
    
    
    
    Lesliearray[,,f] = Leslieout$A #3D array of Leslie matricies
    e1[f] = extract_first_eigen_value(Leslieout$A)
    e2[f] = extract_second_eigen_value(Leslieout$A)
    e12[f] = eigenvals2[f] / eigenvals1[f]
    }
    
  Aarray[[i]]= Lesliearray #store Leslie 3D array in list of all pops
  eigenvals1[,i] = e1 #store lambda1 
  eigenvals2[,i] = e2 #store lambda2
  eigenvals12[,i] = e12 #store inverse of damping ratio
}
rm(e1,e2,e12,i,f,Lesliearray) #clean up


# *************************************** #
# (3) Simulate pops. Loop over Aarray list to simulate using different Leslie matrices 
# *************************************** #
# set params for simulation:
timesteps = 1000 #need this now to create
rm_first_timesteps = 200
alpha = 1.2
betas = c(10^7,20^7,30^7,40^7)
sig_r = 0.3
span.multiplier = 1 # adjusting the span in spec.prgm()
alphas <- rep(alpha, length=length(codNames)) #alpha could be diff for pops
#alphas <- c(2.02,8.38,4.1,2.10,3.77,2.10,3.24,4.3,5.54,0.48,1.69,0.77,2.31,0.42,1.22)

output.3d.list <- as.list(rep(NA,length=length(codNames))) #store timeseries here
names(output.3d.list) <- codNames

for (i in 1:length(Aarray)) { #step through each pop
  Leslie3d = Aarray[[i]] #select the 3d array of Leslie matricies
  # array dims: row=ts length, col=4 is number of ts (eggs,recruits,Nt,Nsize), depth=F vals
  output.matrix <- array(NA,c(timesteps-2,4,length(betas))) 
  
  for (b in 1:length(betas)) { #step through each Leslie matrix (for each F value)
    output = sim_model(A=Leslie3d[,,1], timesteps=timesteps, 
                       alpha=alphas[i], beta=betas[b], 
                       sig_r=sig_r, initial_eggs=betas[b])
    
    length(output$Nsize) <- length(output$N_t) #trim Nsize ts vector, -2 elements
    output.matrix[,,b] <- do.call(cbind,output) #fill in array for pop i
    #colnames(output.matrix) <- names(output)
  }
  
  output.3d.list[[i]] <- output.matrix
}
rm(i,b,Leslie3d,output.matrix,output) #clean up

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
  colnames(aa) <- c(betas,"year")
  aa1 <- aa %>% gather(betavalue,value,1:length(betas))
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
  colnames(aa) <- c(betas,"year")
  aa1 <- aa %>% gather(betavalue,value,1:length(betas))
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
  colnames(aa) <- c(betas,"year")
  aa1 <- aa %>% gather(betavalue,value,1:length(betas))
  aa1$variable <- rep(variable_type[var.number],length=length(aa1[,1]))
  aa1$codNames <- rep(codNames[i],length=length(aa[,1]))
  df.list[[i]] <- aa1
  rm(aa1,aa)}
nsize.ts <- do.call(rbind,df.list)
rm(df.list) #clean up

ts.data <- rbind(eggs.ts,recruits.ts,nsize.ts) #combine data
rownames(ts.data) <- NULL


# *************************************** #
# (4) Now that timeseries data is formated, let's plot! --- This section is option, skip to section 5 to calculate spectra.
# *************************************** #

#plot recruitment - one plot per pop, similar to egg plots
prec <- list()
codNames_ordered_by_peak <- eigentable %>% arrange(mode_age) %>% pull(codNames)
#sd_ordered <- rep(NA,length=length(codNames_ordered_by_peak))
#mean_ordered <- rep(NA,length=length(codNames_ordered_by_peak))
#CV_ordered <- rep(NA,length=length(codNames_ordered_by_peak))

for (i in 1:length(codNames_ordered_by_peak)){
  dd <- ts.data[ts.data$variable == "recruits" & 
                  ts.data$codNames == codNames_ordered_by_peak[i] &
                  ts.data$year %in% seq(from=rm_first_timesteps,to=(timesteps-2),by=1),]
  
  prec[[i]] <- ggplot(dd, #aes(x=year,y=value,color=Fval)) +
                   aes(x=year,y=value,color=betavalue)) +
    xlab("year") + ylab("recruits (before noise)") +
    geom_line() + theme_classic() + #ylim(c(1000,3000)) +
    #scale_color_brewer(palette = "Reds") +
    ggtitle(paste(codNames_ordered_by_peak[i]))
  
  #sd_ordered[i] <- sd(dd$value)
  #mean_ordered[i] <- mean(dd$value)
  #CV_ordered[i] <- sd_ordered[i]/mean_ordered[i]
}
#pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/recruitment_beforenoise_for_diff_FLEPs_alpha50_meanvsSTDEV.pdf', width=7, height=10) #note: file name specifies the alpha used in simluation model
do.call(grid.arrange,c(prec,ncol=3))
#dev.off()



# plot eggs - one plot per pop, on each plot 5 lines for different F levels
peggs <- list()
# reorder codNames by peak spawning age (increasing)
#codNames_ordered_by_peak <- eigentable %>% arrange(mode_age) %>% pull(codNames)
#sd_ordered <- rep(NA,length=length(codNames_ordered_by_peak))
#mean_ordered <- rep(NA,length=length(codNames_ordered_by_peak))
#CV_ordered <- rep(NA,length=length(codNames_ordered_by_peak))

for (i in 1:length(codNames_ordered_by_peak)){
  
  dd <- ts.data[ts.data$variable == "eggs" & 
                  ts.data$codNames == codNames_ordered_by_peak[i] &
                  ts.data$year %in% seq(from=rm_first_timesteps,to=(timesteps-2),by=1),]
  
  peggs[[i]] <- ggplot(dd, #aes(x=year,y=value,color=Fval)) +
                   aes(x=year,y=value,color=betavalue)) +
    xlab("year") + ylab("egg production") +
    geom_line() + theme_classic() + #ylim(c(0,600)) +
    #scale_color_brewer(palette = "Reds") +
    ggtitle(paste(codNames_ordered_by_peak[i],"var=",
                  round(var(dd$value),digits=0)))
    
  #sd_ordered[i] <- sd(dd$value)
  #mean_ordered[i] <- mean(dd$value)
  #CV_ordered[i] <- sd_ordered[i]/mean_ordered[i]
}
#pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/egg_production_for_diff_FLEPs_alpha50_meanvsStdev.pdf', width=7, height=10) #note: file name specifies the alpha used in simluation model
do.call(grid.arrange,c(peggs,ncol=3))


# *************************************** #
# (5) Calculate frequency content from timeseries
# *************************************** #
# Plan:
# 1. Walk through each cod pop, do spectral analysis at beta levels
# 2. Store spec values for eggs, recruits, and Nsize

# 1. Walk through each cod pop, do spectral analysis at beta levels
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
  spsaveL <- as.list(rep(NA,length=length(betas)))
  names(spsaveL) <- betas
  for (b in 1:length(betas)){
    xx = ts[ts$variable == "eggs" & ts$betavalue == betas[b],]$value[rm_first_timesteps:(timesteps-2)] - mean(ts[ts$variable == "eggs" & ts$betavalue == betas[b],]$value[rm_first_timesteps:(timesteps-2)])
    sp = spectrum(x=xx,spans=c(m,m),plot = FALSE)
    #save spec output for plotting, Helen Wearing says to multiply by 2
    spsaveL[[b]] <- sp$spec*2
  }
  spsave <- as.data.frame(do.call(cbind,spsaveL))
  spsave$freq <- sp$freq
  spsavelong <- spsave %>% gather(betavalue, value, 1:length(betas))
  spsavelong$codNames <- rep(codNames[i],length=length(spsavelong[,1]))
  sp.eggsL[[i]] <- spsavelong
  rm(spsave)
  
  # --- spectral analysis on RECRUIT --- #
  spsaveL <- as.list(rep(NA,length=length(betas)))
  names(spsaveL) <- betas
  for (b in 1:length(betas)){
    yy = ts[ts$variable == "recruits" & ts$betavalue == betas[b],]$value[rm_first_timesteps:(timesteps-2)] - mean(ts[ts$variable == "recruits" & ts$betavalue == betas[b],]$value[rm_first_timesteps:(timesteps-2)])
    sp = spec.pgram(yy,spans=c(m,m),plot = FALSE)
    spsaveL[[b]] = 2*sp$spec # save matrix of spec values for different FLEP, index by pop i
  }
  spsave <- as.data.frame(do.call(cbind,spsaveL))
  spsave$freq <- sp$freq
  spsavelong <- spsave %>% gather(betavalue, value, 1:length(betas))
  spsavelong$codNames <- rep(codNames[i],length=length(spsavelong[,1]))
  sp.recruitL[[i]] <- spsavelong
  rm(spsave)
  
}
eggs <- do.call(rbind,sp.eggsL)
eggs$variable.type <- rep("eggs",length=length(eggs$freq))
rec <- do.call(rbind,sp.recruitL)
rec$variable.type <- rep("recruits",length=length(rec$freq))
specdatalong <- rbind(eggs,rec)
head(specdatalong)

# *************************************** #
# (6) Plot spectral analysis 
# *************************************** #
# --- egg spectra --- # 
peggs_sp <- list()
for (i in 1:length(codNames_ordered_by_peak)){
  
  dataforplot <- specdatalong[specdatalong$variable.type == "eggs" 
                              & specdatalong$codNames==codNames_ordered_by_peak[i],]
  # store plots in list
  peggs_sp[[i]] <- ggplot(data=dataforplot, aes(x=freq,y=value,color=betavalue)) + 
    geom_line() + #ylim(1,10) +
    geom_vline(xintercept = (1/eigentable[eigentable$codNames == codNames_ordered_by_peak[i],]$mode_age),
             linetype="dotted") +
    ylab("log(eggs)") +
    ggtitle(paste(codNames_ordered_by_peak[i]," mode=",
                eigentable[eigentable$codNames == codNames_ordered_by_peak[i],]$mode_age," max=",
                eigentable[eigentable$codNames == codNames_ordered_by_peak[i],]$max_ages)) + 
    theme(plot.title = element_text(size = 6)) + theme_classic()
}
#pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/spec_eggs_alpha50_log10_nofishing.pdf', width=7, height=14) #note: file name specifies the alpha used in simluation model
do.call(grid.arrange,c(peggs_sp,ncol=3))
#dev.off()
rm(i,dataforplot)

# --- recruits spectra --- #
prec_sp <- list()
for (i in 1:length(codNames_ordered_by_peak)){
  dataforplot <- specdatalong[specdatalong$variable.type == "recruits" 
                              & specdatalong$codNames==codNames_ordered_by_peak[i],]
  # store plots in list
  prec_sp[[i]] <- ggplot(data=dataforplot,aes(x=freq,y=value,color=betavalue)) + 
    geom_line() + #ylim(1,5) +
    geom_vline(xintercept = (1/eigentable[eigentable$codNames == codNames_ordered_by_peak[i],]$mode_age),
               linetype="dotted") +
    ylab("log(recruits)") +
    ggtitle(paste(codNames_ordered_by_peak[i]," mode=",
                  eigentable[eigentable$codNames == codNames_ordered_by_peak[i],]$mode_age," max=",
                  eigentable[eigentable$codNames == codNames_ordered_by_peak[i],]$max_ages)) + 
    theme(plot.title = element_text(size = 10)) + theme_classic()
}
#pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/spec_recruitsbeforenoise_alpha50_log10.pdf', width=7, height=14) #note: file name specifies the alpha used in simluation model
do.call(grid.arrange,c(prec_sp,ncol=3))
#dev.off()
rm(i,dataforplot)

# --- recruits spectra: all on one plot ---#
dataforplot <- specdatalong[specdatalong$variable.type == "recruit" &
                              specdatalong$codNames == "NE_Arctic",]
dataforplot$peak <- eigentable[match(dataforplot$codNames, eigentable$codNames),"mode_age"]
dataforplot$maxage <- eigentable[match(dataforplot$codNames, eigentable$codNames),"max_ages"]
dataforplot$sd_mode <- eigentable[match(dataforplot$codNames, eigentable$codNames),"sd_mode"]
dataforplot$cvs_mode <- eigentable[match(dataforplot$codNames, eigentable$codNames),"cvs_mode"]
j <- ggplot(dataforplot, aes(x=freq,y=value,group=codNames)) + 
  geom_line(aes(color=peak)) + theme_classic() + 
  ggtitle("All pops: recruit spectra (color = peak spawning age)") +
  ylim(0,250000)

jj <- ggplot(dataforplot, aes(x=freq,y=value,group=codNames)) + 
  geom_line(aes(color=maxage)) + theme_classic() + 
  ggtitle("All pops: recruit spectra (color = max age)") +
  ylim(0,250000)

jjj <- ggplot(dataforplot, aes(x=freq,y=value,group=codNames)) + 
  geom_line(aes(color=cvs_mode)) + theme_classic() + 
  ggtitle("All pops: recruit spectra (color = CV of spawning distribution)") +
  ylim(0,250000)

plotj <- list(j,jj,jjj)
do.call(grid.arrange,c(plotj,ncol=1))

# *************************************** #
# (10) Calculate & plot k values
# *************************************** #
# In my analysis I want to evaluate population spectra when
# populations sit on different slopes of the BH curve. 
# To change the position on the curve I change beta while keeping
# alpha constant. As beta increases, the slope at equilibrium
# increases (the population moves down the curve). 
# Since LEP is 1 for all populations, equilibrium occurs where
# recruits = egg production (where the 1:1 line crosses the BH.

# For each beta value, plot the BH curve and the 1:1 line
BH <- function(alpha,beta,E){ E/((1/alpha)+(E/beta))}
betas <- c(10^7,10^8,10^9,10^10)
ee <- seq(from=0,to=max(betas),by=100) #range of egg values
a = 1.2

Rlist <- as.list(rep(NA,length=length(betas)))
names(Rlist) <- betas
for(b in 1:length(betas)){ Rlist[[b]] <- BH(alpha=a,beta=betas[b],E=ee) }
RvE <- as.data.frame(do.call(cbind,Rlist))
RvE$eggs <- ee
RvElong <- RvE %>% gather(betavalue,value,1:length(betas))
# plot BH curves
ggplot(RvElong,aes(x=eggs,y=value,color=betavalue)) + geom_line() + 
  geom_abline(intercept=0, slope=1, color="black") #+ xlim(0,1000) + ylim(0,1000)

# For each beta, calculate the intersection with the 1:1 line. Find when R=E
eqpoint <- rep(NA,length=length(betas))
slopes <- rep(NA,length=length(betas))
for(b in 1:length(betas)){
  ee <- seq(from=0,to=max(betas),by=100)
  rr <- BH(alpha=a,beta=betas[b],E=ee)
  dat <- as.data.frame(cbind(ee,rr))
  names(dat) <- c('eggs','recruits')
  dat$diff <- dat$eggs - dat$recruits #diff between R and E
  eqpoint[b] <- which(dat$diff == 0 )
  # calculate slope around equal point
  x1 <- dat[(eqpoint[b]-2),]$eggs
  y1 <- dat[(eqpoint[b]-2),]$recruits
  x2 <- dat[(eqpoint[b]+2),]$eggs
  y2 <- dat[(eqpoint[b]+2),]$recruits
  slopes[b] <- (y2-y1)/(x2-x1)
}
intersets <- cbind(betas,slopes)
plot(x=betas,y=slopes,type="l")
dat[497:503,]
# *************************************** #
# (10) Area under curve
# *************************************** #
head(specdatalong)
# Plan:
# 1. calculate AUC at high and low frequencies (threshold=peak spawning age, others)
# 2. plot %AUClow and %AUChigh vs peak spawning age 
# 3. plot %AUClow and %AUChigh vs max age
# 4. plot %AUClow and %AUChigh vs spawning distribution sd & CV
AUC_greater <- rep(NA,length=length(codNames_ordered_by_peak))
AUC_less <- rep(NA,length=length(codNames_ordered_by_peak))
AUC_total <- rep(NA,length=length(codNames_ordered_by_peak))
AUCperlow <- rep(NA,length=length(codNames_ordered_by_peak))
AUCperhigh <- rep(NA,length=length(codNames_ordered_by_peak))

AUCthreshold <- c(0.05,0.1,0.2,1)
plottitles <- c("threshold freq=0.5","threshold freq=0.1","threshold freq=0.2","threshold freq=peak spawning age")

par(mfrow=c(4,3))
for (j in 1:length(AUCthreshold)){ #step through different AUC thresholds
  
for (i in 1:length(codNames_ordered_by_peak)){
  if(AUCthreshold[j] == 1) 
  {AUCthreshold[j] <- (1/eigentable[eigentable$codNames == codNames_ordered_by_peak[i],]$mode_age)} else
  {AUCthreshold[j]}
  
  AUC_greater[i] <- sum(freq[1]*specdatalong[specdatalong$variable.type == "recruit" 
                                     & specdatalong$codNames==codNames_ordered_by_peak[i]
                                     & specdatalong$freq > AUCthreshold[j],]$value)
  
  AUC_less[i] <- sum(freq[1]*specdatalong[specdatalong$variable.type == "recruit" 
                                  & specdatalong$codNames==codNames_ordered_by_peak[i]
                                  & specdatalong$freq <= AUCthreshold[j],]$value)
  if(j==3){AUC_less[i] <- AUC_less[i]*AUCthreshold[j]} else {AUC_less[i]}
  if(j==3){AUC_greater[i] <- AUC_greater[i]*(0.5-AUCthreshold[j])} else {AUC_greater[i]}
  
  AUC_total[i] <- sum(freq[1]*specdatalong[specdatalong$variable.type == "recruit" 
                                   & specdatalong$codNames==codNames_ordered_by_peak[i],]$value)
  AUCperlow[i] <- AUC_less[i]/AUC_total[i]
  AUCperhigh[i] <- AUC_greater[i]/AUC_total[i]
  print(AUCperlow[i]+AUCperhigh[i])
}
plot(x=eigentable[order(eigentable$mode_age),]$mode, 
               xlab="peak spawning age",y=AUCperlow,ylab="Percent AUC at low freq",main=plottitles[j])
plot(x=eigentable[order(eigentable$mode_age),]$max_ages, 
               xlab="max age",y=AUCperlow,ylab="Percent AUC at low freq",main=plottitles[j])
plot(x=eigentable[order(eigentable$mode_age),]$cvs_mode, 
              xlab="Spawning biomass distribution CV",y=AUCperlow,
             ylab="Percent AUC at low freq",main=plottitles[j])
#plist <- as.list(p1peak,p2max,p3CV)
#pList[[j]] <- plist

}
par(mfrow=c(1,1))

par(mfrow=c(2,2))
plot(x=eigentable[order(eigentable$mode_age),]$mode, 
     xlab="peak spawning age",y=AUC_total,ylab="Total AUC")
plot(x=eigentable[order(eigentable$mode_age),]$max_ages, 
     xlab="max age",y=AUC_total,ylab="Total AUC")
plot(x=eigentable[order(eigentable$mode_age),]$cvs_mode, 
     xlab="Spawning biomass distribution CV",y=AUC_total,
     ylab="Total AUC")
par(mfrow=c(1,1))

# *************************************** #
# (10b) Area under curve - remove pops with truncated distributions
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

