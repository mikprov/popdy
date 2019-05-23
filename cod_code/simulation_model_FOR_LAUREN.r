# Run the simulation model
# by: mikaela provost
# last edited: Feb 19, 2019
# ===================================================================

library(tidyr)
library(dplyr)
library(gridExtra)
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
codNames <- c("Northsea","Coas","W_Baltic",
              "Faroe","NE_Arctic","Celtic",
              "Iceland","Kat","W_Scotland",
              "NGulf","GB","GM",
              "cod3NO","cod3M","cod2J3KL",
              "cod3Ps")



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

# **************** #
# Test Drive Celtic
# **************** #
# modify the Celtic Leslie matrix so that peak age is 5, 10
cel <- Aarray$Celtic[,,1]
celtop <- cel[1,]
addages <- c(2,7,22)
celtopplus <- c(rep(0,length=addages[3]),celtop) 
celLeslie <- matrix(0,nrow=length(celtopplus),ncol=length(celtopplus))
celLeslie[1,] <- celtopplus
for(i in 1:(length(celLeslie[,1])-1)){
  celLeslie[1+i,i] <- 1
}
celarray <- array(NA,dim=c(length(celLeslie[1,]),length(celLeslie[1,]),1))
celarray[,,1] <- celLeslie
# replace modified Leslie into Aarray
Aarray$Celtic <- celarray
rm(cel,celtop,addages,celtopplus,celLeslie,celarray) #clean up

# set params for simulation:
timesteps = 1000 #need this now to create
rm_first_timesteps = 200
alpha = 1.2
beta = 10000 # -----> testing different beta values, re-run script from here
initial_eggs = beta
sig_r = 0.3
span.multiplier = 1 # adjusting the span in spec.prgm()
alphas <- rep(alpha, length=length(codNames)) #alpha could be diff for pops
output.3d.list <- as.list(rep(NA,length=length(codNames))) #store timeseries here
names(output.3d.list) <- codNames

for (i in 1:length(Aarray)) { #step through each pop
  Leslie3d = Aarray[[i]] #select the 3d array of Leslie matricies
  # array dims: row=ts length, col=4 is number of ts (eggs,recruits,Nt,Nsize), depth=F vals
  output.matrix <- array(NA,c(timesteps-2,4,length(Fvalues))) 
  
  for (f in 1:length(Fvalues)) { #step through each Leslie matrix (for each F value)
    output = sim_model(A=Leslie3d[,,f], timesteps=timesteps, 
                       alpha=alphas[i], beta=beta, 
                       sig_r=sig_r, initial_eggs=initial_eggs)
    
    length(output$Nsize) <- length(output$N_t) #trim Nsize ts vector, -2 elements
    output.matrix[,,f] <- do.call(cbind,output) #fill in array for pop i
    #colnames(output.matrix) <- names(output)
  }
  output.3d.list[[i]] <- output.matrix
}
rm(i,f,Leslie3d,output.matrix,output) #clean up

# (3) Format output ts for plotting simulations using output.3d.list
variable_type <- c("Nt","eggs","recruits","Nsize")
# --- reorganize egg timeseries data --- #
var.number <- 2 # eggs
df.list <- as.list(rep(NA,length=length(codNames)))
names(df.list) <- codNames
for (i in 1:length(output.3d.list)) {
  # first, reformat data to work with ggplot
  aa <- as.data.frame(output.3d.list[[i]][,var.number,])
  aa$year <- seq(from=1, to=length(aa[,1]))
  colnames(aa) <- c(Fvalues,"year")
  aa1 <- aa %>% gather(Fval,value,1:length(Fvalues))
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
  colnames(aa) <- c(Fvalues,"year")
  aa1 <- aa %>% gather(Fval,value,1:length(Fvalues))
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
  colnames(aa) <- c(Fvalues,"year")
  aa1 <- aa %>% gather(Fval,value,1:length(Fvalues))
  aa1$variable <- rep(variable_type[var.number],length=length(aa1[,1]))
  aa1$codNames <- rep(codNames[i],length=length(aa[,1]))
  df.list[[i]] <- aa1
  rm(aa1,aa)}
nsize.ts <- do.call(rbind,df.list)
rm(df.list) #clean up

ts.data <- rbind(eggs.ts,recruits.ts,nsize.ts) #combine data
rownames(ts.data) <- NULL

# Plan:
# 1. Walk through each cod pop
# 2. Store spec values for eggs, recruits, and Nsize

# 1. Walk through each cod pop
sp.eggsL <- as.list(rep(NA,length=length(codNames))) #object for spec analysis  
sp.recruitL <- as.list(rep(NA,length=length(codNames))) 
ts.for.spec.eg <- as.list(rep(NA,length=length(codNames)))
ts.for.spec.re <- as.list(rep(NA,length=length(codNames)))

for (i in 1:length(codNames)){
  
  ts <- ts.data[ts.data$codNames == codNames[i],] #subset data for pop i 
  
  # setting 'span' - a vector of odd integers to specify the smoothers
  tmp <- ceiling(sqrt(length(1:(timesteps-rm_first_timesteps-1)))) #sq root of timeseries lgth, rounded
  if (tmp %% 2 == 0) {m <- tmp+1} else {m <- tmp} #make it odd, if the square root is even
  m = m * span.multiplier
  
  xx = ts[ts$variable == "eggs",]$value[rm_first_timesteps:(timesteps-2)] - mean(ts[ts$variable == "eggs",]$value[rm_first_timesteps:(timesteps-2)])
  ts.for.spec.eg[[i]] <- xx #save time series that goes into spec.pgram()
  sp = spec.pgram(x=xx,spans=c(m,m),plot = FALSE)
  #save spec output for plotting pops together, Helen Wearing says to multiply by 2
  sp.eggsL[[i]] <- sp$spec*2
  rm(xx,sp)
  
  yy = ts[ts$variable == "recruits",]$value[rm_first_timesteps:(timesteps-2)] - mean(ts[ts$variable == "recruits",]$value[rm_first_timesteps:(timesteps-2)])
  ts.for.spec.re[[i]] <- yy #save time series that goes into spec.prgam()
  sp = spec.pgram(yy,spans=c(m,m),plot = FALSE)
  sp.recruitL[[i]] = sp$spec*2 # save matrix of spec values for different FLEP, index by pop i
  rm(yy)
  
}
freq <- sp$freq
sp.eggs <- as.data.frame(do.call(cbind,sp.eggsL))
sp.recruit <- as.data.frame(do.call(cbind,sp.recruitL))
names(sp.eggs) <- codNames
names(sp.recruit) <- codNames
sp.eggs$freq <- freq
sp.recruit$freq <- freq

# Plan
# 1. re-arrange dataframes for plotting (sp.eggs, sp.recruit)
# 2. combine 2 variable type data frames to make one df 
# 3. plot

# 1. re-arrange egg & recruit dataframes
egglong <- sp.eggs %>% gather("codNames","value",1:length(codNames))
egglong$variable.type <- rep("eggs",length=length(egglong$freq))
head(egglong)
recruitlong <- sp.recruit %>% gather("codNames","value",1:length(codNames))
recruitlong$variable.type <- rep("recruit",length=length(recruitlong$freq))
head(recruitlong)
# 2. combine 2 variable type data frames to make one df 
specdatalong <- rbind(egglong,recruitlong)
head(specdatalong)

# 3. save specdatalong for each variation of Celtic Leslie matrix
specdatalong3 
specdatalong5 
specdatalong10 
specdatalong25 
# add peak age column
specdatalong3$peak <- rep(3,length=length(specdatalong3$freq))
specdatalong5$peak <- rep(5,length=length(specdatalong5$freq))
specdatalong10$peak <- rep(10,length=length(specdatalong10$freq))
specdatalong25$peak <- rep(25,length=length(specdatalong25$freq))


# combine specdatalong dfs
specdatalong <- rbind(specdatalong3,specdatalong5,specdatalong10,specdatalong25)
specdatalong$peak <- factor(specdatalong$peak)
levels(specdatalong$peak)

specdatalong3510_x2 
specdatalong3510 
specdatalong3510_span3 <- specdatalong
specdatalong351025 <- specdatalong

specdatalong <- specdatalong351025

# probability of spawning
source(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_pops/',"Celtic", '.r', sep=''))
# calculate LEP at each age
lsb.at.k = calculate_LSB_at_age_by_F(maxage=maxage,L_inf=L_inf,K=K,TEMP=TEMP,
                                     F.halfmax=0,B0=B0,B1=B1)
Ages = seq(from=1,to=length(lsb.at.k[,1]),by=1)
# calculate probability of spawning at age
p_spawn = as.data.frame(lsb.at.k[,1] / sum(lsb.at.k[,1])) 
colnames(p_spawn) <- "p_spawn"
p_spawnT= cbind(p_spawn,Ages)
# additional rows for adjusting peak age
addspeak5 <- as.data.frame(cbind(rep(0,length=2),rep(0,length=2)))
addspeak10 <- as.data.frame(cbind(rep(0,length=7),rep(0,length=7)))
addspeak25 <- as.data.frame(cbind(rep(0,length=22),rep(0,length=22)))
names(addspeak5) <- c("p_spawn","Ages")
names(addspeak10) <- c("p_spawn","Ages")
names(addspeak25) <- c("p_spawn","Ages")
# create dfs for 3,5,10 -- for plotting
p_spawn3 <- p_spawnT
p_spawn3$peak <- rep(3,length=length(p_spawn3$p_spawn))
p_spawn5 <- rbind(addspeak5,p_spawnT)
p_spawn5$Ages <- seq(from=1,to=length(p_spawn5$p_spawn),by=1)
p_spawn5$peak <- rep(5,length=length(p_spawn5$p_spawn))

p_spawn10 <- rbind(addspeak10,p_spawnT)
p_spawn10$Ages <- seq(from=1,to=length(p_spawn10$p_spawn),by=1)
p_spawn10$peak <- rep(10,length=length(p_spawn10$p_spawn))

p_spawn25 <- rbind(addspeak25,p_spawnT)
p_spawn25$Ages <- seq(from=1,to=length(p_spawn25$p_spawn),by=1)
p_spawn25$peak <- rep(25,length=length(p_spawn25$p_spawn))

p_spawnplot <- rbind(p_spawn3,p_spawn5,p_spawn10,p_spawn25)
p_spawnplot$peak <- factor(p_spawnplot$peak)
spawnplot <- ggplot(data=p_spawnplot,aes(x=Ages,y=p_spawn,color=peak)) +
  geom_line() + ggtitle("Celtic probability of spawning at age")

# 3. plot frequency content of eggs

# --- recruits spectra: all on one plot ---#
dataforplot <- specdatalong[specdatalong$variable.type == "recruit" &
                              specdatalong$codNames == "Celtic",]
c <- ggplot(dataforplot, aes(x=freq,y=value,color=peak)) +
  geom_line() + geom_vline(xintercept = c(1/3,1/5,1/10),linetype="dashed") +
  ggtitle("Celtic test drive: peak age at 3, 5, 10") + ylab("value")

clog <- ggplot(dataforplot, aes(x=freq,y=log10(value),color=peak)) +
  geom_line() + geom_vline(xintercept = c(1/3,1/5,1/10),linetype="dashed") +
  ggtitle("Celtic test drive: peak age at 3, 5, 10") + ylab("log(value)")

clist <- list(c,clog,spawnplot)
do.call(grid.arrange,c(clist,ncol=1))




# *************************************** #
# (2) Simulate pops. Loop over Aarray list to simulate using different Leslie matrices 
# *************************************** #
# set params for simulation:
timesteps = 1000 #need this now to create
rm_first_timesteps = 200
alpha = 1.2
beta = 10000 # -----> testing different beta values, re-run script from here
initial_eggs = beta
sig_r = 0.3
span.multiplier = 1 # adjusting the span in spec.prgm()
alphas <- rep(alpha, length=length(codNames)) #alpha could be diff for pops
#alphas <- c(2.02,8.38,4.1,2.10,3.77,2.10,3.24,4.3,5.54,0.48,1.69,0.77,2.31,0.42,1.22)

output.3d.list <- as.list(rep(NA,length=length(codNames))) #store timeseries here
names(output.3d.list) <- codNames

for (i in 1:length(Aarray)) { #step through each pop
  Leslie3d = Aarray[[i]] #select the 3d array of Leslie matricies
  # array dims: row=ts length, col=4 is number of ts (eggs,recruits,Nt,Nsize), depth=F vals
  output.matrix <- array(NA,c(timesteps-2,4,length(Fvalues))) 
  
  for (f in 1:length(Fvalues)) { #step through each Leslie matrix (for each F value)
    output = sim_model(A=Leslie3d[,,f], timesteps=timesteps, 
                       alpha=alphas[i], beta=beta, 
                       sig_r=sig_r, initial_eggs=initial_eggs)
    
    length(output$Nsize) <- length(output$N_t) #trim Nsize ts vector, -2 elements
    output.matrix[,,f] <- do.call(cbind,output) #fill in array for pop i
    #colnames(output.matrix) <- names(output)
  }
  output.3d.list[[i]] <- output.matrix
}
rm(i,f,Leslie3d,output.matrix,output) #clean up

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
  colnames(aa) <- c(Fvalues,"year")
  aa1 <- aa %>% gather(Fval,value,1:length(Fvalues))
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
  colnames(aa) <- c(Fvalues,"year")
  aa1 <- aa %>% gather(Fval,value,1:length(Fvalues))
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
  colnames(aa) <- c(Fvalues,"year")
  aa1 <- aa %>% gather(Fval,value,1:length(Fvalues))
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
sd_ordered <- rep(NA,length=length(codNames_ordered_by_peak))
mean_ordered <- rep(NA,length=length(codNames_ordered_by_peak))
CV_ordered <- rep(NA,length=length(codNames_ordered_by_peak))

for (i in 1:length(codNames_ordered_by_peak)){
  dd <- ts.data[ts.data$variable == "recruits" & 
                  ts.data$codNames == codNames_ordered_by_peak[i] &
                  ts.data$year %in% seq(from=rm_first_timesteps,to=(timesteps-2),by=1),]
  
  prec[[i]] <- ggplot(dd, #aes(x=year,y=value,color=Fval)) +
                   aes(x=year,y=value)) +
    xlab("year") + ylab("recruits (before noise)") +
    geom_line() + theme_classic() + #ylim(c(1000,3000)) +
    #scale_color_brewer(palette = "Reds") +
    ggtitle(paste(codNames_ordered_by_peak[i],"var=",
                  round(var(dd$value),digits=1)))
  
  sd_ordered[i] <- sd(dd$value)
  mean_ordered[i] <- mean(dd$value)
  CV_ordered[i] <- sd_ordered[i]/mean_ordered[i]
}
#pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/recruitment_beforenoise_for_diff_FLEPs_alpha50_meanvsSTDEV.pdf', width=7, height=10) #note: file name specifies the alpha used in simluation model
do.call(grid.arrange,c(prec,ncol=3))
#dev.off()

plot(x=sd_ordered^2, y=AUC_total)
abline(a=0,b=1)


# plot eggs - one plot per pop, on each plot 5 lines for different F levels
peggs <- list()
# reorder codNames by peak spawning age (increasing)
codNames_ordered_by_peak <- eigentable %>% arrange(mode_age) %>% pull(codNames)
sd_ordered <- rep(NA,length=length(codNames_ordered_by_peak))
mean_ordered <- rep(NA,length=length(codNames_ordered_by_peak))
CV_ordered <- rep(NA,length=length(codNames_ordered_by_peak))

for (i in 1:length(codNames_ordered_by_peak)){
  
  dd <- ts.data[ts.data$variable == "eggs" & 
                  ts.data$codNames == codNames_ordered_by_peak[i] &
                  ts.data$year %in% seq(from=rm_first_timesteps,to=(timesteps-2),by=1),]
  
  peggs[[i]] <- ggplot(dd, #aes(x=year,y=value,color=Fval)) +
                   aes(x=year,y=value)) +
    xlab("year") + ylab("egg production") +
    geom_line() + theme_classic() + #ylim(c(0,600)) +
    #scale_color_brewer(palette = "Reds") +
    ggtitle(paste(codNames_ordered_by_peak[i],"var=",
                  round(var(dd$value),digits=0)))
    
  sd_ordered[i] <- sd(dd$value)
  mean_ordered[i] <- mean(dd$value)
  CV_ordered[i] <- sd_ordered[i]/mean_ordered[i]
}
#pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/egg_production_for_diff_FLEPs_alpha50_meanvsStdev.pdf', width=7, height=10) #note: file name specifies the alpha used in simluation model
do.call(grid.arrange,c(peggs,ncol=3))


# *************************************** #
# (5) Calculate frequency content from timeseries
# *************************************** #
# Plan:
# 1. Walk through each cod pop, do spectral analysis at F levels
# 2. Store spec values for eggs, recruits, and Nsize

# 1. Walk through each cod pop, do spectral analysis at F levels
sp.eggsL <- as.list(rep(NA,length=length(codNames))) #object for spec analysis  
sp.recruitL <- as.list(rep(NA,length=length(codNames))) 
ts.for.spec.eg <- as.list(rep(NA,length=length(codNames)))
ts.for.spec.re <- as.list(rep(NA,length=length(codNames)))

for (i in 1:length(codNames)){
  
  ts <- ts.data[ts.data$codNames == codNames[i],] #subset data for pop i 
  
  # setting 'span' - a vector of odd integers to specify the smoothers
  tmp <- ceiling(sqrt(length(1:(timesteps-rm_first_timesteps-1)))) #sq root of timeseries lgth, rounded
  if (tmp %% 2 == 0) {m <- tmp+1} else {m <- tmp} #make it odd, if the square root is even
  m = m * span.multiplier
  
  # --- spectral analysis on EGGS --- #
  # Problem I ran into: spectra shape is the same for different beta values
  # When I change beta, I am chaning 1) equilibirum and 2) position
  # Divide time series by the equilibrium value -- isolate just position
  #equilibrium_egg = mean(ts[ts$variable == "eggs",]$value[rm_first_timesteps:(timesteps-2)])
  #xxx = ts[ts$variable == "eggs",]$value[rm_first_timesteps:(timesteps-2)] / equilibrium_egg
  #xx = xxx-mean(xxx) #after dividing by equilibrium, subtract mean from time series 
  
  xx = ts[ts$variable == "eggs",]$value[rm_first_timesteps:(timesteps-2)] - mean(ts[ts$variable == "eggs",]$value[rm_first_timesteps:(timesteps-2)])
  ts.for.spec.eg[[i]] <- xx #save time series that goes into spec.pgram()
  sp = spectrum(x=xx,spans=c(m,m),plot = FALSE)
  #save spec output for plotting pops together, Helen Wearing says to multiply by 2
  sp.eggsL[[i]] <- sp$spec*2
  rm(xx,sp)
  
  # --- spectral analysis on RECRUIT --- #
  # Follow same method as egg production
  # Divide time series by the equilibrium value
  #equilibrium_rec <- mean(ts[ts$variable == "recruits",]$value[rm_first_timesteps:(timesteps-2)])
  #yyy = ts[ts$variable == "recruits",]$value[rm_first_timesteps:(timesteps-2)] / equilibrium_rec
  #yy = yyy - mean(yyy) #after dividing by equilibrium value -- isolate just position
  
  yy = ts[ts$variable == "recruits",]$value[rm_first_timesteps:(timesteps-2)] - mean(ts[ts$variable == "recruits",]$value[rm_first_timesteps:(timesteps-2)])
  ts.for.spec.re[[i]] <- yy #save time series that goes into spec.prgam()
  sp = spec.pgram(yy,spans=c(m,m),plot = FALSE)
  sp.recruitL[[i]] = 2*sp$spec # save matrix of spec values for different FLEP, index by pop i
  rm(yy)
  
}
freq <- sp$freq
sp.eggs <- as.data.frame(do.call(cbind,sp.eggsL))
sp.recruit <- as.data.frame(do.call(cbind,sp.recruitL))
names(sp.eggs) <- codNames
names(sp.recruit) <- codNames
sp.eggs$freq <- freq
sp.recruit$freq <- freq

# format time series for spec.prgam() for plotting
ts.for.spec.eg.df <- as.data.frame(do.call(cbind,ts.for.spec.eg))
names(ts.for.spec.eg.df) <- codNames
ts.for.spec.eg.df$year <- seq(from=1,to=length(ts.for.spec.eg.df[,1]),by=1)
ts.for.spec.eg.long <- ts.for.spec.eg.df %>% gather("codNames","value",1:length(codNames))
ts.for.spec.eg.long$variable.type <- rep("eggs",length=length(ts.for.spec.eg.long$codNames))

ts.for.spec.re.df <- as.data.frame(do.call(cbind,ts.for.spec.re))
names(ts.for.spec.re.df) <- codNames
ts.for.spec.re.df$year <- seq(from=1,to=length(ts.for.spec.re.df[,1]),by=1)
ts.for.spec.re.long <- ts.for.spec.re.df %>% gather("codNames","value",1:length(codNames))
ts.for.spec.re.long$variable.type <- rep("recruit",length=length(ts.for.spec.re.long$codNames))

ts.for.spec <- rbind(ts.for.spec.eg.long, ts.for.spec.re.long)

# plot time series that was input into spec.prgam()
p <- list()
for (i in 1:length(codNames_ordered_by_peak)){
  
  dataforplot <- ts.for.spec[ts.for.spec$variable.type == "eggs" 
                              & ts.for.spec$codNames==codNames_ordered_by_peak[i],]
  # store plots in list
  p[[i]] <- ggplot(data=dataforplot, aes(x=year,y=value)) + 
    geom_line() + #ylim(1,10) +
    ylab("eggs") +
    ggtitle(paste(codNames_ordered_by_peak[i])) + 
    theme_classic()
}
do.call(grid.arrange,c(p,ncol=3))
rm(dataforplot)

# *************************************** #
# (6) Plot spectral analysis 
# *************************************** #

# Plan
# 1. re-arrange dataframes for plotting (sp.eggs, sp.recruit)
# 2. combine 2 variable type data frames to make one df 
# 3. plot

# 1. re-arrange egg dataframe
egglong <- sp.eggs %>% gather("codNames","value",1:length(codNames))
egglong$variable.type <- rep("eggs",length=length(egglong$freq))
head(egglong)
recruitlong <- sp.recruit %>% gather("codNames","value",1:length(codNames))
recruitlong$variable.type <- rep("recruit",length=length(recruitlong$freq))
head(recruitlong)
# 2. combine 2 variable type data frames to make one df 
specdatalong <- rbind(egglong,recruitlong)
head(specdatalong)
# 3. plot frequency content of eggs
# here I can change the order I want plots to appear:
codNames_ordered_by_peak <- eigentable %>% arrange(mode_age) %>% pull(codNames)

# --- egg spectra --- # 

peggs_sp <- list()
for (i in 1:length(codNames_ordered_by_peak)){
  
  dataforplot <- specdatalong[specdatalong$variable.type == "eggs" 
                              & specdatalong$codNames==codNames_ordered_by_peak[i],]
  # store plots in list
  peggs_sp[[i]] <- ggplot(data=dataforplot, aes(x=freq,y=log10(value))) + 
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
  dataforplot <- specdatalong[specdatalong$variable.type == "recruit" 
                              & specdatalong$codNames==codNames_ordered_by_peak[i],]
  # store plots in list
  prec_sp[[i]] <- ggplot(data=dataforplot,aes(x=freq,y=value)) + 
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

