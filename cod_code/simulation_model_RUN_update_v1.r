# Run the simulation model
# by: mikaela provost
# last edited: Feb 19, 2019
# ===================================================================

library(tidyr)
library(dplyr)
# ---
# load the simulation model
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/simulation_model_cod_v3.r")

# load functions
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/2_cod_functions.r")

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

# ---
# create empty objects for simulation
freq = seq(from=0.001111111, to=0.5, by=0.001111111) #all frequencies, for plotting later
spall = matrix(NA, length(freq),length(eigentable$codNames)) #set up empty matrix, spec for all plots
outputL = vector("list", length(eigentable$codNames)) #create empty list to store sim output
names(outputL) = eigentable$codNames
# ---




# *************************************** #
# (1) calc LEP from range of F values (create LEP matrix)
# *************************************** #
Fvalues = seq(0,3,by=0.01) #Check max F value (some pops can withstand high F)
FLEP = matrix(NA,nrow=length(Fvalues),ncol=length(codNames)) #store FLEP values here
colnames(FLEP) = codNames
LEPmatrix = matrix(NA,nrow=length(Fvalues),ncol=length(codNames))
for (i in 1:length(codNames)) { #step through each cod population
  
  # load parms for cod pop i: L_inf, K (for vonB), TEMP, maxage
  source(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_pops/',codNames[i], '.r', sep=''))
  
  #for (f in 1:length(Fvalues)){ #step through each F value
  # calculate_LSB_at_age_by_F() creates a matrix: col=F vals, row=age
  LSB_at_age = calculate_LSB_at_age_by_F(maxage=maxage, 
                                             L_inf=L_inf, 
                                             K=K, 
                                             TEMP=TEMP, 
                                             F.halfmax=Fvalues, #here are f values
                                             B0=B0,B1=B1)
  LEPmatrix[,i] <- colSums(LSB_at_age) #store LEP at different F values for pop i
}


# *************************************** #
# (2) calc FLEP from LEPmatrix (create FLEP matrix)
# *************************************** #
# here I calculate FLEP: LSB at all F levels / LSB at F=0
FLEP = t(t(LEPmatrix)/LEPmatrix[1,])
colnames(FLEP) <- codNames #assign pop names
#FLEP <- round(FLEP, digits=2)
FLEP <- as.data.frame(cbind(Fvalues,FLEP)) #first col=F values, rest are FLEP vals for pops

# [make plot: F vs FLEP]
FLEPplot <- FLEP %>% gather(pop,value,2:17)
FLEPvsF <- ggplot(FLEPplot,aes(x=Fvalues,y=log(value),color=pop)) +
  geom_line() + theme_classic() +
  xlab("F") + ylab("ln(FLEP)")


# *************************************** #
# (3) Generate Leslie matricies for diff F values (create Leslie arrays)
# *************************************** #
Aarray = as.list(rep(NA,length(codNames))) #Leslie matrix storage for each F value for pop i
eigenvals1 = matrix(NA,nrow=length(FLEP$Fvalues),ncol=length(codNames)) 
eigenvals2 = matrix(NA,nrow=length(FLEP$Fvalues),ncol=length(codNames))
eigenvals12 = matrix(NA,nrow=length(FLEP$Fvalues),ncol=length(codNames))


for (i in 1:length(codNames)){ #for each pop i
  # load parms for cod pop i: L_inf, K (for vonB), TEMP, maxage,B0,B1 (matur)
  source(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_pops/',codNames[i], '.r', sep=''))
  
  Lesliearray <- array(NA,c(maxage,maxage,length(Fvalues))) #store Leslie matricies
  e1 = rep(NA,length=length(Fvalues)) #store lambda1
  e2 = rep(NA,length=length(Fvalues)) #store lambda2
  e12 = rep(NA,length=length(Fvalues)) #store inverse damping ratio
  
  for (f in 1:length(FLEP$Fvalues)){ #step through F values 
    # create Leslie matrix:
    Leslieout = assemble_Leslie(maxage=maxage, K=K, L_inf=L_inf, TEMP=TEMP,
                                F.halfmax=FLEP$Fvalues[f], tknot=0, B0=B0, B1=B1)
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
rm(e1,e2,e12,i,f,Leslieout,Lesliearray) #clean up


# *************************************** #
# (4) Simulate pops. Loop over Aarray list to simulate using different Leslie matrices 
# *************************************** #
# set params for simulation:
timesteps = 200 #need this now to create
rm_first_timesteps = 50
beta = 100 
initial_eggs = 100
sig_r = 0.7
span.multiplier = 1 # what is this again?
alphas <- rep(50, length=length(codNames)) #alpha could be diff for pops

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
# (5) Format output ts for plotting simulations using output.3d.list
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

# add FLEP levels from FLEP df to ts.data (matching on Fvalues)
final <- matrix(NA,nrow = 0, ncol = (ncol(ts.data) + 1))
final <- as.data.frame(final)
names(final) <- c(names(ts.data),'FLEP')

for (i in 1:length(codNames)){ 
  print(i)
  sub <- ts.data[ts.data$codNames == codNames[i],]
  sub2 <- FLEP[,c('Fvalues', codNames[i])]
  merged <- merge(x = sub, y = sub2, by.x = 'Fval', by.y = 'Fvalues', sort = F)
  names(merged)[ncol(merged)] <- 'FLEP'
  final <- rbind(final, merged)
}
rm(i,sub,sub2,merged)
head(final)

# *************************************** #
# (6) Subset final to have F values associated with tarvel FLEP levels 
# *************************************** #
# find the nearest F values that correspond to the desired FLEP value
target_Fs <- rbind(
  
  as.data.frame(melt(FLEP,id="Fvalues") %>% 
                  group_by(variable) %>%
                  arrange(abs(value-1)) %>% #desired FLEP=1
                  slice(1) %>%
                  #mutate(value=round(value,2)) %>%
                  mutate(FLEPlevel=rep(1))),
  
  as.data.frame(melt(FLEP,id="Fvalues") %>% 
                  group_by(variable) %>%
                  arrange(abs(value-0.8)) %>% #desired FLEP=0.8
                  slice(1) %>%
                  #mutate(value=round(value,2)) %>%
                  mutate(FLEPlevel=rep(0.8))),
  
  as.data.frame(melt(FLEP,id="Fvalues") %>% 
                  group_by(variable) %>%
                  arrange(abs(value-0.5)) %>% #desired FLEP=0.5
                  slice(1) %>%
                  #mutate(value=round(value,2)) %>%
                  mutate(FLEPlevel=rep(0.5))),
  
  as.data.frame(melt(FLEP,id="Fvalues") %>% 
                  group_by(variable) %>%
                  arrange(abs(value-0.35)) %>% #desired FLEP=0.35
                  slice(1) %>%
                  #mutate(value=round(value,2)) %>%
                  mutate(FLEPlevel=rep(0.35))),
  
  as.data.frame(melt(FLEP,id="Fvalues") %>% 
                  group_by(variable) %>%
                  arrange(abs(value-0.2)) %>% #desired FLEP=0.2
                  slice(1) %>%
                  #mutate(value=round(value,2)) %>%
                  mutate(FLEPlevel=rep(0.2))) 
  
)
#target_FLEP = c(1,0.8,0.5,0.35,0.2)

sublist = as.list(NA, rep(length(codNames)))
for (i in 1:length(codNames)){
  sub <- final[final$codNames == codNames[i],]
  sub2 <- target_Fs[target_Fs$variable == codNames[i],]
  sub3 <- sub[sub$Fval %in% sub2$Fvalues,]
  sublist[[i]] <- sub3
}
finaltargets <- do.call(rbind,sublist)
rm(i,sub,sub2,sub3,sublist) #clean up



# *************************************** #
# (7) Now that timeseries data is formated, let's plot!
# *************************************** #

# plot eggs - one plot per pop, on each plot 5 lines for different F levels
p <- list()
for (i in 1:length(codNames)){
  
  p[[i]] <- ggplot(finaltargets[finaltargets$variable == "eggs" & 
                                  finaltargets$codNames == codNames[i],], 
                   aes(x=year,y=value,color=Fval)) +
    xlab("year") + ylab("egg production") +
    geom_line() +
    scale_color_brewer(palette = "Reds") +
    ggtitle(paste(codNames[i])) +
    theme_classic()
}
pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/egg_production_for_diff_FLEPs_alpha50.pdf', width=7, height=10) #note: file name specifies the alpha used in simluation model
do.call(grid.arrange,c(p,ncol=2))
dev.off()
rm(p,i)



#plot recruitment - one plot per pop, similar to egg plots
p <- list()
for (i in 1:length(codNames)){
  p[[i]] <- ggplot(finaltargets[finaltargets$variable == "recruits" & 
                                  finaltargets$codNames == codNames[i],], 
                   aes(x=year,y=value,color=Fval)) + 
    xlab("year") + ylab("recruits") +
    geom_line() + scale_color_brewer(palette = "Reds") +
    ggtitle(paste(codNames[i])) +
    theme_classic()
}
pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/recruitment_beforenoise_for_diff_FLEPs_alpha50.pdf', width=7, height=10) #note: file name specifies the alpha used in simluation model
do.call(grid.arrange,c(p,ncol=2))
dev.off()
rm(p,i)

#plot pop size - one plot per pop, similar to egg & recruitment plots
p <- list()
for (i in 1:length(codNames)){
  p[[i]] <- ggplot(finaltargets[finaltargets$variable == "Nsize" & 
                                  finaltargets$codNames == codNames[i],], 
                   aes(x=year,y=value,color=Fval)) + 
    xlab("year") + ylab("Nsize") +
    geom_line() + scale_color_brewer(palette = "Reds") +
    ggtitle(paste(codNames[i])) +
    theme_classic()
}
pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/Nsize_for_diff_FLEPs_alpha50.pdf', width=7, height=10) #note: file name specifies the alpha used in simluation model
do.call(grid.arrange,c(p,ncol=2))
dev.off()
rm(p,i)





# ------------------------------------------------------------------------- #
#  Extra code below
# ------------------------------------------------------------------------- #

    
    
  # ---
  # Define alpha for pop i -- this is a big deal!
  alpha <- 1/FLEPinfo[FLEPinfo$FLEPlevel == 0.35,]$LEP
  alphas[i] <- alpha #save alphas for popultions. 
  # --- #
  # make sure enough of a fishing range to get pops down to 0.2 FLEP. Look at: 
#FLEPinfoDF <- do.call("rbind",FLEPinfolist)
#alphatable <- cbind(codNames,as.data.frame(alphas))

# CHECK THIS: is alpha related to high fishing at FLEP=0.2?
#plot(x = FLEPinfoDF[FLEPinfoDF$FLEPlevel == 0.35,]$Fvalues,
#     y = alphas, ylab="alpha = 1/(LEP*0.35)", xlab="F at FLEP=0.2")
#text(alphas~FLEPinfoDF[FLEPinfoDF$FLEPlevel == 0.35,]$Fvalues, 
#     labels=codNames,cex=0.9, font=2, pos=3)
  rm(L_inf,K,TEMP,B0,B1) #rm these variables after each loop, seems safer
 
rm(i,eigenvals1,eigenvals2,eigenvals12,Leslieout,A,FLEPinfo,
   target_Fs,target_LEPs,Alist,LEP,FLEP,FLEP.F,FLEP.LEP,
   MG,Mp,name,rm_first_timesteps,S50,f) #clean up



# ------------------------------------------------------------------------- #
# --- 
# ------------------------------------------------------------------------- #






# ------------------------------------------------------------------------- #
# --- Calculate frequency content from timeseries
# ------------------------------------------------------------------------- #

# From earlier above, I have:
# 1. [output.3d.list] a list of 3d arrays. Each array is timeseries output
#    from simulations at different F levels. 

# Plan:
# 1. Walk through each cod pop, calculate spectral analysis at different F levels
# 2. Store spec values for eggs, recruits, and Nsize at different FLEP levels
dataNAs <- rep(NA,length(freq)*5*length(codNames))
sp.eggs <- array(dataNAs,c(length(freq),5,length(codNames)))# 3D object with spectral analysis for 
sp.recruit <- array(dataNAs,c(length(freq),5,length(codNames)))# 3D object with spectral analysis for 
sp.nsize <- array(dataNAs,c(length(freq),5,length(codNames)))# 3D object with spectral analysis for 

for (i in 1:length(codNames)){

  ts3d <- output.3d.list[[i]] #subset list to data for pop i 
  
  # setting 'span' - a vector of odd integers to specify the smoothers
  tmp <- ceiling(sqrt(length(1:(timesteps-rm_first_timesteps-1)))) #sq root of timeseries lgth, rounded
  if (tmp %% 2 == 0) {m <- tmp+1} else {m <- tmp} #make it odd, if the square root is even
  m = m * span.multiplier
  
  # --- spectral analysis on EGGS --- #
  e = matrix(NA, length(freq),length(ts3d[1,1,])) #store egg sp vals at FLEP levels here
  for (f in 1:length(ts3d[1,1,])){ #step through FLEP levels
    sp = spec.pgram(x=ts3d[,2,f][rm_first_timesteps:(timesteps-2)],spans=c(m,m),plot = FALSE)
    e[,f] = 2*sp$spec} #save spec output for plotting pops together, Helen Wearing says to multiply by 2
  
  # --- spectral analysis on RECRUIT --- #
  r = matrix(NA, length(freq),length(ts3d[1,1,])) #store recruit sp vals at FLEP levels here
  for (f in 1:length(ts3d[1,1,])){ #step through FLEP levels
    sp = spec.pgram(x=ts3d[,3,f][rm_first_timesteps:(timesteps-2)],spans=c(m,m),plot = FALSE)
    r[,f] = 2*sp$spec} #save spec output for plotting pops together, Helen Wearing says to multiply by 2 
  
  # --- spectral analysis on NSIZE --- #
  n = matrix(NA, length(freq),length(ts3d[1,1,])) #store nsize sp vals at FLEP levels here
  for (f in 1:length(ts3d[1,1,])){ #step through FLEP levels
    sp = spec.pgram(x=ts3d[,4,f][rm_first_timesteps:(timesteps-2)],spans=c(m,m),plot = FALSE)
    n[,f] = 2*sp$spec} #save spec output for plotting pops together, Helen Wearing says to multiply by 2 
  
  sp.eggs[,,i] = e # save matrix of spec values for different FLEP, index by pop i
  sp.recruit[,,i] = r # save matrix of spec values for different FLEP, index by pop i
  sp.nsize[,,i] = n # save matrix of spec values for different FLEP, index by pop i
  
}
rm(i,f,sp,e,r,n,tmp,m) #clean up



# ------------------------------------------------------------------------- #
# --- Plot spectral analysis 
# ------------------------------------------------------------------------- #

# Plan
# 1. re-arrange 3d dataframes for plotting (sp.eggs, sp.recruit, sp.nsize)
# 2. combine all 3 variable type data frames to make one df 
# 3. plot

# 1. re-arrange egg dataframe
templist <- as.list(rep(NA,length=length(codNames)))
names(templist) <- codNames
for (i in 1:length(codNames)){
  df <- as.data.frame(sp.eggs[,,i])
  colnames(df) <- c("FLEP1","FLEP0.8","FLEP0.5","FLEP0.35","FLEP0.2")
  df$freq <- freq
  df$variable.type <- rep("eggs",length=length(freq))
  df$codNames <- rep(paste(codNames[i]),length=length(freq))
  templist[[i]] <- df
}
specdata.e <- do.call("rbind",templist)
rownames(specdata.e) <- NULL
rm(i,templist,df)

# re-arrange recruit dataframe
templist <- as.list(rep(NA,length=length(codNames)))
names(templist) <- codNames
for (i in 1:length(codNames)){
  df <- as.data.frame(sp.recruit[,,i])
  colnames(df) <- c("FLEP1","FLEP0.8","FLEP0.5","FLEP0.35","FLEP0.2")
  df$freq <- freq
  df$variable.type <- rep("recruit",length=length(freq))
  df$codNames <- rep(paste(codNames[i]),length=length(freq))
  templist[[i]] <- df
}
specdata.r <- do.call("rbind",templist)
rownames(specdata.r) <- NULL
rm(i,templist,df)

# re-arrange nsize dataframe
templist <- as.list(rep(NA,length=length(codNames)))
names(templist) <- codNames
for (i in 1:length(codNames)){
  df <- as.data.frame(sp.nsize[,,i])
  colnames(df) <- c("FLEP1","FLEP0.8","FLEP0.5","FLEP0.35","FLEP0.2")
  df$freq <- freq
  df$variable.type <- rep("nsize",length=length(freq))
  df$codNames <- rep(paste(codNames[i]),length=length(freq))
  templist[[i]] <- df
}
specdata.n <- do.call("rbind",templist)
rownames(specdata.n) <- NULL
rm(i,templist,df)

# 2. combine all 3 data frames
specdata <- rbind(specdata.e, specdata.r, specdata.n)
specdatalong <- specdata %>% gather(FLEP,value,1:5) #ready for ggplot

# 3. plot

# 
ggplot(data=specdatalong[specdatalong$variable.type == "eggs" & specdatalong$codNames=="GM",],
       aes(x=freq,y=value,color=FLEP)) + geom_line() + 
       scale_color_brewer(palette = "Reds") + theme_classic()



ggplot(forplot,aes(x=age,y=value,color=variable)) +
  ylab("egg production") +
  xlab("age") +
  geom_line() +
  scale_color_brewer(palette = "Reds") +
  ggtitle(paste(names(LSBlist)[i])) +
  theme_classic()










# ===================================================================
# 5) calculate area under curve for all pops (data is not normalized, not log):
AUC = rep(NA,length(eigentable$codNames))
for (i in 1:(length(eigentable$codNames))) {
  x = sp$freq
  y = sp[,i+1]
  id = order(x)
  AUC[i] = sum(diff(x[id]) * rollmean(y[id],2))
}
AUC = cbind.data.frame(eigentable$codNames,AUC)
colnames(AUC) <- c("codNames","varAUC")
rm(id,x,y) #clean up
# ---
# calculate variance in time steps, compare to AUC
varpops = rep(NA, length(eigentable$codNames))
for (i in 1:length(outputL)) {
  varpops[i] = var(outputL[[i]]$recruits[rm_first_timesteps:(timesteps-2)])}
dd <- cbind.data.frame(AUC,varpops)
colnames(dd) <- c("codNames","varAUC","varTS")
dd$varAUC - dd$varTS
#rm(dd)
# ---
# calculate CV of timeseries
cvpops = rep(NA, length(eigentable$codNames))
for (i in 1:length(outputL)) {
  cvpops[i] = (sd(outputL[[i]]$recruits[rm_first_timesteps:(timesteps-2)])/
                 mean(outputL[[i]]$recruits[rm_first_timesteps:(timesteps-2)]))*100}
cvpops = cbind.data.frame(eigentable$codNames,cvpops)
colnames(cvpops) <- c("codNames", "cvTS")


# ===================================================================
# 6) normalize the variance (spec)
maxs_rec <- apply(sp_rec[2:18],2,max) #find max variance for each pop
spp_rec = cbind( seq(from=0.001111111,to=0.5,by=0.001111111),
             t(t(sp_rec[2:18])/maxs_rec)) #divide var at all freq by max
colnames(spp_rec) = c("freq",eigentable$codNames)
spp_rec = as.data.frame(spp_rec) #use this dataframe to plot 

maxs_egg <- apply(sp_egg[2:18],2,max) #find max variance for each pop
spp_egg = cbind( seq(from=0.001111111,to=0.5,by=0.001111111),
                 t(t(sp_egg[2:18])/maxs_egg)) #divide var at all freq by max
colnames(spp_egg) = c("freq",eigentable$codNames)
spp_egg = as.data.frame(spp_egg) #use this dataframe to plot 

# take log() of spec before normalizing, then normalize
# ---
# recruits
sp_rec_log <- log(sp_rec[,2:18])
sp_rec_log$freq <- freq
maxs_rec_log <- apply(sp_rec_log[2:18],2,max) #find max variance for each pop
spp_rec_log = cbind( seq(from=0.001111111,to=0.5,by=0.001111111),
                 t(t(sp_rec_log[2:18])/maxs_rec_log)) #divide var at all freq by max
colnames(spp_rec_log) = c("freq",eigentable$codNames)
spp_rec_log = as.data.frame(spp_rec_log) #use this dataframe to plot 
# ---
# eggs
sp_egg_log <- log(sp_egg[,2:18])
sp_egg_log$freq <- freq
maxs_egg_log <- apply(sp_egg_log[2:18],2,max) #find max variance for each pop
spp_egg_log = cbind( seq(from=0.001111111,to=0.5,by=0.001111111),
                     t(t(sp_egg_log[2:18])/maxs_egg_log)) #divide var at all freq by max
colnames(spp_egg_log) = c("freq",eigentable$codNames)
spp_egg_log = as.data.frame(spp_egg_log) #use this dataframe to plot 


# ===================================================================
# 6) multiply the freq values for each population by the mode (shifts the frequency  
# plot so that cohort bump lines up at freq=1, in theory this should work)
mode_table <- as.data.frame(subset(eigentable, select=c("codNames","mode","temp")))
# eggs
spp <- spp_egg
long <- melt(spp,id="freq")
colnames(long) <- c("freq","codNames","spec")
splong_egg <- long %>% 
  left_join(mode_table) %>%
  mutate(freq.adjust = freq*mode) 
# log(eggs)
spp <- sp_egg_log
long <- melt(spp,id="freq")
colnames(long) <- c("freq","codNames","spec")
splong_egg_log <- long %>% 
  left_join(mode_table) %>%
  mutate(freq.adjust = freq*mode)
# recruits
spp <- spp_rec
long <- melt(spp,id="freq")
colnames(long) <- c("freq","codNames","spec")
splong_rec <- long %>% 
  left_join(mode_table) %>%
  mutate(freq.adjust = freq*mode)
# log(recruits)
spp <- sp_rec_log
long <- melt(spp,id="freq")
colnames(long) <- c("freq","codNames","spec")
splong_rec_log <- long %>% 
  left_join(mode_table) %>%
  mutate(freq.adjust = freq*mode)
rm(long,mode_table,spp) # clean up


# ===================================================================
# 7) plot all frequencies on one plot (each line color coded w/mode): 
# get ready to plot (default is 'spp' for normalized, choose 'sp' for original)

# ---
pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/freq_content_oneplot_eggs_centered.pdf', width=8, height=9)

# plot: eggs, color=mode, plots for Symposium presentation
splong <- splong_egg
ggplot(data=splong, aes(x=freq.adjust,y=spec, group=codNames)) + 
  geom_line(aes(color=mode)) +
  scale_x_continuous(limits=c(0,1.5) ) + 
                     #breaks = round(seq(min(splong$freq.adjust), 1.5, by = 0.5),1)) +
  #scale_y_log10(limits=c(2e-03,1)) + #adjusting scale
  coord_cartesian(xlim = c(0, 2)) + #adjusting scale
  geom_text_repel(data=as.data.frame(splong %>%  
                                       group_by(codNames) %>% 
                                       arrange(abs(freq.adjust-1.5)) %>% 
                                       slice(1) %>%
                                       mutate(freq.adjust=round(freq.adjust,1))),
                  aes(label = codNames,color=mode),
                  nudge_x = 1,segment.color = "grey",
                  size = 6,
                  na.rm = TRUE) +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("spectral analysis (eggs)") +
  geom_vline(xintercept=1,linetype="dashed")
dev.off()

# 2J3KL & Celtic pops only:
pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/symposium_figs/freq_content_oneplot_eggs_centered_extremes.pdf', width=8, height=9)
ggplot(data=splong[splong$codNames %in% c("2J3KL","Celtic"),], aes(x=freq.adjust,y=spec, group=codNames)) + 
  geom_line(aes(color=mode)) +
  scale_x_continuous(limits=c(0,1.5), 
                     breaks = round(seq(min(splong$freq.adjust), 1.5, by = 0.5),1)) +
  scale_y_log10(limits=c(2e-03,1)) + #adjusting scale
  coord_cartesian(xlim = c(0, 2)) + #adjusting scale
  geom_text_repel(data=as.data.frame(splong[splong$codNames %in% c("2J3KL","Celtic"),] %>%  
                                       group_by(codNames) %>% 
                                       arrange(abs(freq.adjust-1.5)) %>% 
                                       slice(1) %>%
                                       mutate(freq.adjust=round(freq.adjust,1))),
                  aes(label = codNames, color=mode),
                  nudge_x = 1, segment.color = "grey",
                  size = 6,
                  na.rm = TRUE) +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("spectral analysis in egg production") 
#geom_vline(xintercept=1,linetype="dashed")
dev.off()

# ===================================================================
# plot: recruits, normalized, color=mode, freq adjusted
pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/freq_content_oneplot_recruits_eggs.pdf', width=8, height=9)
splong <- splong_rec
rec.p <- ggplot(data=splong, aes(x=freq.adjust,y=spec, group=codNames)) + 
  geom_line(aes(color=mode)) +
  scale_x_continuous(limits=c(0,1.5), 
                     breaks = round(seq(min(splong$freq.adjust), 1.5, by = 0.5),1)) +
  scale_y_log10(limits=c(2e-03,1)) + #adjusting scale
  coord_cartesian(xlim = c(0, 2)) + #adjusting scale
  geom_text_repel(data=as.data.frame(splong %>%  
                                       group_by(codNames) %>% 
                                       arrange(abs(freq.adjust-1.5)) %>% 
                                       slice(1) %>%
                                       mutate(freq.adjust=round(freq.adjust,1))),
                  aes(label = codNames,color=mode),
                  nudge_x = 1,segment.color = "grey",
                  size = 4,
                  na.rm = TRUE) +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("spectral analysis (recruits)") +
  geom_vline(xintercept=1,linetype="dashed")
# --------- color=temp plots not working ------------ #
# plot: plot: recruits, normalized, color=temp, freq adjusted
ggplot(data=splong, aes(x=freq.adjust,y=spec, group=codNames)) + 
  geom_line(aes(color=temp)) +
  #scale_color_gradient2(midpoint=mean(splong$temp,low="blue",high="red")) +
  scale_color_gradientn(colors=rev(rainbow(n=16,start=0,end=0.7))) +
  scale_x_continuous(limits=c(0,1.5), 
                     breaks = round(seq(min(splong$freq.adjust), 1.5, by = 0.5),1)) +
  scale_y_log10(limits=c(2e-03,1)) + #adjusting scale
  coord_cartesian(xlim = c(0, 2)) + #adjusting scale
  geom_text_repel(data=as.data.frame(splong %>%  
                                       group_by(codNames) %>% 
                                       arrange(abs(freq.adjust-1.5)) %>% 
                                       slice(1) %>%
                                       mutate(freq.adjust=round(freq.adjust,1))),
                  aes(label = codNames,color=temp),
                  nudge_x = 1,segment.color = "grey",
                  size = 6,
                  na.rm = TRUE) +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("spectral analysis (recruits)") +
  geom_vline(xintercept=1,linetype="dashed")
# -------------------- #

# plot: eggs, normalized, color=mode, freq adjusted
splong <- splong_egg
egg.p <- ggplot(data=splong, aes(x=freq.adjust,y=spec, group=codNames)) + 
  geom_line(aes(color=mode)) +
  scale_x_continuous(limits=c(0,1.5), 
                     breaks = round(seq(min(splong$freq.adjust), 1.5, by = 0.5),1)) +
  scale_y_log10(limits=c(2e-03,1)) + #adjusting scale
  coord_cartesian(xlim = c(0, 2)) + #adjusting scale
  geom_text_repel(data=as.data.frame(splong %>%  
                                       group_by(codNames) %>% 
                                       arrange(abs(freq.adjust-1.5)) %>% 
                                       slice(1) %>%
                                       mutate(freq.adjust=round(freq.adjust,1))),
                  aes(label = codNames,color=mode),
                  nudge_x = 1,segment.color = "grey",
                  size = 4,
                  na.rm = TRUE) +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("spectral analysis (eggs)") +
  geom_vline(xintercept=1,linetype="dashed")

# plot: recruits, normalized, color=mode, freq adjusted, log(spec)
splong <- splong_rec_log
reclog.p <- ggplot(data=splong, aes(x=freq.adjust,y=spec, group=codNames)) + 
  geom_line(aes(color=mode)) +
  scale_x_continuous(limits=c(0,1.5), 
                     breaks = round(seq(min(splong$freq.adjust), 1.5, by = 0.5),1)) +
 # scale_y_log10(limits=c(2e-03,1)) + #adjusting scale
  coord_cartesian(xlim = c(0, 2)) + #adjusting scale
  geom_text_repel(data=as.data.frame(splong %>%  
                                       group_by(codNames) %>% 
                                       arrange(abs(freq.adjust-1.5)) %>% 
                                       slice(1) %>%
                                       mutate(freq.adjust=round(freq.adjust,1))),
                  aes(label = codNames,color=mode),
                  nudge_x = 1,segment.color = "grey",
                  size = 4,
                  na.rm = TRUE) +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("spectral analysis log(recruits)") +
  geom_vline(xintercept=1,linetype="dashed")

splong <- splong_egg_log
egglog.p <- ggplot(data=splong, aes(x=freq.adjust,y=spec, group=codNames)) + 
  geom_line(aes(color=mode)) +
  scale_x_continuous(limits=c(0,1.5), 
                     breaks = round(seq(min(splong$freq.adjust), 1.5, by = 0.5),1)) +
  scale_y_log10(limits=c(2e-03,1)) + #adjusting scale
  coord_cartesian(xlim = c(0, 2)) + #adjusting scale
  geom_text_repel(data=as.data.frame(splong %>%  
                                       group_by(codNames) %>% 
                                       arrange(abs(freq.adjust-1.5)) %>% 
                                       slice(1) %>%
                                       mutate(freq.adjust=round(freq.adjust,1))),
                  aes(label = codNames,color=mode),
                  nudge_x = 1,segment.color = "grey",
                  size = 4,
                  na.rm = TRUE) +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("spectral analysis log(eggs)") +
  geom_vline(xintercept=1,linetype="dashed")
# export to PDFs
# 2 plots: spectral analysis of eggs and recruits
grid.arrange(rec.p,egg.p,nrow=2)
pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/freq_content_oneplot_recruits_eggs.pdf', width=6, height=10)
grid.arrange(rec.p,egg.p,nrow=2)
dev.off()

pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/freq_content_oneplot_recruits_eggs_log.pdf', width=6, height=10)
grid.arrange(reclog.p,egglog.p,nrow=2)
dev.off()









# ===================================================================
# these are extra plots, just in case for future
# ---
# plot: total variance as a function of max age, mode, and CV in spawing age distribution
# gather relevant data:
ddplot <- cbind.data.frame(eigentable,varpops)
# mode
mode <- ggplot(ddplot, aes(x=mode,y=varpops)) +
  theme_set(theme_classic(base_size = 8)) +
  #theme(axis.title.y=element_text(size=rel(1))) +
  #theme(axis.title.x=element_text(size=rel(1))) +
  geom_point(aes(color=temp), size=3) +
  scale_color_gradientn(colors=rev(rainbow(n=10,start=0,end=0.7))) +
  geom_text(aes(label=codNames),size = 2,hjust=-0.1,vjust=-0.2,check_overlap = F) +
  ylab("total variance") +
  xlab("mode of the spawning age distribution") +
  theme(plot.title = element_text(hjust = 0.5)) 
# max age
max <- ggplot(ddplot, aes(x=max_ages,y=varpops)) +
  theme_set(theme_classic(base_size = 8)) +
  #theme(axis.title.y=element_text(size=rel(1))) +
  #theme(axis.title.x=element_text(size=rel(1))) +
  geom_point(aes(color=temp), size=3) +
  scale_color_gradientn(colors=rev(rainbow(n=10,start=0,end=0.7))) +
  geom_text(aes(label=codNames),size = 2,hjust=-0.1,vjust=-0.2,check_overlap = F) +
  ylab("total variance") +
  xlab("max age") +
  theme(plot.title = element_text(hjust = 0.5)) 
# stdev about the mode
stdev <- ggplot(ddplot, aes(x=sd_mode,y=varpops)) +
  theme_set(theme_classic(base_size = 8)) +
  #theme(axis.title.y=element_text(size=rel(1))) +
  #theme(axis.title.x=element_text(size=rel(1))) +
  geom_point(aes(color=temp), size=3) +
  scale_color_gradientn(colors=rev(rainbow(n=10,start=0,end=0.7))) +
  geom_text(aes(label=codNames),size = 2,hjust=-0.1,vjust=-0.2,check_overlap = F) +
  ylab("total variance") +
  xlab("stdev about the mode") +
  theme(plot.title = element_text(hjust = 0.5)) 
# temp
temp <- ggplot(ddplot, aes(x=temp,y=varpops)) +
  theme_set(theme_classic(base_size = 8)) +
  #theme(axis.title.y=element_text(size=rel(1))) +
  #theme(axis.title.x=element_text(size=rel(1))) +
  geom_point(aes(color=temp), size=3) +
  scale_color_gradientn(colors=rev(rainbow(n=10,start=0,end=0.7))) +
  geom_text(aes(label=codNames),size = 2,hjust=-0.1,vjust=-0.2,check_overlap = F) +
  ylab("total variance") +
  xlab("temp") +
  theme(plot.title = element_text(hjust = 0.5)) 
pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/totalvariance_multiplot.pdf', width=11, height=8)
grid.arrange(mode,max,stdev,temp,ncol=2)
dev.off()
# ---
dampratio <- ggplot(ddplot, aes(x=dampratio,y=varpops)) +
  theme_set(theme_classic(base_size = 8)) +
  #theme(axis.title.y=element_text(size=rel(1))) +
  #theme(axis.title.x=element_text(size=rel(1))) +
  geom_point(aes(color=temp), size=3) +
  scale_color_gradientn(colors=rev(rainbow(n=10,start=0,end=0.7))) +
  geom_text(aes(label=codNames),size = 2,hjust=-0.1,vjust=-0.2,check_overlap = F) +
  ylab("total variance") +
  xlab(expression(paste(lambda[1],"/","|",lambda[2],"|"," (damping ratio, ",lambda[1]," adjusted to 1)"))) +
  theme(plot.title = element_text(hjust = 0.5)) 
pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/totalvariance_dampratio.pdf', width=6, height=7)
dampratio
dev.off()


# ===================================================================
# I think the code below doesn't work, but saving just in case
# ---
# plot timeseries, can be eggs, recruits, or N
#pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/timeseries1.pdf', width=7, height=10)
#par(mfrow=c(5,3))
# ---
for (i in 1:length(eigentable$codNames)) { # step through each cod population
  A = read.table(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/matrix_maxages/'
                              ,eigentable$codNames[i], '.txt', sep=''))
  A = as.matrix(A)
  # run simulation 
  timesteps = 1000
  rm_first_timesteps = 900
  alpha = 100
  beta = 1000
  initial_eggs = 1
  sig_r = 0.1
  output = sim_model(A=A, timesteps=timesteps, 
                     alpha=alpha, beta=beta, 
                     sig_r=sig_r, initial_eggs=initial_eggs)  
  # plot Nt over time
  #plot(x=1:length(output[[1]]),
  #y=output[[1]],
  #main=eigentable$codNames[i],
  #type="l", ylab="eggs")
  tspt[,i] = output[[2]][rm_first_timesteps:(timesteps-2)]
  plot(x=1:(length(output[[2]][rm_first_timesteps:(timesteps-2)])),
       y=output[[2]][rm_first_timesteps:(timesteps-2)],
       main=eigentable$codNames[i],
       type="l", ylab="eggs")
}
dev.off()
par(mfrow=c(1,1))
tspt = cbind(c(1:899),tspt)
colnames(tspt) <- c("time",eigentable$codNames)
tslong = melt(as.data.frame(tspt),id="time") #long form for ggplot
mycolors = c(brewer.pal(name="Set1", n = 8), brewer.pal(name="Paired", n = 8)) #16 pops, but pallettes only have 8

ggplot(data=tslong, aes(x=time,y=value,color=variable,linetype=variable)) + 
  geom_line(lwd=1) + theme_classic() +
  scale_linetype_manual(values=c(rep("solid",10),rep("dashed",6))) +
  scale_color_manual(values=mycolors) +
  ggtitle("N")


# plot eggs, recruits, Nsize
par(mfrow = c(2,2))
plot(outputL$Northsea[[4]],type="l",main="eggs")
plot(x=1:timesteps,y=output[[3]],type="l",main="recruits (before variability)")
plot(x=1:timesteps,y=output[[4]],type="l",main="N")
plot(x=1:timesteps,y=output[[1]][1,],type="l",main="recruits (after variability)")
par(mfrow = c(1,1))

output$Nt
output$recruits[1:20]
plot(output$eggs[100:999],type="l")
matplot(t(output$Nt[,1:20]), type="l")



