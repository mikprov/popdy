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

# load maxages
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/6_find_and_export_oldest_fish_in_each_pop.r")
max_ages_table
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
# (2) Simulate pops. Loop over Aarray list to simulate using different Leslie matrices 
# *************************************** #
# set params for simulation:
timesteps = 1000 #need this now to create
rm_first_timesteps = 200
alpha = 1.2
beta = 700
initial_eggs = beta
sig_r = 0.3
span.multiplier = 1 # what is this again?
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
# (4) Now that timeseries data is formated, let's plot!
# *************************************** #

#plot recruitment - one plot per pop, similar to egg plots
p <- list()
codNames_ordered_by_peak <- eigentable %>% arrange(mode_age) %>% pull(codNames)
sd_ordered <- rep(NA,length=length(codNames_ordered_by_peak))
mean_ordered <- rep(NA,length=length(codNames_ordered_by_peak))
CV_ordered <- rep(NA,length=length(codNames_ordered_by_peak))

for (i in 1:length(codNames)){
  dd <- ts.data[ts.data$variable == "recruits" & 
                  ts.data$codNames == codNames_ordered_by_peak[i] &
                  ts.data$year %in% seq(from=rm_first_timesteps,to=(timesteps-2),by=1),]
  
  p[[i]] <- ggplot(dd, #aes(x=year,y=value,color=Fval)) +
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
do.call(grid.arrange,c(p,ncol=3))

#par(mfrow=c(1,1))
#plot(x=eigentable %>% arrange(mode_age) %>% pull(mode_age),y=mean_ordered,
#     xlab="peak spawning age", ylab="recruits timseries mean")
#plot(x=eigentable %>% arrange(mode_age) %>% pull(mode_age),y=sd_ordered,
#     xlab="peak spawning age", ylab="recruits timseries Stdev")
#plot(x=eigentable %>% arrange(mode_age) %>% pull(mode_age),y=CV_ordered,
#     xlab="peak spawning age", ylab="recruits timseries CV")
#par(mfrow=c(1,1))
#dev.off()
#rm(p,i,dd)


# plot eggs - one plot per pop, on each plot 5 lines for different F levels
p <- list()
# reorder codNames by peak spawning age (increasing)
codNames_ordered_by_peak <- eigentable %>% arrange(mode_age) %>% pull(codNames)
sd_ordered <- rep(NA,length=length(codNames_ordered_by_peak))
mean_ordered <- rep(NA,length=length(codNames_ordered_by_peak))
CV_ordered <- rep(NA,length=length(codNames_ordered_by_peak))

for (i in 1:length(codNames_ordered_by_peak)){
  
  dd <- ts.data[ts.data$variable == "eggs" & 
                  ts.data$codNames == codNames_ordered_by_peak[i] &
                  ts.data$year %in% seq(from=rm_first_timesteps,to=(timesteps-2),by=1),]
  
  p[[i]] <- ggplot(dd, #aes(x=year,y=value,color=Fval)) +
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
do.call(grid.arrange,c(p,ncol=3))

#par(mfrow=c(2,2))
#plot(x=eigentable %>% arrange(mode_age) %>% pull(mode_age),y=mean_ordered,
#           xlab="peak spawning age", ylab="eggs timseries mean")
#plot(x=eigentable %>% arrange(mode_age) %>% pull(mode_age),y=sd_ordered,
#           xlab="peak spawning age", ylab="eggs timseries Stdev")
#plot(x=eigentable %>% arrange(mode_age) %>% pull(mode_age),y=CV_ordered,
#           xlab="peak spawning age", ylab="eggs timseries CV")
#par(mfrow=c(1,1))
#dev.off()
#rm(p,i,dd)



#plot pop size - one plot per pop, similar to egg & recruitment plots
p <- list()
for (i in 1:length(codNames)){
  p[[i]] <- ggplot(ts.data[ts.data$variable == "Nsize" & 
                             ts.data$codNames == codNames[i],], 
                   aes(x=year,y=value,color=Fval)) + 
    xlab("year") + ylab("Nsize") +
    geom_line() + scale_color_brewer(palette = "Reds") +
    ggtitle(paste(codNames[i])) +
    theme_classic()
}
#pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/Nsize_for_diff_FLEPs_alpha50.pdf', width=7, height=10) #note: file name specifies the alpha used in simluation model
do.call(grid.arrange,c(p,ncol=2))
#dev.off()
rm(p,i)


# *************************************** #
# (5) Calculate frequency content from timeseries
# *************************************** #
# Plan:
# 1. Walk through each cod pop, do spectral analysis at F levels
# 2. Store spec values for eggs, recruits, and Nsize
#del <- 0.001111111
#freq = seq(from=0.0025, to=0.5, by=0.0025) #frequencies 
#length(freq)
#dataNAs <- rep(NA,length=length(freq)*length(Fvalues)*length(codNames)) #
sp.eggsL <- as.list(rep(NA,length=length(codNames))) #object for spec analysis  
sp.recruitL <- as.list(rep(NA,length=length(codNames))) 
#sp.nsize <- as.data.frame(matrix(NA,nrow=0,ncol=length(codNames)))#2D object for spec analysis 

for (i in 1:length(codNames)){
  
  ts <- ts.data[ts.data$codNames == codNames[i],] #subset data for pop i 
  
  # setting 'span' - a vector of odd integers to specify the smoothers
  tmp <- ceiling(sqrt(length(1:(timesteps-rm_first_timesteps-1)))) #sq root of timeseries lgth, rounded
  if (tmp %% 2 == 0) {m <- tmp+1} else {m <- tmp} #make it odd, if the square root is even
  m = m * span.multiplier
  
  # --- spectral analysis on EGGS --- #
  xx = ts[ts$variable == "eggs",]$value[rm_first_timesteps:(timesteps-2)]-mean(ts[ts$variable == "eggs",]$value[rm_first_timesteps:(timesteps-2)])
  sp = spec.pgram(x=xx,spans=c(m,m),plot = FALSE)
  #save spec output for plotting pops together, Helen Wearing says to multiply by 2
  sp.eggsL[[i]] <- sp$spec*2
  rm(xx,sp)
  
  # --- spectral analysis on RECRUIT --- #
  yy = ts[ts$variable == "recruits",]$value[rm_first_timesteps:(timesteps-2)]-mean(ts[ts$variable == "recruits",]$value[rm_first_timesteps:(timesteps-2)])
  sp = spec.pgram(yy,spans=c(m,m),plot = FALSE)
  sp.recruitL[[i]] = 2*sp$spec # save matrix of spec values for different FLEP, index by pop i
  rm(yy)
  
  # --- spectral analysis on NSIZE --- #
  #n = matrix(NA, nrow=length(freq), ncol=length(Fvalues)) #store nsize sp vals at FLEP levels here
  #zz = ts[ts$variable == "eggs" & ts$Fval == Fvalues[f],]$value[rm_first_timesteps:(timesteps-2)]-mean(ts[ts$variable == "eggs" & ts$Fval == Fvalues[f],]$value[rm_first_timesteps:(timesteps-2)])
  
  #for (f in 1:length(Fvalues)){ #step through FLEP levels
    #sp = spec.pgram(zz,spans=c(m,m),plot = FALSE)
    #nn = 2*sp$spec
    #names(nn) <- Fvalues[f]
    #n[,f] <- nn
  #} #save spec output for plotting pops together, Helen Wearing says to multiply by 2
  #rm(zz)
  
  
  
  #sp.nsize[,,i] = n # save matrix of spec values for different FLEP, index by pop i
  
}
freq <- sp$freq
sp.eggs <- as.data.frame(do.call(cbind,sp.eggsL))
sp.recruit <- as.data.frame(do.call(cbind,sp.recruitL))
names(sp.eggs) <- codNames
names(sp.recruit) <- codNames
sp.eggs$freq <- freq
sp.recruit$freq <- freq

rm(i,f,sp,e,r,n,tmp,m) #clean up


# *************************************** #
# (9) Plot spectral analysis 
# *************************************** #

# Plan
# 1. re-arrange 3d dataframes for plotting (sp.eggs, sp.recruit, sp.nsize)
# 2. combine all 3 variable type data frames to make one df 
# 3. plot

egglong <- sp.eggs %>% gather("codName","value",1:length(codNames))
egglong$variable.type <- rep("eggs",length=length(egglong$freq))
head(egglong)
recruitlong <- sp.recruit %>% gather("codName","value",1:length(codNames))
recruitlong$variable.type <- rep("recruit",length=length(recruitlong$freq))
head(recruitlong)
specdatalong <- rbind(egglong,recruitlong)
head(specdatalong)


# 1. re-arrange egg dataframe
templist <- as.list(rep(NA,length=length(codNames)))
names(templist) <- codNames
for (i in 1:length(codNames)){
  df <- as.data.frame(sp.eggs[,,i])
  colnames(df) <- c("value")
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
  colnames(df) <- c("value")
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
  colnames(df) <- c("value")
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
specdatalong <- specdata
#specdatalong <- specdata %>% gather(FLEP,value,1) #ready for ggplot

# here I can change the order I want plots to appear:
eigentable
codNames_ordered_by_peak <- eigentable %>% arrange(mode_age) %>% pull(codNames)


# 3. plot frequency content of eggs
p <- list()
for (i in 1:length(codNames_ordered_by_peak)){
  p[[i]] <- ggplot(data=specdatalong[specdatalong$variable.type == "eggs" 
                                   & specdatalong$codNames==codNames_ordered_by_peak[i],],
                 #aes(x=freq,y=log10(value),color=FLEP)) + 
                 aes(x=freq,y=log10(value))) + 
    geom_line() + #ylim(9,13) +
    geom_vline(xintercept = (1/eigentable[eigentable$codNames == codNames_ordered_by_peak[i],]$mode_age),
             linetype="dotted") +
    ylab("log(eggs)") +
    ggtitle(paste(codNames_ordered_by_peak[i]," mode=",
                eigentable[eigentable$codNames == codNames_ordered_by_peak[i],]$mode_age," max=",
                eigentable[eigentable$codNames == codNames_ordered_by_peak[i],]$max_ages)) + 
    theme(plot.title = element_text(size = 10)) + theme_classic()
#+ scale_color_brewer(palette = "Reds")

}
pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/spec_eggs_alpha50_log10_nofishing.pdf', width=7, height=14) #note: file name specifies the alpha used in simluation model
do.call(grid.arrange,c(p,ncol=3))
dev.off()
rm(p,i)

# -- recruits plotting -- #
p <- list()
for (i in 1:length(codNames_ordered_by_peak)){
  p[[i]] <- ggplot(data=specdatalong[specdatalong$variable.type == "recruit" 
                                     & specdatalong$codNames==codNames_ordered_by_peak[i],],
                  # aes(x=freq,y=log10(value),color=FLEP)) + 
                   aes(x=freq,y=log10(value))) + 
    geom_line()  + theme_classic() + ylim(1,5) +
    geom_vline(xintercept = (1/eigentable[eigentable$codNames == codNames_ordered_by_peak[i],]$mode_age),
               linetype="dotted") +
    ggtitle(paste(codNames_ordered_by_peak[i]," mode=",
                  eigentable[eigentable$codNames == codNames_ordered_by_peak[i],]$mode_age)) + 
    ylab("log(recruits)")
    #+ scale_color_brewer(palette = "Reds")
  
}
pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/spec_recruitsbeforenoise_alpha50_log10.pdf', width=7, height=14) #note: file name specifies the alpha used in simluation model
do.call(grid.arrange,c(p,ncol=3))
dev.off()
rm(p,i)

# -- Nsize plotting -- #
p <- list()
for (i in 1:length(codNames)){
  p[[i]] <- ggplot(data=specdatalong[specdatalong$variable.type == "nsize" 
                                     & specdatalong$codNames==codNames[i],],
                   aes(x=freq,y=value,color=FLEP)) + 
    geom_line() + scale_color_brewer(palette = "Reds") + theme_classic() + 
    ggtitle(paste(codNames[i]))
  
}
pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/spec_Nsize_alpha50.pdf', width=7, height=14) #note: file name specifies the alpha used in simluation model
do.call(grid.arrange,c(p,ncol=2))
dev.off()
rm(p,i)


# *************************************** #
# (10) Area under curve
# *************************************** #
head(specdatalong)
AUC_greater <- rep(NA,length=length(codNames_ordered_by_peak))
AUC_less <- rep(NA,length=length(codNames_ordered_by_peak))
AUC_total <- rep(NA,length=length(codNames_ordered_by_peak))
AUCperlow <- rep(NA,length=length(codNames_ordered_by_peak))
AUCperhigh <- rep(NA,length=length(codNames_ordered_by_peak))

AUCthreshold <- 0.3

for (i in 1:length(codNames_ordered_by_peak)){
  AUC_greater[i] <- sum(specdatalong[specdatalong$variable.type == "eggs" 
                             & specdatalong$codNames==codNames_ordered_by_peak[i]
                             & specdatalong$freq > AUCthreshold,]$value)
  AUC_less[i] <- sum(specdatalong[specdatalong$variable.type == "eggs" 
                                     & specdatalong$codNames==codNames_ordered_by_peak[i]
                                     & specdatalong$freq < AUCthreshold,]$value)
  AUC_total[i] <- sum(specdatalong[specdatalong$variable.type == "eggs" 
                                & specdatalong$codNames==codNames_ordered_by_peak[i],]$value)
  AUCperlow[i] <- AUC_less[i]/AUC_total[i]
  AUCperhigh[i] <- AUC_greater[i]/AUC_total[i]
}
par(mfrow=c(2,2))
plot(x=eigentable$mode_age, y=AUCperhigh,main="threshold=0.05",xlab="peak spawning age")
plot(x=eigentable$mode_age, y=AUCperhigh,main="threshold=0.1",xlab="peak spawning age")
plot(x=eigentable$mode_age, y=AUCperhigh,main="threshold=0.2",xlab="peak spawning age")
plot(x=eigentable$mode_age, y=AUCperhigh,main="threshold=0.3",xlab="peak spawning age")


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



