# Run the simulation model
# by: mikaela provost
# last edited: Jan 25, 2019
# ===================================================================

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
# read in eigentable - I'm using some information from the table
eigentable = read.csv("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLSB/eigentable.csv",
                      header=TRUE,stringsAsFactors = FALSE)
eigentable = as.data.frame(eigentable)

# ---
# create empty objects for simulation
freq = seq(from=0.001111111, to=0.5, by=0.001111111) #all frequencies, for plotting later
spall = matrix(NA, length(freq),length(eigentable$codNames)) #set up empty matrix, spec for all plots
outputL = vector("list", length(eigentable$codNames)) #create empty list to store sim output
names(outputL) = eigentable$codNames

# ---
# set params for simulation:
timesteps = 1000
rm_first_timesteps = 100
beta = 10000 #note: alpha is different for each pop
initial_eggs = 1000
sig_r = 0.1
span.multiplier = 1 # what is this again?

# ---
# store 3d arrays of Leslie matricies for each Fvalues for each pop
Aarray <- as.list(rep(NA,length=length(codNames)))
names(Aarray) <- codNames
# store FLEP and LEP information
FLEPinfolist <- as.list(rep(NA,length=length(codNames)))
# store alpha values at 1/(LEP*.35) for each pop i
alphas <- rep(NA, length=length(codNames))

# --- #
# choose the F values associated with FLEP values
Fvalues = seq(0,200,by=0.01) #I may need to change the max F value (some pops can withstand high F)
FLEP = matrix(NA,nrow=length(Fvalues),ncol=length(codNames)) #store FLEP values here
colnames(FLEP) = codNames

# ---
# run simulation - EGGS --> that means output[[2]]
for (i in 1:length(eigentable$codNames)) { # step through each cod population
  
  # load parms for cod pop i
  source(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_pops/',codNames[i], '.r', sep=''))
  # this should load parms: L_inf, K (vonB), TEMP, maxage
  LSBmatrix = calculate_LSB_at_age_by_F(maxage=maxage, 
                                             L_inf=L_inf, 
                                             K=K, 
                                             TEMP=TEMP, 
                                             F.halfmax=Fvalues,
                                             B0=B0,B1=B1)
  # here I calculate LEP by summing LSB at age across all ages
  LEP <- colSums(LSBmatrix)
  # here I calculate FLEP: LSB at all F levels / LSB at F=0
  FLEP[,i] <- colSums(LSBmatrix) / sum(LSBmatrix[,1]) 
  FLEP.F <- as.data.frame(cbind(Fvalues,FLEP[,i])) #first col=F values, 2nd col=pop i, rows=FLEP values
  colnames(FLEP.F) <- c("Fvalues","FLEP") 
  
  # --- #
  # find the nearest F values that correspond to the desired FLEP value
  target_Fs <- rbind(
    
    as.data.frame(melt(FLEP.F,id="Fvalues") %>% 
                    group_by(variable) %>%
                    arrange(abs(value-1)) %>% #desired FLEP=1
                    slice(1) %>%
                    #mutate(value=round(value,2)) %>%
                    mutate(FLEPlevel=rep(1))),
    
    as.data.frame(melt(FLEP.F,id="Fvalues") %>% 
                    group_by(variable) %>%
                    arrange(abs(value-0.8)) %>% #desired FLEP=0.8
                    slice(1) %>%
                    #mutate(value=round(value,2)) %>%
                    mutate(FLEPlevel=rep(0.8))),
    
    as.data.frame(melt(FLEP.F,id="Fvalues") %>% 
                    group_by(variable) %>%
                    arrange(abs(value-0.5)) %>% #desired FLEP=0.5
                    slice(1) %>%
                    #mutate(value=round(value,2)) %>%
                    mutate(FLEPlevel=rep(0.5))),
    
    as.data.frame(melt(FLEP.F,id="Fvalues") %>% 
                    group_by(variable) %>%
                    arrange(abs(value-0.35)) %>% #desired FLEP=0.35
                    slice(1) %>%
                    #mutate(value=round(value,2)) %>%
                    mutate(FLEPlevel=rep(0.35))),
    
    as.data.frame(melt(FLEP.F,id="Fvalues") %>% 
                    group_by(variable) %>%
                    arrange(abs(value-0.2)) %>% #desired FLEP=0.2
                    slice(1) %>%
                    #mutate(value=round(value,2)) %>%
                    mutate(FLEPlevel=rep(0.2))) 
    
  )
  # Check: did the range of F values deplete the population low enough to 
  # reach FLEP value of 0.2? AKA: does calcFLEP = targetFLEP? If not, 
  # then I need to increase the F range (...fish harder!)
  #target_Fs
  
  # --- #
  # here, I extract the LEP values associated with each FLEP level
  # I could probably do this faster, but this works for now. 
  FLEP.LEP <- as.data.frame(cbind(LEP,FLEP[,i]))
  colnames(FLEP.LEP) <- c("LEP","FLEP")
  # find the nearest LEP values that correspond to the desired FLEP value
  target_LEPs <- rbind(
    
    as.data.frame(melt(FLEP.LEP,id="LEP") %>% 
                    group_by(variable) %>%
                    arrange(abs(value-1)) %>% #desired FLEP=1
                    slice(1) %>%
                    #mutate(value=round(value,2)) %>%
                    mutate(FLEPlevel=rep(1))),
    
    as.data.frame(melt(FLEP.LEP,id="LEP") %>% 
                    group_by(variable) %>%
                    arrange(abs(value-0.8)) %>% #desired FLEP=0.8
                    slice(1) %>%
                    #mutate(value=round(value,2)) %>%
                    mutate(FLEPlevel=rep(0.8))),
    
    as.data.frame(melt(FLEP.LEP,id="LEP") %>% 
                    group_by(variable) %>%
                    arrange(abs(value-0.5)) %>% #desired FLEP=0.5
                    slice(1) %>%
                    #mutate(value=round(value,2)) %>%
                    mutate(FLEPlevel=rep(0.5))),
    
    as.data.frame(melt(FLEP.LEP,id="LEP") %>% 
                    group_by(variable) %>%
                    arrange(abs(value-0.35)) %>% #desired FLEP=0.35
                    slice(1) %>%
                    #mutate(value=round(value,2)) %>%
                    mutate(FLEPlevel=rep(0.35))),
    
    as.data.frame(melt(FLEP.LEP,id="LEP") %>% 
                    group_by(variable) %>%
                    arrange(abs(value-0.2)) %>% #desired FLEP=0.2
                    slice(1) %>%
                    #mutate(value=round(value,2)) %>%
                    mutate(FLEPlevel=rep(0.2)))
    
  )
  target_Fs$LEP <- target_LEPs$LEP
  FLEPinfo <- target_Fs %>% select("Fvalues","FLEPlevel","LEP")
  FLEPinfo$codNames <- rep(codNames[i],length=length(FLEPinfo$Fvalues))
  #rm(target_Fs,target_LEPs,FLEP,FLEP.F,FLEP.LEP)
  
  # --- #
  # generate Leslie matricies for different F values
  # note: parms for cod pop i should already be loaded.
  Alist = as.list(rep(NA,length(FLEPinfo$Fvalues))) #Leslie matrix storage for each F value for pop i
  eigenvals1 = rep(NA,length=length(FLEPinfo$Fvalues))
  eigenvals2 = rep(NA,length=length(FLEPinfo$Fvalues))
  eigenvals12 = rep(NA,length=length(FLEPinfo$Fvalues))
  
  for(f in 1:length(FLEPinfo$Fvalues)){ #for each F value (associated with a FLEP level)
  
    # calculate the Leslie matrix:
    Leslieout = assemble_Leslie(maxage=maxage, K=K, L_inf=L_inf, TEMP=TEMP,
                      F.halfmax=FLEPinfo$Fvalues[f], tknot=0, B0=B0, B1=B1)
    Alist[[f]]=Leslieout$A #store Leslie in a list just in case I need it later
    A = Leslieout$A #store Leslie as A for later use
    
    # using the Leslie, calculate first and second eigenvalues
    eigenvals1[f] <- extract_first_eigen_value(A) #calc first eigenvalue
    eigenvals2[f] <- extract_second_eigen_value(A) #calc second eigenvalue
    eigenvals12[f] <- eigenvals2[f] / eigenvals1[f] #calc damp ratio
  
  }
  # collapse list of Leslie matricies into 3d array for pop i
  Aarray[[i]] <- array(unlist(Alist),dim=c(length(A[,1]),length(A[1,]),length(FLEPinfo$Fvalues)))
  FLEPinfo$lambda1 <- eigenvals1
  FLEPinfo$lambda2 <- eigenvals2
  FLEPinfo$lambda12 <- eigenvals12
  
  FLEPinfolist[[i]] <- FLEPinfo
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
} 
rm(i,eigenvals1,eigenvals2,eigenvals12,Leslieout,A,FLEPinfo,
   target_Fs,target_LEPs,Alist,LEP,FLEP,FLEP.F,FLEP.LEP,
   L_inf,K,TEMP,B0,B1,MG,Mp,name,rm_first_timesteps,S50,f) #clean up

# At this point I have some important objects:
# 1. [alphas] a vector of alpha values (defined as 1/(LEP*0.35)) for each pop
# 2. [Aarray] a list of arrays: each element in list is a stack of Leslie matrcies at each F
# 3. [FLEPinfolist] a list of dataframes with FLEP info: F values, lambda 1 and 2, LEP


# ------------------------------------------------------------------------- #
# --- Loop over Aarray list to simulate using different Leslie matrices 
# ------------------------------------------------------------------------- #

# There are different ways to run the simulation:
# 1. Use different alpha values for each pop (alpha = 1/(35%LEP)) --> alphas
# 2. Use all the same alpha values for each pop (alpha = pick a value) --> constant alpha
alphas
constant_alpha0.5 <- rep(0.5,length=length(alphas))
constant_alpha1 <- rep(1, length=length(alphas))
constant_alpha2 <- rep(2, length=length(alphas))
constnat_alpha3 <- rep(3, length=length(alphas))

output.3d.list <- as.list(rep(NA,length=length(codNames))) #store timeseries here
names(output.3d.list) <- codNames

for (i in 1:length(Aarray)) { #step through each pop
  Leslie3d = Aarray[[i]] #select the 3d array of Leslie matricies
  # array dims: 5 is number of F values, 4 is number of ts, 998 is length of ts
  output.matrix <- array(NA,c(998,4,5)) 
  
  for (f in 1:5) { #step through each Leslie matrix (for pop i)
    output = sim_model(A=Leslie3d[,,f], timesteps=timesteps, 
                       alpha=alphas[i], beta=1000, 
                       sig_r=sig_r, initial_eggs=initial_eggs)
    length(output$Nsize) <- length(output$N_t) #trim ts vector, -2 elements
    output.matrix[,,f] <- do.call(cbind,output) #fill in array for pop i
    #colnames(output.matrix) <- names(output)
  }
  output.3d.list[[i]] <- output.matrix
}
rm(i,f,Leslie3d,output.matrix,output) #clean up

# At this point I have one important object:
# 1. [output.3d.list] a list of 3d arrays. Each array is timeseries output
#    from simulations at different F levels. 


# ------------------------------------------------------------------------- #
# --- Format output ts for plotting simulations using output.3d.list
# ------------------------------------------------------------------------- #
variable_type <- c("Nt","eggs","recruits","Nsize")

# --- reorganize egg timeseries data --- #
var.number <- 2 # eggs
df.list <- as.list(rep(NA,length=length(codNames)))
names(df.list) <- codNames
for (i in 1:length(output.3d.list)) {
  # first, reformat data to work with ggplot
  aa <- as.data.frame(output.3d.list[[i]][,var.number,])
  aa$year <- seq(from=1, to=length(aa[,1]))
  colnames(aa) <- c("1","0.8","0.5","0.35","0.2","year")
  aa1 <- aa %>% gather(FLEP,value,1:5)
  aa1$variable <- rep(variable_type[var.number],length=length(aa1[,1]))
  aa1$codNames <- rep(codNames[i],length=length(aa[,1]))
  df.list[[i]] <- aa1
  rm(aa1,aa)}
eggs.ts <- do.call(rbind,df.list)
rm(df.list) #clean up

# --- reorganize recruit timeseries data --- #
var.number <- 3 # recruits
df.list <- as.list(rep(NA,length=length(codNames)))
names(df.list) <- codNames
for (i in 1:length(output.3d.list)) {
  # first, reformat data to work with ggplot
  aa <- as.data.frame(output.3d.list[[i]][,var.number,])
  aa$year <- seq(from=1, to=length(aa[,1]))
  colnames(aa) <- c("1","0.8","0.5","0.35","0.2","year")
  aa1 <- aa %>% gather(FLEP,value,1:5)
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
  colnames(aa) <- c("1","0.8","0.5","0.35","0.2","year")
  aa1 <- aa %>% gather(FLEP,value,1:5)
  aa1$variable <- rep(variable_type[var.number],length=length(aa1[,1]))
  aa1$codNames <- rep(codNames[i],length=length(aa[,1]))
  df.list[[i]] <- aa1
  rm(aa1,aa)}
nsize.ts <- do.call(rbind,df.list)
rm(df.list) #clean up

ts.data <- rbind(eggs.ts,recruits.ts,nsize.ts) #combine data
# ready for plotting!


# ------------------------------------------------------------------------- #
# --- Now that data is formated, let's plot!
# ------------------------------------------------------------------------- #
head(ts.data)

p <- list()
# plot eggs - one plot per pop, on each plot 5 lines for different F levels
for (i in 1:length(codNames)){
p[[i]] <- ggplot(ts.data[ts.data$variable == "eggs" & ts.data$codNames == codNames[i],], 
       aes(x=year,y=value,color=FLEP)) +
  xlab("year") + ylab("egg production") +
  geom_line() +
  scale_color_brewer(palette = "Reds") +
  ggtitle(paste(codNames[i])) +
  theme_classic()
}
pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/egg_production_for_diff_FLEPs.pdf', width=7, height=10)
do.call(grid.arrange,c(p,ncol=2))
dev.off()






  plot(output$eggs)

  outputL[[i]] = output #save sim output for each pop 
  # setting 'span' - a vector of odd integers to specify the smoothers
  tmp <- ceiling(sqrt(length(1:(timesteps-rm_first_timesteps-1)))) #square root of timeseries length, rounded
  if (tmp %% 2 == 0) {m <- tmp+1} else {m <- tmp} #make it odd, if the square root is even
  m = m * span.multiplier
  # spectral analysis
  sp = spec.pgram(x=output[[2]][rm_first_timesteps:(timesteps-2)], 
                  spans=c(m,m),plot = FALSE)
  spall[,16] = 2*sp$spec #save spec output for plotting pops together, Helen Wearing says to multiply by 2 ***** change 16 --> i
  #plot(x=sp$freq,y=spall[,i],type="l",main=eigentable$codNames[i])
  #legend("topright",c(paste("mode=",eigentable$mode[i]),
  #                    paste("sd (mode)=",round(x=eigentable$sd_mode[i],digits=2)),
  #                    paste("CV=",round(x=eigentable$cvs_mode[i],digits=2))))
}
rm(sp,tmp,A,output,i) #clean up
sp_egg = cbind(seq(from=0.001111111, to=0.5, by=0.001111111),spall) #save freq with spec
colnames(sp_egg) = c("freq",eigentable$codNames)
sp_egg = as.data.frame(sp_egg)
plot(output$eggs)

# ---
# run simulation - RECRUITS (before noise) --> that means output[[3]]
# create empty objects for simulation
freq = seq(from=0.001111111, to=0.5, by=0.001111111) #all frequencies, for plotting later
spall = matrix(NA, length(freq),length(eigentable$codNames)) #set up empty matrix, spec for all plots
outputL = vector("list", length(eigentable$codNames)) #create empty list to store sim output
names(outputL) = eigentable$codNames

for (i in 1:length(eigentable$codNames)) { # step through each cod population
  #A = read.table(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/matrix_maxages/'
  #                            ,eigentable$codNames[i], '.txt', sep=''))
  
  # --- #
  # run this section if 'base' Leslie matricies need to be generated
  # load parms for cod pop i
  source(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_pops/',eigentable$codNames[i], '.r', sep=''))
  # this should load parms: L_inf, K, TEMP, maxage
  out=assemble_Leslie(data=datalist[[i]], maxage=maxage, K=K, L_inf=L_inf, TEMP=TEMP,
                        F.halfmax=0,tknot=0)
  Alist[[i]]=out$A #store Leslie in a list just in case I need it later
  A = out$A #store Leslie as A for later use
  #Alist[[i]][1,] <- Alist[[i]][1,]*k[j] #this line multiplies fecundities by k
  eigenvals1[i] <- extract_first_eigen_value(Alist[[i]]) #calc first eigenvalue
  eigenvals2[i] <- extract_second_eigen_value(Alist[[i]]) #calc second eigenvalue
  eigenvals1.2[i] <- eigenvals2[i] / eigenvals1[i] #calc damp ratio
  # remove pop parms for next loop 
  rm(K,L_inf,maxage)
  # --- #
  
  A = as.matrix(A)
  output = sim_model(A=A, timesteps=timesteps, 
                     alpha=1, beta=beta, #change alpha to alpha not 1
                     sig_r=sig_r, initial_eggs=initial_eggs)  
  outputL[[i]] = output #save sim output for each pop 
  # setting 'span' - a vector of odd integers to specify the smoothers
  tmp <- ceiling(sqrt(length(1:(timesteps-rm_first_timesteps-1)))) #square root of timeseries length, rounded
  if (tmp %% 2 == 0) {m <- tmp+1} else {m <- tmp} #make it odd, if the square root is even
  m = m * span.multiplier
  # plot frequency content
  sp = spec.pgram(x=output[[3]][rm_first_timesteps:(timesteps-2)], 
                  spans=c(m,m),plot = FALSE)
  spall[,i] = 2*sp$spec #save spec output for plotting pops together, Helen Wearing says to multiply by 2
  #plot(x=sp$freq,y=spall[,i],type="l",main=eigentable$codNames[i])
  #legend("topright",c(paste("mode=",eigentable$mode[i]),
  #                    paste("sd (mode)=",round(x=eigentable$sd_mode[i],digits=2)),
  #                    paste("CV=",round(x=eigentable$cvs_mode[i],digits=2))))
}
rm(sp,tmp,A,output,i) #clean up
sp_rec = cbind(seq(from=0.001111111, to=0.5, by=0.001111111),spall) #save freq with spec
colnames(sp_rec) = c("freq",eigentable$codNames)
sp_rec = as.data.frame(sp_rec)

# ===========================================================
# Let's look at some of the output before doing more analysis
# ---
# total population size over time
n_size <- lapply(outputL, function(x) x['Nsize'])
n_sizedf <- do.call(cbind.data.frame, n_size)
colnames(n_sizedf) <- eigentable$codNames
n_sizedf$time <- seq(1:length(n_sizedf$Northsea))
n_sizedf <- subset(n_sizedf, time > 500)
n_size_long <- melt(n_sizedf, id = "time")
n_size_long$variable <- factor(n_size_long$variable)

ggplot(n_size_long, aes(x=time, y=value, group=variable, color=variable, shape=variable)) +
  scale_shape_manual(values=1:nlevels(n_size_long$variable)) +
  geom_point() +
  geom_line() +
  labs(y="N",x="time")
rm(n_size, n_sizedf, n_size_long) # clean up
# ---
# egg production over time (this should be the same as total population size, 
# but one time step off)
eg <- lapply(outputL, function(x) x['eggs'])
egdf <- do.call(cbind.data.frame, eg)
colnames(egdf) <- eigentable$codNames
egdf$time <- seq(1:length(egdf$Northsea))
egdf <- subset(egdf, time > 500)
eg_long <- melt(egdf, id = "time")
eg_long$variable <- factor(eg_long$variable)

ggplot(eg_long, aes(x=time, y=value, group=variable, color=variable, shape=variable)) +
  scale_shape_manual(values=1:nlevels(eg_long$variable)) +
  geom_point() +
  geom_line() +
  labs(y="eggs",x="time")
rm(eg, egdf, eg_long) # clean up
# ---
# recruitment (before noise) over time
rec <- lapply(outputL, function(x) x['recruits'])
recdf <- do.call(cbind.data.frame, rec)
colnames(recdf) <- eigentable$codNames
recdf$time <- seq(1:length(recdf$Northsea))
recdf <- subset(recdf, time > 500)
rec_long <- melt(recdf, id = "time")
rec_long$variable <- factor(rec_long$variable)

ggplot(rec_long, aes(x=time, y=value, group=variable, color=variable, shape=variable)) +
  scale_shape_manual(values=1:nlevels(rec_long$variable)) +
  geom_point() +
  geom_line() +
  labs(y="recruits (before noise)",x="time")
rm(rec, recdf, rec_long) # clean up



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
# plot: spectral analysis after fishing is varied
# multi-panel plot






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



