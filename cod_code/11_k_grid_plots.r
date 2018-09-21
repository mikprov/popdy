# This code creates grid plots: spectral analysis at different k values
# for each cod pop

# Plan
# 1. read in LSBlist = matricies of LSB at age for diff k, for each pop
# 2. use LSBlist to generate Jacobian matricies for simulations
# 3. get output from simulations ready for plotting
# 4. plot 1 - 4 panels, all pops together
# 5. plot 2 - 4 rows (k vals), columns are different peak spawning ages

# --
# load functions, libraries
library(ggplot2)
library(gridExtra)
library(dplyr)
library(plyr)
library(ggrepel)
library(devtools)
library(broom)
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/2_cod_functions.r") # load functions
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/0_load_cod_data.r") # load cod data
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/simulation_model_cod_v3.r") # load sim model
# read in eigentable - I'm using some information from the table
eigentable = read.csv("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLSB/eigentable5.csv",
                      header=TRUE,stringsAsFactors = FALSE)
eigentable = as.data.frame(eigentable)

# ---
# prep output for new Leslie matricies
# Irish gives errors, remove from below
codNames <- c("Northsea","Coas","W_Baltic",
              "Faroe","NE_Arctic","Celtic",
              "Iceland","Kat","W_Scotland",
              "NGulf","GB","GM",
              "cod3NO","cod3M","cod2J3KL",
              "cod3Ps")
ks = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
F.halfmax = 0 #for now, F is 0
LSBlist <- as.list(rep(NA,length(codNames)))
names(LSBlist) <- codNames
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
rm(i,u)



# Goal: create list of arrays with Leslie matricies for each k-population
A3dlist = as.list(rep(NA,length(ks)))
for (i in 1:length(LSBlist)) {
  lsbatk = LSBlist[[i]]
  Alist = as.list(rep(NA,length(ks)))
  
  for (j in 1:length(ks)) { #this assumes I have LSB at the desired k values
    # start with empty matrix
    A = matrix(0,nrow=length(lsbatk[,1]),ncol=length(lsbatk[,1])) #create an empty 'leslie' matrix
    
    # create new Lelise matricies at each k
    A[1,]=lsbatk[,j]/sum(lsbatk[,j]) #fill in relative LSB at age across top row
    A[1,]=A[1,]*ks[j] #multiply relative fecundities by k
    
    for(a in 2:length(A[,1])-1){ #filling in survivals
      A[a+1,a]=1} #fill in survival of 1 on subdiagonal
    
    Alist[[j]] = A
    
  }
  A3d = array(unlist(Alist),dim=c(length(A[,1]),length(A[,1]),length(ks)))
  A3dlist[[i]] = A3d #store array 
}
rm(i,j)

# simulation using Leslie matricies at each k

timesteps = 1000
rm_first_timesteps = 100
alpha = 100
beta = 10000
initial_eggs = 1000
sig_r = 0.1
span.multiplier = 1 # what is this again?


freq = seq(from=0.001111111, to=0.5, by=0.001111111) #all frequencies, for plotting later
output.norm.3dList = as.list(rep(NA,length(codNames)))
output.raw.3dList = as.list(rep(NA,length(codNames)))

for (y in 1:length(A3dlist)) { #walk through each population
  A3d = A3dlist[[y]] #choose the 3d array of matricies, use these to simulate below
  
  outputL.raw = as.list(rep(NA,length(ks))) #create empty list for storing output at each k 
  outputL.norm = as.list(rep(NA,length(ks)))
                         
  for (k in 1:length(ks)) { #step through k values
    
    A = A3d[,,k] #choose the Leslie matrix for this k value
   
    output = sim_model(A=A, timesteps=timesteps, 
                       alpha=alpha, beta=beta, 
                       sig_r=sig_r, initial_eggs=initial_eggs)  
    #outputL.raw[[k]] = output #save sim output for each pop 
    
    # setting 'span' - a vector of odd integers to specify the smoothers
    tmp <- ceiling(sqrt(length(1:(timesteps-rm_first_timesteps-1)))) #square root of timeseries length, rounded
    if (tmp %% 2 == 0) {m <- tmp+1} else {m <- tmp} #make it odd, if the square root is even
    m = m * span.multiplier
    
    # --- do the spectral analysis
    spec_all = matrix(NA,nrow=3,ncol=length(freq))
    
    norm_to_eq_egg = output[[2]][rm_first_timesteps:(timesteps-2)] / mean(output[[2]][rm_first_timesteps:(timesteps-2)])
    spec_egg_normtoeq = spec.pgram(x=norm_to_eq_egg, 
                          spans=c(m,m),plot = FALSE)
    spec_all[1,] <- 2*spec_egg_normtoeq$spec #save spec output for plotting pops together, Helen Wearing says to multiply by 2
    
    norm_to_eq_rec = output[[3]][rm_first_timesteps:(timesteps-2)] / mean(output[[3]][rm_first_timesteps:(timesteps-2)])
    spec_rec_normtoeq = spec.pgram(x=norm_to_eq_rec, 
                          spans=c(m,m),plot = FALSE)
    spec_all[2,] = 2*spec_rec_normtoeq$spec #save spec output for plotting pops together, Helen Wearing says to multiply by 2
    
    norm_to_eq_adu = output[[4]][rm_first_timesteps:(timesteps-2)] / mean(output[[4]][rm_first_timesteps:(timesteps-2)])
    spec_adu_normtoeq = spec.pgram(x=norm_to_eq_adu, 
                          spans=c(m,m),plot = FALSE)
    spec_all[3,] = 2*spec_adu_normtoeq$spec #save spec output for plotting pops together, Helen Wearing says to multiply by 2
    
    outputL.norm[[k]] = spec_all
  }
  output.norm.3d = array(unlist(outputL.norm),dim=c(3,length(freq),length(ks)))
  output.norm.3dList[[y]] = output.norm.3d #store array 
  #output.raw.3d = array(unlist(outputL.raw), dim=c(3,timesteps,length(ks)))
  #output.raw.3dList[[y]] = output.raw.3d
  
  
}

# convert this into a dataframe
output.norm.3dList # for each pop, 3d object: rows=eggs/rec/adu col=spec output stack=diff k 

simpoplist = as.list(rep(NA,length(output.norm.3dList)))
names(simpoplist) <- codNames

for (e in 1:length(output.norm.3dList)) {
  a = output.norm.3dList[[e]] #convert one pop at a time
  
  simlist = as.list(rep(NA,length(ks)))
  
  for (k in 1:length(ks)) { #for each k, convert sim output into df
  
    simout = as.data.frame(t(a[,,k])) #choose simout for each k
    names(simout) = c("eggs","recruits","adults")
    simout$freq = freq
    simlong = melt(simout,id.vars="freq")
    names(simlong) <- c("freq","out.type","value")
    simlong$kval <- rep(ks[k],length=length(freq))
    simlist[[k]] = simlong #store df for each k
  } 
  df <- ldply(simlist, data.frame)
  df$codNames <- rep(codNames[e],length=length(df$freq))
  simpoplist[[e]] = df
}

spec.out.df <- ldply(simpoplist, data.frame)
table(spec.out.df$codNames)

# Plot
ggplot(spec.out.df[spec.out.df$out.type=="recruits",],
       aes(x=freq,y=value)) +
  geom_line() +
  facet_grid(kval~codNames)







# ------ rewriting this loop above ------------ #
# 1. nested for loop: for each k value, step through each population, simulate
# timeseries of eggs, recruits, and spawner abundance
for (j in 1:length(k)) { #step through k values
  
  for (i in 1:length(names(Alist))) { # step through each dataset in datalist
    
    # --- create Jacobian matrix for pop
    # load parms for cod pop i
    source(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_pops/',names(Alist)[i], '.r', sep=''))
    # this should load parms: L_inf, K, TEMP, maxage
    out=assemble_Leslie(data=datalist[[i]], maxage=maxage, K=K, L_inf=L_inf, TEMP=TEMP,
                        F.halfmax=0,tknot=0)
    Alist[[i]]=out$A
    Alist[[i]][1,] <- Alist[[i]][1,]*k[j] #multiply fecundities by k
    # remove pop parms for next loop 
    rm(K,L_inf,maxage)
    
    # --- using Jacobian matrix (A) simulate timeseries
    A = as.matrix(out$A)
    output = sim_model(A=A, timesteps=timesteps, 
                       alpha=alpha, beta=beta, 
                       sig_r=sig_r, initial_eggs=initial_eggs)  
    outputL[[i]] = output #save sim output for each pop 
    # setting 'span' - a vector of odd integers to specify the smoothers
    tmp <- ceiling(sqrt(length(1:(timesteps-rm_first_timesteps-1)))) #square root of timeseries length, rounded
    if (tmp %% 2 == 0) {m <- tmp+1} else {m <- tmp} #make it odd, if the square root is even
    m = m * span.multiplier
    
    # --- do the spectral analysis
    norm_to_eq_egg = output[[2]][rm_first_timesteps:(timesteps-2)] / mean(output[[2]][rm_first_timesteps:(timesteps-2)])
    spec_egg = spec.pgram(x=output[[2]][rm_first_timesteps:(timesteps-2)], 
                    spans=c(m,m),plot = FALSE)
    spec_egg_all[,i] <- 2*spec_egg$spec #save spec output for plotting pops together, Helen Wearing says to multiply by 2
    
    norm_to_eq_rec = output[[3]][rm_first_timesteps:(timesteps-2)] / mean(output[[3]][rm_first_timesteps:(timesteps-2)])
    spec_rec = spec.pgram(x=norm_to_eq_rec, 
                          spans=c(m,m),plot = FALSE)
    spec_rec_all[,i] = 2*spec_rec$spec #save spec output for plotting pops together, Helen Wearing says to multiply by 2
    
    norm_to_eq_adu = output[[4]][rm_first_timesteps:(timesteps-2)] / mean(output[[4]][rm_first_timesteps:(timesteps-2)])
    spec_adu = spec.pgram(x=norm_to_eq_adu, 
                          spans=c(m,m),plot = FALSE)
    spec_adu_all[,i] = 2*spec_adu$spec #save spec output for plotting pops together, Helen Wearing says to multiply by 2
    
    
  }
  
  # --- normalize the variance to pop with largest variance/freq (spec)
  maxs_egg <- apply(spec_egg_all,2,max) #find max variance for each pop
  normalized_egg = t(t(spec_egg_all)/maxs_egg) #divide var at all freq by max
  
  maxs_rec <- apply(spec_rec_all,2,max) #find max variance for each pop
  normalized_rec = t(t(spec_rec_all)/maxs_rec) #divide var at all freq by max
  
  maxs_adu <- apply(spec_adu_all,2,max)
  normalized_adu = t(t(spec_adu_all)/maxs_adu) #divide var at all freq by max
  
  
  # --- assemble df (with spec and normalized dfs): 
  spec_egg_all <- as.data.frame(spec_egg_all)
  colnames(spec_egg_all) = c(eigentable$codNames)
  spec_egg_all$freq = freq #save freq with spec
  normalized_egg <- as.data.frame(normalized_egg)
  colnames(normalized_egg) = c(eigentable$codNames)
  normalized_egg$freq = freq #save freq with spec
  
  spec_rec_all <- as.data.frame(spec_rec_all)
  colnames(spec_rec_all) = c(eigentable$codNames)
  spec_rec_all$freq = freq #save freq with spec
  normalized_rec <- as.data.frame(normalized_rec)
  colnames(normalized_rec) = c(eigentable$codNames)
  normalized_rec$freq = freq #save freq with spec
  
  spec_adu_all <- as.data.frame(spec_adu_all)
  colnames(spec_adu_all) = c(eigentable$codNames)
  spec_adu_all$freq = freq #save freq with spec
  normalized_adu = as.data.frame(normalized_adu)
  colnames(normalized_adu) = c(eigentable$codNames)
  normalized_adu$freq = freq #save freq with spec
  
  # --- add these cols to df
  add_these_cols <- as.data.frame(subset(eigentable, select=c("codNames","mode","temp")))
  
  
  # --- melt spec and normalized dfs & add extra columns
  # eggs
  p <- melt(spec_egg_all,id="freq")
  colnames(p) <- c("freq",'codNames','spec')
  pp <- melt(normalized_egg, id="freq")
  colnames(pp) <- c('freq','codNames','normalized')
  long <- p %>% left_join(pp, by=c('freq','codNames'))
  
  splong_egg <- long %>% 
    left_join(add_these_cols,by='codNames') %>%
    mutate(freq.adjust = freq*mode) 
  splong_egg$kval <- rep(k[j],length(splong_egg[,1]))
  rm(p,pp,long)
  
  # recruits
  p <- melt(spec_rec_all,id="freq")
  colnames(p) <- c("freq",'codNames','spec')
  pp <- melt(normalized_rec, id="freq")
  colnames(pp) <- c('freq','codNames','normalized')
  long <- p %>% left_join(pp, by=c('freq','codNames'))
  
  splong_rec <- long %>% 
    left_join(add_these_cols) %>%
    mutate(freq.adjust = freq*mode) 
  splong_rec$kval <- rep(k[j],length(splong_rec[,1]))
  rm(p,pp,long)
  
  # spawn adults
  p <- melt(spec_adu_all,id="freq")
  colnames(p) <- c("freq",'codNames','spec')
  pp <- melt(normalized_adu, id="freq")
  colnames(pp) <- c('freq','codNames','normalized')
  long <- p %>% left_join(pp, by=c('freq','codNames'))
  
  splong_adu <- long %>% 
    left_join(add_these_cols) %>%
    mutate(freq.adjust = freq*mode) 
  splong_adu$kval <- rep(k[j],length(splong_adu[,1]))
  rm(p,pp,long)
  
  # store dfs in list 
  kdfs_egg[[j]] <- splong_egg
  kdfs_rec[[j]] <- splong_rec
  kdfs_adu[[j]] <- splong_adu
  
}

# 3. get output from simulations ready for plotting
# combine dfs associated with diff k values
df_egg <- ldply(kdfs_egg, data.frame)
df_rec <- ldply(kdfs_rec, data.frame)
df_adu <- ldply(kdfs_adu, data.frame)

# weird cleaning up... why extra rows? codNames=NA ??
df_egg1 <- df_egg[!is.na(df_egg$codNames),]
df_rec1 <- df_rec[!is.na(df_rec$codNames),]
df_adu1 <- df_adu[!is.na(df_adu$codNames),]




# 4. plot 1: 4 panels, all pops together

