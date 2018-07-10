# Spectral analysis for different F values
# April 20, 2018

# Plan
# 1. for each population:
# 2. select the F values for that population (target_Fs)
# 3. generate new Leslie matricies for every F value.
# 4. run sim_model on each Leslie matrix for different F values, generates time series
# 5. run spectral analysis to get spec vs freq: store spec for different Fs in one df 
# 6. repeat the above steps for each population, store spec at F in an array
# 7. query the array to plot spectral analysis for each population at different Fs. 
     #use df that goes along depth of array: one pop, specs at different Fs.

# load packages
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)
library(reshape2)
library(gridExtra)

# load functions
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/2_cod_functions.r")
# assemble_Leslie()
# extract_first_eigen_value()
# calc_LSB_at_age_by_F()
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/0_load_cod_data.r")
# this loads 'datalist' which has the data for each population
# combined into one list

# read in eigentable - I'm using some information from the table
eigentable = read.csv("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLSB/eigentable5.csv",
                      header=TRUE,stringsAsFactors = FALSE)
eigentable = as.data.frame(eigentable)
codNames <- eigentable$codNames


# ---
# 2. define F values of interest
# ---
# we set the values of FLEP constant, but this means there are different values of 
# F that correspond to each of the FLEP values for different populations. 
# read in df of different FLEPs associated with F values
FLEP.F <- read.csv('C:/Users/provo/Documents/GitHub/popdy/cod_code/FLEP.F.csv')
FLEP.F <- FLEP.F[,-1]
head(FLEP.F)
#colnames(FLEP.F) <- c("F.halfmax",codNames)
# these are the FLEP values we are using:
# 1, 0.64, 0.46, 0.28, 0.145
# they correspond to EI values of:
# 0, 0.4, 0.6, 0.8, 0.95

# find the nearest F values that correspond to the desired FLEP value
target_Fs <- rbind(
  
  as.data.frame(melt(FLEP.F,id="F.halfmax") %>% 
                  group_by(variable) %>%
                  arrange(abs(value-1)) %>% #desired FLEP=1
                  slice(1) %>%
                  #mutate(value=round(value,2)) %>%
                  mutate(FLEPlevel=rep(1))),
  
  as.data.frame(melt(FLEP.F,id="F.halfmax") %>% 
                  group_by(variable) %>%
                  arrange(abs(value-0.64)) %>% #desired FLEP=0.64
                  slice(1) %>%
                  #mutate(value=round(value,2)) %>%
                  mutate(FLEPlevel=rep(0.64))),
  
  as.data.frame(melt(FLEP.F,id="F.halfmax") %>% 
                  group_by(variable) %>%
                  arrange(abs(value-0.46)) %>% #desired FLEP=0.46
                  slice(1) %>%
                  #mutate(value=round(value,2)) %>%
                  mutate(FLEPlevel=rep(0.46))),
  
  as.data.frame(melt(FLEP.F,id="F.halfmax") %>% 
                  group_by(variable) %>%
                  arrange(abs(value-0.28)) %>% #desired FLEP=0.28
                  slice(1) %>%
                  #mutate(value=round(value,2)) %>%
                  mutate(FLEPlevel=rep(0.28))),
  
  as.data.frame(melt(FLEP.F,id="F.halfmax") %>% 
                  group_by(variable) %>%
                  arrange(abs(value-0.145)) %>% #desired FLEP=0.145
                  slice(1) %>%
                  #mutate(value=round(value,2)) %>%
                  mutate(FLEPlevel=rep(0.145)))
  
)
target_Fs <- melt(FLEP.F, id="F.halfmax")

# ---
# 3. generate new Leslie matricies for every F value for a single population. 
# ---

# define parms for sim_model()
timesteps = 1000
rm_first_timesteps = 100
alpha = 100
beta = 1000
initial_eggs = 1000
tknot = 0
k = 1
sig_r = 0.1
span.multiplier = 1
freq = seq(from=0.001111111, to=0.5, by=0.001111111) #all frequencies, for plotting later
numFs = length(target_Fs[target_Fs$variable == codNames[1],]$F.halfmax)


# prep for spectral analysis
# setting 'span' - a vector of odd integers to specify the smoothers
tmp <- ceiling(sqrt(length(1:(timesteps-rm_first_timesteps-1)))) #square root of timeseries length, rounded
if (tmp %% 2 == 0) {m <- tmp+1} else {m <- tmp} #make it odd, if the square root is even
m = m * span.multiplier

# create empty objects for loop
Alist_F = vector("list", numFs)
outputL_F = vector("list", numFs)
spec_F = as.data.frame(matrix(NA, nrow = length(freq), ncol = numFs))
NAvector = rep(NA,length=length(freq)*numFs*length(codNames))
speclist = as.list(rep(NA,length(codNames)))
names(speclist) = codNames


for (i in 1:length(codNames)) {# for each population 
  
  for (j in 1:length(target_Fs[target_Fs$variable == codNames[i],]$F.halfmax)) { # step through the different F values
    
    # for some reason 
    
    # pull out the F values for this pop
    Fs <- target_Fs[target_Fs$variable == codNames[i],]$F.halfmax
    
    # load parms for cod pop i: L_inf, K, TEMP, maxage
    source(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_pops/',codNames[i], '.r', sep=''))
    
    # assemble Leslie for one population, at diff F values
    aL = assemble_Leslie(data=datalist[[i]], littlek=k, maxage=maxage, K=K, L_inf=L_inf, TEMP=TEMP, 
                         F.halfmax = Fs[j], tknot=tknot)
    Alist_F[[j]] = aL$A # store the Leslie matricies for diff F values in a list, Alist_F
    
    # DEFINE PARMS OUTSIDE store output from sim_model in a list for diff F values
    outputL_F[[j]] = sim_model(A=Alist_F[[j]],timesteps,alpha,beta,sig_r,initial_eggs) 
    
    # do spectral analysis
    sp = spec.pgram(x=outputL_F[[j]]$eggs[rm_first_timesteps:(timesteps-2)], 
                    spans=c(m,m),plot = FALSE)
    # ===================================================================
    # normalize the variance at frequency
    #spp <- cbind( freq,(sp$spec/max(sp$spec))) #divide var at all freq by max
    spp <- cbind( freq,sp$spec) #use this to not normalize var
    colnames(spp) = c("freq",eigentable$codNames[i])
    spp = as.data.frame(spp) #use this dataframe to plot 
    # ===================================================================
    # 6) multiply the freq values for each population by the mode (shifts the frequency  
    # plot so that cohort bump lines up at freq=1, in theory this should work)
    splong <- melt(spp,id="freq")
    colnames(splong) <- c("freq","codNames","spec")
    spec_F[,j] = splong$spec #save spec output 
    colnames(spec_F)[j] <- paste("F",Fs[j],sep="_")
    
  }
  
  speclist[[i]] = spec_F #store df of spec at diff Fs in array

}
rm(i,j,K,L_inf,TEMP,sp,aL,Fs,splong,spp,spec_F) #clean up
# ---- add names to columns

head(speclist)
# ---
# plot spectral analysis for different F values, one plot for each pop
# ---
mode_table <- as.data.frame(subset(eigentable, select=c("codNames","mode")))
mode_table <- mode_table[order(mode_table$mode),]
p <- list()
for (i in 1:length(codNames)){
  df <- speclist[[i]] # freq vs. spec df for one pop
  df = as.data.frame(cbind(freq, df))
  dflong = melt(df, id="freq")
  colnames(dflong) <- c("freq","Fishing_Level","spec")
  dflong.m <-
    dflong %>% mutate(mode=rep(mode_table[mode_table$codNames==codNames[i],]$mode, length(freq)))
  
  p[[i]] <- ggplot(dflong.m,aes(x=freq,y=spec,color=Fishing_Level)) +
    ylab("spec") +
    xlab("frequency") +
    geom_line() +
    scale_color_brewer(palette = "Reds") +
    scale_y_log10() + #limits=c(1e-03,1)) +
    ggtitle(paste(codNames[i],mode_table[mode_table$codNames==codNames[i],]$mode,sep="_")) +
    theme_classic() +
    geom_vline(xintercept=(1/dflong.m$mode[1]),linetype="dashed")
}
names(p) <- codNames
p <- p[mode_table$codNames]

pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/spectral_analysis_for_diff_Fs_multiplot_log.pdf', width=7, height=15)
do.call(grid.arrange,c(p,ncol=2))
dev.off()

# plot CVs for all F values -----

coef.variation <- function(x) {
  sqrt(var(x))/mean(x)
}

eggproduction <- matrix(NA,nrow=length(FLEP.F$F.halfmax),ncol=length(codNames))
colnames(eggproduction) <- codNames
numFs = length(FLEP.F$F.halfmax)
Alist_F = vector("list", numFs)
for (i in 1:(length(codNames))) {# for each population 
  
  for (j in 1:length(FLEP.F$F.halfmax)) { # step through the different F values
    
    # pull out the F values for this pop
    Fs <- FLEP.F[,i+1]
    
    # load parms for cod pop i: L_inf, K, TEMP, maxage
    source(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_pops/',codNames[i], '.r', sep=''))
    
    # assemble Leslie for one population, at diff F values
    aL = assemble_Leslie(data=datalist[[i]], 
                         littlek=k, maxage=maxage, K=K, L_inf=L_inf, TEMP=TEMP, 
                         F.halfmax = Fs[j], tknot=tknot)
    Alist_F[[j]] = aL$A # store the Leslie matricies for diff F values in a list, Alist_F
    
    # DEFINE PARMS OUTSIDE store output from sim_model in a list for diff F values
    output = sim_model(A=Alist_F[[j]],timesteps,alpha,beta,sig_r,initial_eggs) 
    
    eggproduction[j,i] <- coef.variation(output$eggs[200:length(output$eggs)])
    
  }
}
FLEP <- FLEP.F$F.halfmax
edf <- cbind(as.data.frame(eggproduction),FLEP)
head(edf)
edf_long <- melt(edf,id="FLEP")
colnames(edf_long) <- c("FLEP","codNames","CV")
mode_table <- as.data.frame(subset(eigentable, select=c("codNames","mode","temp","cvs_mode","max_ages")))

edf_long <- edf_long %>% 
  left_join(mode_table) 
head(edf_long)

eggs_Fs_temp <- ggplot(data=edf_long, aes(x=FLEP,y=CV,group=codNames)) +
  geom_line(aes(color=temp)) +
  scale_color_gradientn(colors=rev(rainbow(n=16,start=0,end=0.7))) +
  coord_cartesian(xlim = c(min(edf_long$FLEP), max(edf_long$FLEP) + 1)) +
  geom_text_repel(
    data = subset(edf_long, FLEP == max(FLEP)),
    aes(label = codNames, color = temp),
    size = 3,
    nudge_x = 3,
    #nudge_y = 0.5,
    show.legend = TRUE,
    segment.color = "gray",
    segment.alpha = 0.5
  ) + ggtitle("eggs") + theme_classic()

eggs_Fs_mode <- ggplot(data=edf_long, aes(x=FLEP,y=CV,group=codNames)) +
  geom_line(aes(color=max_ages)) +
  # scale_color_gradientn(colors=rev(rainbow(n=16,start=0,end=0.7))) +
  scale_color_gradient(low="black",high="red")+
  coord_cartesian(xlim = c(min(edf_long$FLEP), max(edf_long$FLEP) + 1)) +
  geom_text_repel(
    data = subset(edf_long, FLEP == max(FLEP)),
    aes(label = codNames, color = max_ages),
    size = 3,
    nudge_x = 3,
    #nudge_y = 0.5,
    show.legend = TRUE,
    segment.color = "gray",
    segment.alpha = 0.5
  ) + ggtitle("eggs") + theme_classic()

# export eggs: FLEP vs CV
plist <- list(eggs_Fs_mode,eggs_Fs_temp)
pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/CV_vs_FLEP_eggs_mode&temp.pdf', width=7, height=15)
do.call("grid.arrange",c(plist,ncol=1))
dev.off()
rm(plist)


# recruits
# load the simulation model
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/simulation_model_cod_v3.r")

recruits_atF <- matrix(NA,nrow=length(FLEP.F$F.halfmax),ncol=length(codNames))
colnames(recruits_atF) <- codNames
numFs = length(FLEP.F$F.halfmax)
Alist_F = vector("list", numFs)
for (i in 1:(length(codNames))) {# for each population 
  
  for (j in 1:length(FLEP.F$F.halfmax)) { # step through the different F values
    
    # pull out the F values for this pop
    Fs <- FLEP.F[,i+1]
    
    # load parms for cod pop i: L_inf, K, TEMP, maxage
    source(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_pops/',codNames[i], '.r', sep=''))
    
    # assemble Leslie for one population, at diff F values
    aL = assemble_Leslie(data=datalist[[i]], 
                         littlek=k, maxage=maxage, K=K, L_inf=L_inf, TEMP=TEMP, 
                         F.halfmax = Fs[j], tknot=tknot)
    Alist_F[[j]] = aL$A # store the Leslie matricies for diff F values in a list, Alist_F
    
    # DEFINE PARMS OUTSIDE store output from sim_model in a list for diff F values
    output = sim_model(A=Alist_F[[j]],timesteps,alpha,beta,sig_r,initial_eggs) 
    
    recruits_atF[j,i] <- coef.variation(output$recruits[200:length(output$recruits)])
    
  }
}
FLEP <- FLEP.F$F.halfmax
rdf <- cbind(as.data.frame(recruits_atF),FLEP)
head(rdf)
rdf_long <- melt(rdf,id="FLEP")
colnames(rdf_long) <- c("FLEP","codNames","CV")
mode_table <- as.data.frame(subset(eigentable, select=c("codNames","mode","temp")))

rdf_long <- rdf_long %>% 
  left_join(mode_table) 
head(rdf_long)

recruits_Fs_temp <- ggplot(data=rdf_long, aes(x=FLEP,y=CV,group=codNames)) +
  geom_line(aes(color=temp)) +
  scale_color_gradientn(colors=rev(rainbow(n=16,start=0,end=0.7))) +
  coord_cartesian(xlim = c(min(rdf_long$FLEP), max(rdf_long$FLEP) + 1)) +
  geom_text_repel(
    data = subset(rdf_long, FLEP == max(FLEP)),
    aes(label = codNames, color = temp),
    size = 3,
    nudge_x = 3,
    #nudge_y = 0.5,
    show.legend = TRUE,
    segment.color = "gray",
    segment.alpha = 0.5
  ) + ggtitle("recruits") + theme_classic()

recruits_Fs_mode <- ggplot(data=rdf_long, aes(x=FLEP,y=CV,group=codNames)) +
  geom_line(aes(color=mode)) +
 # scale_color_gradientn(colors=rev(rainbow(n=16,start=0,end=0.7))) +
  coord_cartesian(xlim = c(min(rdf_long$FLEP), max(rdf_long$FLEP) + 1)) +
  geom_text_repel(
    data = subset(rdf_long, FLEP == max(FLEP)),
    aes(label = codNames, color = mode),
    size = 3,
    nudge_x = 3,
    #nudge_y = 0.5,
    show.legend = TRUE,
    segment.color = "gray",
    segment.alpha = 0.5
  ) + ggtitle("recruits") + theme_classic()


# export recruits: FLEP vs CV
plist <- list(recruits_Fs_mode,recruits_Fs_temp)
pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/CV_vs_FLEP_recruits_mode&temp.pdf', width=7, height=15)
do.call("grid.arrange",c(plist,ncol=1))
dev.off()
