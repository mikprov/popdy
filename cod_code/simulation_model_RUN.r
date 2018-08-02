# Run the simulation model
# ===================================================================

# ---
# load the simulation model
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/simulation_model_cod_v3.r")

# load functions
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/2_cod_functions.r")

# load cod data, break into separate populations
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/0_load_cod_data.r")

# ---
# read in eigentable - I'm using some information from the table
eigentable = read.csv("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLSB/eigentable5.csv",
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
alpha = 100
beta = 10000
initial_eggs = 1000
sig_r = 0.1
span.multiplier = 1 # what is this again?
Alist = as.list(rep(NA,length(datalist))) # store Leslie matrix
names(Alist) = codNames
eigenvals1 = rep(NA,length(codNames))
eigenvals2 = rep(NA,length(codNames))
eigenvals1.2 = rep(NA,length(codNames))

# ---
# run simulation - EGGS --> that means output[[2]]
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
                     alpha=alpha, beta=beta, 
                     sig_r=sig_r, initial_eggs=initial_eggs)  
  outputL[[i]] = output #save sim output for each pop 
  # setting 'span' - a vector of odd integers to specify the smoothers
  tmp <- ceiling(sqrt(length(1:(timesteps-rm_first_timesteps-1)))) #square root of timeseries length, rounded
  if (tmp %% 2 == 0) {m <- tmp+1} else {m <- tmp} #make it odd, if the square root is even
  m = m * span.multiplier
  # plot frequency content
  sp = spec.pgram(x=output[[2]][rm_first_timesteps:(timesteps-2)], 
                  spans=c(m,m),plot = FALSE)
  spall[,i] = 2*sp$spec #save spec output for plotting pops together, Helen Wearing says to multiply by 2
  #plot(x=sp$freq,y=spall[,i],type="l",main=eigentable$codNames[i])
  #legend("topright",c(paste("mode=",eigentable$mode[i]),
  #                    paste("sd (mode)=",round(x=eigentable$sd_mode[i],digits=2)),
  #                    paste("CV=",round(x=eigentable$cvs_mode[i],digits=2))))
}
rm(sp,tmp,A,output,i) #clean up
sp_egg = cbind(seq(from=0.001111111, to=0.5, by=0.001111111),spall) #save freq with spec
colnames(sp_egg) = c("freq",eigentable$codNames)
sp_egg = as.data.frame(sp_egg)

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
                     alpha=alpha, beta=beta, 
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



