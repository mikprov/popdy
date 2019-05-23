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

