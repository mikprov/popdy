# North Atlantic Cod Populations - simulate timeseries
# by: Mikaela Provost

# Plan
# 

#library(foreach)
rm(list=ls())
# ===================================================================
# 1) specify parms 

# If I add populations, I need to specify their names in codPop
# And I need to make sure a script is saved in [path] and is labeled w/the name in codPop
codPop <- c("Celtic_Sea","3NO","North_Sea","West_Scotland",
            "West_Baltic","Kattegat","Coastal_Norway",
            "Faroes","Iceland","NE_Arctic","Georges_Bank",
            "Gulf_of_Maine","N_Gulf_St_Lawrence")

# ===================================================================
# 2) define function
library(ggplot2)
library(reshape2)
library(RColorBrewer)

sim_model <- function(A,timesteps,alpha,beta,sig_r,initial_eggs) {
  
  maxage = length(A[,1])
  ages = length(seq(1,maxage))
  N0 = c(initial_eggs, rep(0,ages-1))
  #sig_r=0.1 #Sigma multiplies variability term
  set.seed(2) #Set the seed so that every simulation uses same random sequence
  
  Nt = matrix(0,ages,timesteps) 
    #Initialize vector of population sizes with extra 
    #rows for egg production (top value in age vector) & 
    #recruitment before variability (output from BH)
  Nt[,1] = N0 #Put in initial values
  Nt[,2] = A %*% Nt[,1]  # multiply initial age vector with Leslie to get 2nd age vector
  eggs = c(Nt[1,1], rep(NA,timesteps-1)) #will save egg production here, this will be the input in BH
  recruits = c(initial_eggs, rep(0, timesteps-1)) #will save recruits here (output from BH)
 
    for(t in 1:(timesteps-2)){ #step through time
      eggs[t+1] = Nt[1,t+1] 
        #Save egg production, this is new number of age 1 individuals
      recruits[t+1] = eggs[t+1]/( (1/alpha) + (eggs[t+1]/beta) ) #((alpha*eggs[t+1])/(1+beta*eggs[t+1])) 
        #save recruits, treat egg production as new spawners in BH
      Nt[1,t+1] = (recruits[t+1])*exp(sig_r*rnorm(1,mean=0,sd=1)) 
        #replace age 1 with recruits from BH, add noise
      Nt[,t+2] = A %*% Nt[,t+1] 
        #perform population projection for one time step
    }
  Nsize = colSums(Nt)
  N_t = Nt[1,][2:(timesteps-1)] #Nt is number of adults
  eggs = eggs[2:(timesteps-1)]
  recruits = recruits[2:(timesteps-1)]
  
  return(list(N_t=N_t, eggs=eggs, recruits=recruits, Nsize=Nsize))
}

# read in eigentable - I'm using some information from the table
# to include on frequency content plots
eigentable = read.csv("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLSB/eigentable5.csv",
                      header=TRUE,stringsAsFactors = FALSE)
eigentable = as.data.frame(eigentable)


# Create multi-panel plot of frequency content plots
pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/freq_content_plots_N.pdf', width=7, height=10)
par(mfrow=c(5,3))
freq = seq(from=0.001111111, to=0.5, by=0.001111111) #all frequencies, for plotting later
spall = matrix(NA, length(freq),length(eigentable$codNames)) #set up empty matrix, for plotting all pops together

# read in each Leslie matrix, simulate, then plot frequency content
for (i in 1:length(eigentable$codNames)) { # step through each cod population
  A = read.table(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/matrix_maxages/'
                              ,eigentable$codNames[i], '.txt', sep=''))
  A = as.matrix(A)
  # run simulation 
  timesteps = 1000
  rm_first_timesteps = 100
  alpha = 100
  beta = 1000
  initial_eggs = 1000
  sig_r = 0.1
  output = sim_model(A=A, timesteps=timesteps, 
                     alpha=alpha, beta=beta, 
                     sig_r=sig_r, initial_eggs=initial_eggs)  
  # setting the 'span' - a vector of odd integers to specify the smoothers
  tmp <- ceiling(sqrt(length(1:(timesteps-rm_first_timesteps-1)))) #square root of timeseries length, rounded
  if (tmp %% 2 == 0) {m <- tmp+1} else {m <- tmp} #make it odd, if the square root is even
  
  # plot frequency content
  sp = spec.pgram(x=output[[1]][rm_first_timesteps:(timesteps-2)], 
             spans=c(m,m),plot = FALSE)
  spall[,i] = sp$spec #save spec output for plotting pops together
  plot(x=sp$freq,y=sp$spec,type="l",main=eigentable$codNames[i])
  legend("topright",c(paste("mode=",eigentable$mode[i]),
                      paste("sd (mode)=",round(x=eigentable$sd_mode[i],digits=2)),
                      paste("CV=",round(x=eigentable$cvs_mode[i],digits=2))))
}
dev.off()
par(mfrow=c(1,1))

# plot freq content on one figure
spall = cbind(seq(from=0.001111111, to=0.5, by=0.001111111),spall)
colnames(spall) = c("freq",eigentable$codNames)
splong = melt(as.data.frame(spall),id="freq") #long form for ggplot
mycolors = c(brewer.pal(name="Set1", n = 8), brewer.pal(name="Paired", n = 8)) #16 pops, but pallettes only have 8

pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/freq_content_oneplot_N.pdf', width=7, height=10)
ggplot(data=splong, aes(x=freq,y=value,color=variable,linetype=variable)) + 
  geom_line(lwd=1) + theme_classic() +
  scale_linetype_manual(values=c(rep("solid",10),rep("dashed",6))) +
  scale_color_manual(values=mycolors) +
  ggtitle("N")
dev.off()


# plot timeseries, can be eggs, recruits, or N
pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/timeseries1.pdf', width=7, height=10)
par(mfrow=c(5,3))
tspt = matrix(NA,nrow=899,ncol=length(eigentable$codNames))
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
plot(x=1:timesteps,y=output[[2]],type="l",main="eggs")
plot(x=1:timesteps,y=output[[3]],type="l",main="recruits (before variability)")
plot(x=1:timesteps,y=output[[4]],type="l",main="N")
plot(x=1:timesteps,y=output[[1]][1,],type="l",main="recruits (after variability)")
par(mfrow = c(1,1))
 
output$Nt
output$recruits[1:20]
plot(output$eggs[100:999],type="l")
matplot(t(output$Nt[,1:20]), type="l")



