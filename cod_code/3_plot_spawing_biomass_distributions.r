# Plot spawning biomass distributions
# by: Mikaela Provost
# Goals:
# 1. calculate probability of spawning biomass at age distributions
# 2. create and export eigentables for MG and MP parms

library(ggplot2)
library(gridExtra)
library(tidyr)
library(dplyr)
# ---
# Load functions:
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/2_cod_functions.r")

# ---
# load cod data, break into separate populations
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/0_load_cod_data.r")

# ---
# load max ages, table showing max age in each population
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/6_find_and_export_oldest_fish_in_each_pop.r")
# change max age for 23LKJ to 17
# this needs to be changed because 2J3KL is an outlier if max age is 
# left at 20y. The assessment (Brattey et al. 2010, p28) says that most
# ages are 1-10y, with maximum age reaching 17.
max_ages_table[max_ages_table$codNames == "cod2J3KL", ]$max_ages <- 17
max_ages_table

# ---
# reorder pops by peak spawning age
# load peak spawning age info
#eigentable = read.csv("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLSB/eigentable.csv",
#                      header=TRUE,stringsAsFactors = FALSE)
#eigentable = as.data.frame(eigentable)
#codNames_ordered_by_peak <- eigentable %>% arrange(mode_age) %>% pull(codNames)
#codNames_ordered_by_peak_plot <- eigentable %>% arrange(mode_age) %>% pull(codNames_plot)
# ---
# Plot spawning biomass distribution -- new way: treat as probability distribution
# y axis = probability of spawning
# x axis = age
cvs_modeMG = rep(NA, length=length(codNames))
mode_ageMG = rep(NA, length=length(codNames))
sd_modeMG = rep(NA, length=length(codNames))

cvs_modeMP = rep(NA, length=length(codNames))
mode_ageMP = rep(NA, length=length(codNames))
sd_modeMP = rep(NA, length=length(codNames))

max_ages <- rep(NA, length=length(codNames))
temp <- rep(NA, length=length(codNames))
F.halfmax = 0
codNames_plot <- rep(NA, length=length(codNames))
# note: need to recalculate sd with mode, instead of mean

pMG <- as.list(rep(NA,length=length(codNames)))
names(pMG) <- codNames
pMP <- as.list(rep(NA,length=length(codNames)))
names(pMP) <- codNames
for (i in 1:length(codNames)) { # step through each cod population
 
  # this should load parms: L_inf, K, TEMP, maxage
  source(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_pops/',codNames[i], '.r', sep=''))
  # calculate LEP at each age
  out = calculate_LSB_at_age_by_F(maxage=maxage,L_inf=L_inf,K=K,TEMP=TEMP,
                                       F.halfmax=0,B0=B0,B1=B1)
                                  
  Ages = seq(from=1,to=maxage,by=1)
  codNames_plot[i] <- name
  max_ages[i] <- maxage
  temp[i] <- TEMP
  # calculate probability of spawning at age
  #p_spawnMG = as.data.frame((out$LEP_MG/sum(out$LEP_MG))*(15/sum(out$LEP_MG))*(1/length(out$LEP_MG)))
  p_spawnMG = as.data.frame(out$LEP_MG/sum(out$LEP_MG))
  #p_spawnMG = as.data.frame(out$LEP_MG/sum(out$LEP_MG))
  colnames(p_spawnMG) <- "p_spawnMG"
  #p_spawnMP = as.data.frame((out$LEP_MP/sum(out$LEP_MP))*(15/sum(out$LEP_MP))*(1/length(out$LEP_MP)))
  p_spawnMP = as.data.frame(out$LEP_MP/sum(out$LEP_MP))
  #p_spawnMP = as.data.frame(out$LEP_MP/sum(out$LEP_MP))
  colnames(p_spawnMP) <- "p_spawnMP"
  keep= cbind(p_spawnMG,p_spawnMP,Ages)
  
  # using mode in sd
  mode_ageMG[i] = keep$Age[which.max(keep$p_spawnMG)] # what is the age with highest probability?
  sd_modeMG[i] = sqrt(sum(keep$p_spawnMG*(keep$Age-mode_ageMG[i])^2) ) # stdev
  cvs_modeMG[i] = sd_modeMG[i]/mode_ageMG[i] # coefficient of variation 
  
  mode_ageMP[i] = keep$Age[which.max(keep$p_spawnMP)] # what is the age with highest probability?
  sd_modeMP[i] = sqrt( sum(keep$p_spawnMP*(keep$Age-mode_ageMP[i])^2) ) # stdev
  cvs_modeMP[i] = sd_modeMP[i]/mode_ageMP[i] # coefficient of variation 
  
  # Plot spawning distribution for each population:
  pMG[[i]] <- ggplot(keep,aes(x=Ages,y=p_spawnMG)) +
    geom_line() + theme_classic() + xlab("") + 
    ylab("") + 
    ggtitle(paste(codNames[i])) +
    scale_y_continuous(limits = c(0,0.35)) + #y axis for not adjusted
    xlim(0,20) +
    geom_vline(xintercept=mode_ageMG[i],linetype="dashed") +
    geom_text(x=(mode_ageMG[i]+2), y=0.3, label=mode_ageMG[i], size=4) +
    theme(text = element_text(size = 10))
  
  pMP[[i]] <- ggplot(keep,aes(x=Ages,y=p_spawnMP)) +
    geom_line() + theme_classic() + xlab("Age") + 
    ylab("") + 
    ggtitle(paste(codNames[i])) +
    scale_y_continuous(limits = c(0,0.35)) + #y axis for not adjusted
    xlim(0,20) +
    geom_vline(xintercept=mode_ageMP[i],linetype="dashed") +
    geom_text(x=(mode_ageMP[i]+2), y=0.3, label=mode_ageMP[i], size=4) +
    theme(text = element_text(size = 10))
}

eigentable_MM <- as.data.frame(cbind(codNames,codNames_plot,max_ages,temp,mode_ageMG,mode_ageMP,
                                     sd_modeMG,sd_modeMP,cvs_modeMG,cvs_modeMP))
str(eigentable_MM)
eigentable_MM$max_ages <- as.numeric(levels(eigentable_MM$max_ages))[eigentable_MM$max_ages]
eigentable_MM$temp <- as.numeric(levels(eigentable_MM$temp))[eigentable_MM$temp]
eigentable_MM$mode_ageMG <- as.numeric(levels(eigentable_MM$mode_ageMG))[eigentable_MM$mode_ageMG]
eigentable_MM$mode_ageMP <- as.numeric(levels(eigentable_MM$mode_ageMP))[eigentable_MM$mode_ageMP]
eigentable_MM$sd_modeMG <- as.numeric(levels(eigentable_MM$sd_modeMG))[eigentable_MM$sd_modeMG]
eigentable_MM$sd_modeMP <- as.numeric(levels(eigentable_MM$sd_modeMP))[eigentable_MM$sd_modeMP]
eigentable_MM$cvs_modeMG <- as.numeric(levels(eigentable_MM$cvs_modeMG))[eigentable_MM$cvs_modeMG]
eigentable_MM$cvs_modeMP <- as.numeric(levels(eigentable_MM$cvs_modeMP))[eigentable_MM$cvs_modeMP]

plot(eigentable_MM$cvs_modeMG,eigentable_MM$cvs_modeMP)
# export eigentable_MM
write.csv(eigentable_MM,file='C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLSB/eigentable_MM.csv')
# Export high res fig
tiff(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/manuscript3/SI/figS_spawning_distributions_MG.tiff', units="in", width=7, height=7, res=300)
#do.call(grid.arrange,c(pMG,ncol=4,top="Gislason, LEP-at-age/total LEP, where LEP-at-age = f*survival-to-that-age",left="Pr(spawning)"))
do.call(grid.arrange,c(pMG,ncol=4,bottom="Age (years)",left="Pr(spawning)"))
dev.off()

# tiff(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/manuscript/fig1_spawning_distributions_MP_oneoverLEP.tiff', units="in", width=7, height=7, res=300)
# do.call(grid.arrange,c(pMP,ncol=4,top="Pauly, LEP-at-age/total LEP, where LEP-at-age = f*survival-to-that-age",left="Pr(spawning)"))
# dev.off()

tiff(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/manuscript3/fig3a_spawning_distributions_subplot.tiff', units="in", width=2, height=6, res=300)
do.call(grid.arrange,c(pMG[c(6,3,15,5)],ncol=1))
dev.off()




# ****************************************
# code below: produced older plots, likely needs updating
# ****************************************


# ---
# Plot CV vs damping ratio:
# Read in each Leslie matrix for different cod pops (read in from folder),
# Output from functions is two vectors:
# the magnitude of the first and second eigenvalues for 
# each population.

# what's happening in the for loop?
# first, read in Leslie table for one pop,
# then calculate eigenvalues & store these 
# in firsts and seconds. After the loop, 
# create vector of k values (here, k=1) and 
# combine firsts, seconds, k vector.
firsts = rep(NA,length(codNames))
seconds = rep(NA,length(codNames))
multiplier_values <- read.csv(file='C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/multiplier_values.csv')
multiplier_values <- as.vector(multiplier_values['x'])
multiplier_values <- multiplier_values[['x']]

for(i in 1:length(codNames)){
  A = read.table(paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/k1/",
                       codNames[i],".txt",sep=""))
  AA = A[1,]*multiplier_values[i] # multiply fecundities by 1/lambda
  A[1,] <- AA # replace the first row, will new fecundities
  firsts[i] = extract_first_eigen_value(Lesliematrix=A)
  seconds[i] = extract_second_eigen_value(Lesliematrix=A)
  }

# combine eigenvalues datasets for plotting
firsts <- data.frame(firsts)
seconds <- data.frame(seconds)
eigentable <- cbind(firsts,seconds,cvs)
eigentable$dampratio <- eigentable$firsts/eigentable$seconds
eigentable$temp = c(8.4,7.13,7,7.4,4,11.9,5.8,6.5,9.57,1,8,8,1.75,3.5,0,2.5)
eigentable$codNames = codNames
#eigentable$mean_age = mean_age
#eigentable$sd_age = sd_age
eigentable$mode = mode_age
eigentable$sd_mode = sd_mode
eigentable$cvs_mode = cvs_mode
eigentable$max_ages = max_ages_table$max_ages


# export the eigentable
write.csv(eigentable,file='C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLSB/eigentable5.csv')


ggplot(eigentable, aes(x=cvs_mode,y=dampratio)) +
  geom_point(aes(color=temp), size=3) +  
  geom_text(aes(label=codNames),hjust=-0.1,vjust=-0.2,check_overlap = F) +
  scale_color_gradientn(colors=rev(rainbow(n=10,start=0,end=0.7)))

#eigentable <- eigentable[!(eigentable$codNames == "Celtic"),] # not sure why I removed this earlier.

# ---
# plot A: CV vs lambda1/lambda2
# ---
pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/cv_dampratio_adjLeslie_mode_v3_oldestage.pdf', width=8, height=6.5)
ggplot(eigentable, aes(x=cvs_mode,y=dampratio)) +
  theme_classic() +
  theme(axis.title.y=element_text(size=rel(1))) +
  theme(axis.title.x=element_text(size=rel(1))) +
  geom_point(aes(color=temp), size=3) +
  #scale_color_gradient(low="blue",high="red") +
  scale_color_gradientn(colors=rev(rainbow(n=10,start=0,end=0.7))) +
  geom_text(aes(label=codNames),hjust=-0.1,vjust=-0.2,check_overlap = F) +
  #ylim(0,0.75) +
  #xlim(0.3,0.8) +
  ylab(expression(paste(lambda[1],"/","|",lambda[2],"|"," (damping ratio, ",lambda[1]," adjusted to 1)"))) +
  xlab("CV of spawning distribution (sd(mode)/mode)") +
  ggtitle("spawning age distribution tails truncated at oldest age observed in each population \n(23LKJ at 17y)") +
  theme(plot.title = element_text(hjust = 0.5))
  
dev.off()

# ---
# plot B: mode vs stdev -- multiplot starts here
# ---
pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/multiplot_ModeStdevDamp_v3_oldestage.pdf', width=10, height=8)
a <- ggplot(eigentable, aes(x=mode,y=sd_mode)) +
  theme_classic() +
  theme(axis.title.y=element_text(size=rel(1))) +
  theme(axis.title.x=element_text(size=rel(1))) +
  geom_point(aes(color=temp), size=1) +
  #scale_color_gradient(low="blue",high="red") +
  scale_color_gradientn(colors=rev(rainbow(n=10,start=0,end=0.7))) +
  geom_text(aes(label=codNames),hjust=-0.1,vjust=-0.2,check_overlap = F, size=2) +
  ylim(1,5) +
  xlim(3,9) +
  ylab("stdev about the mode") +
  xlab("mode") 
  #ggtitle("spawning age distribution tails truncated at p(0.01)")


# ---
# plot C: mode vs damping ratio
# ---
b <- ggplot(eigentable, aes(x=mode,y=dampratio)) +
  theme_classic() +
  theme(axis.title.y=element_text(size=rel(1))) +
  theme(axis.title.x=element_text(size=rel(1))) +
  geom_point(aes(color=temp), size=1) +
  #scale_color_gradient(low="blue",high="red") +
  scale_color_gradientn(colors=rev(rainbow(n=10,start=0,end=0.7))) +
  geom_text(aes(label=codNames),hjust=-0.1,vjust=-0.2,check_overlap = F, size=2) +
  #ylim(0,0.75) +
  #xlim(1.2,3) +
  ylab("damping ratio") +
  xlab("mode") 


# ---
# plot D: stdev vs damping ratio
# ---
c <- ggplot(eigentable, aes(x=sd_mode,y=dampratio)) +
  theme_classic() +
  theme(axis.title.y=element_text(size=rel(1))) +
  theme(axis.title.x=element_text(size=rel(1))) +
  geom_point(aes(color=temp), size=1) +
  #scale_color_gradient(low="blue",high="red") +
  scale_color_gradientn(colors=rev(rainbow(n=10,start=0,end=0.7))) +
  geom_text(aes(label=codNames),hjust=-0.1,vjust=-0.2,check_overlap = F, size=2) +
  #ylim(0,0.75) +
  xlim(1,5) +
  ylab("damping ratio") +
  xlab("stdev about the mode") 

# ---
# plot E: max age vs dampratio
# ---
d <- ggplot(eigentable, aes(x=max_ages,y=dampratio)) +
  theme_classic() +
  theme(axis.title.y=element_text(size=rel(1))) +
  theme(axis.title.x=element_text(size=rel(1))) +
  geom_point(aes(color=temp), size=1) +
  #scale_color_gradient(low="blue",high="red") +
  scale_color_gradientn(colors=rev(rainbow(n=10,start=0,end=0.7))) +
  geom_text(aes(label=codNames),hjust=-0.1,vjust=-0.2,check_overlap = F, size=2) +
  #ylim(0,0.75) +
  xlim(7,20) +
  ylab("damping ratio") +
  xlab("max age") 

grid.arrange(a, b, c, d, ncol=2)
dev.off()

# ---
# plot E: regression 
# ---
myfit <- lm(dampratio ~ sd_mode + mode, data=eigentable)
summary(myfit)
plot(myfit)


#https://cran.r-project.org/web/packages/stargazer/stargazer.pdf
#install.packages("stargazer")
library(stargazer)
stargazer(myfit, type="text")
stargazer(myfit)

