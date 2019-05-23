# Plot spawning biomass distributions
# by: Mikaela Provost

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
eigentable = read.csv("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLSB/eigentable.csv",
                      header=TRUE,stringsAsFactors = FALSE)
eigentable = as.data.frame(eigentable)
codNames_ordered_by_peak <- eigentable %>% arrange(mode_age) %>% pull(codNames)
codNames_ordered_by_peak_plot <- eigentable %>% arrange(mode_age) %>% pull(codNames_plot)
# ---
# Plot spawning biomass distribution -- new way: treat as probability distribution
# y axis = probability of spawning
# x axis = age
cvs = rep(NA,length(codNames)) # empty vector to store cv for populations
cvs_mode = rep(NA, length(codNames))
mean_age = rep(NA, length(codNames))
mode_age = rep(NA, length(codNames))
sd_age = rep(NA, length(codNames))
sd_mode = rep(NA, length(codNames))
F.halfmax = 0
# note: need to recalculate sd with mode, instead of mean

p <- list()
names(p) <- codNames_ordered_by_peak
for (i in 1:length(codNames_ordered_by_peak)) { # step through each cod population
 
  # this should load parms: L_inf, K, TEMP, maxage
  source(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_pops/',codNames_ordered_by_peak[i], '.r', sep=''))
  # calculate LEP at each age
  lsb.at.k = calculate_LSB_at_age_by_F(maxage=maxage,L_inf=L_inf,K=K,TEMP=TEMP,
                                       F.halfmax=F.halfmax,B0=B0,B1=B1)
  Ages = seq(from=1,to=length(lsb.at.k[,1]),by=1)
  # calculate probability of spawning at age
  p_spawn = as.data.frame(lsb.at.k[,1] / sum(lsb.at.k[,1])) 
  colnames(p_spawn) <- "p_spawn"
  keep= cbind(p_spawn,Ages)
  
  # using mode in sd
  mode_age[i] = keep$Age[which.max(keep$p_spawn)] # what is the age with highest probability?
  sd_mode[i] = sqrt( sum(keep$p_spawn*(keep$Age-mode_age[i])^2) ) # stdev
  cvs_mode[i] = sd_mode[i]/mode_age[i] # coefficient of variation 
  
  # Plot spawning distribution for each population:
  p[[i]] <- ggplot(keep,aes(x=Ages,y=p_spawn)) +
    geom_line() + theme_classic() + xlab("Age") + 
    ylab("Pr(spawning)") + ggtitle(paste(codNames_ordered_by_peak_plot[i])) +
    scale_x_continuous(limits = c(0,25)) +
    scale_y_continuous(limits = c(0,0.35)) +
    geom_vline(xintercept=mode_age[i],linetype="dashed") +
    geom_text(x=(mode_age[i]+2), y=0.30, label=mode_age[i], size=4) +
    theme(text = element_text(size = 10))
    
}
# Export high res fig
tiff(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/fig1_spawning_distributions.tiff', units="in", width=7, height=7, res=300)
do.call(grid.arrange,c(p,ncol=4))
dev.off()
tiff(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/fig1_spawning_distributions_subplot.tiff', units="in", width=6, height=1.75, res=300)
do.call(grid.arrange,c(p[c(1,9,11,16)],ncol=4))
dev.off()
# Export pdf
pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/SBplot_nofishing_probofspawning_v5.pdf', width=7, height=8)
do.call(grid.arrange,c(p,ncol=4))
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

