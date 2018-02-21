# Plot spawning biomass distributions
# by: Mikaela Provost

library(ggplot2)
library(gridExtra)
# ---
# Load functions:
#  extract_first_eigen_value()
#  extract_second_eigen_value()
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
max_ages_table[max_ages_table$codNames == "2J3KL", ]$max_ages <- 17
max_ages_table

# ---
# Plot spawning biomass distribution -- new way: treat as probability distribution
# y axis = probability of spawning
# x axis = age
Age = 1:40
cvs = rep(NA,length(codNames)) # empty vector to store cv for populations
cvs_mode = rep(NA, length(codNames))
mean_age = rep(NA, length(codNames))
mode_age = rep(NA, length(codNames))
sd_age = rep(NA, length(codNames))
sd_mode = rep(NA, length(codNames))
# note: need to recalculate sd with mode, instead of mean

pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/SBplot_nofishing_probofspawning_v3_oldestage.pdf', width=7, height=10)
par(mfrow=c(5,3))

for (i in 1:length(codNames)) { # step through each cod population
  datax <- read.table(file = 
                        paste('C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLSB/k1/'
                              ,codNames[i], '.txt', sep=''),header=T)
  p_spawn = datax[,1] / sum(datax[,1]) # datax[,1] is LSB at age
                                       # probability of spawning at age = LSB at age/total LSB
  p_table = data.frame(cbind(Age,p_spawn))
  #keep <- p_table[which(p_table$p_spawn > 0.01),] # remove probabilities less than 0.01
  keep <- p_table[which(p_table$Age <= 
                          max_ages_table[max_ages_table$codNames == codNames[i],]$max_ages),
                  ] # removes ages older than the oldest age found in the population
  
  # using mode
  #mode_age[i] = which.max(keep$p_spawn) # mode, testing out instead of mean in CV
  mode_age[i] = keep$Age[which.max(keep$p_spawn)] # what is the age with highest probability?
  sd_mode[i] = sqrt( sum(keep$p_spawn*(keep$Age-mode_age[i])^2) ) # stdev
  cvs_mode[i] = sd_mode[i]/mode_age[i] # coefficient of variation 
  
  # using mean -- needs fixing if using mean
  #mean_age[i] = sum(p_spawn * Age) # mean age = sum (probability of spawning at age * age )
  #sd_age[i] = sqrt( sum( p_spawn*(Age-mean_age[i])^2)) # standard deviation 
  #cvs[i] = sd_age[i]/mean_age[i] # coefficient of variation 
  
  # Plot spawning distribution for each population:
  plot(x=keep$Age, y=keep$p_spawn,type="l",
       main=codNames[i], ylab="probability of spawning",
       ylim=c(0,0.3), xlim=c(1,40))
  abline(v=mode_age[i],col="red",lty=2)
  # paste the mean, sd, cv on each population plot
  legend("topright",c(paste("mode=",round(x=mode_age[i],digits=2)),
                      paste("sd (mode)=",round(x=sd_mode[i],digits=2)),
                      paste("CV=",round(x=cvs_mode[i],digits=2))))
}
par(mfrow=c(1,1))
dev.off()



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

