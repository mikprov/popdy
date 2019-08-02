# New simulation method
# Instead of using the BH model, this code simulates
# populations using the Jacobian matrix as if it is 
# a Leslie. We use the Jacobian instead of a Leslie
# because we are evaluating the dynamics about the 
# equilibrium point. 

# Plan:
# 1. create fake Jacobian matrix
# 2. can I simulate displacement from equilibrium when noise is added?


set.seed(2)


toprow <- c(0.1,0.1,0.2,0.4,0.6,0.9)
toprow <- toprow * 0.01
sum(toprow) #should equal 1
length(toprow) #num age classes

jac <- matrix(0,nrow=length(toprow),ncol=length(toprow))
jac[1,] <- toprow
jac[2,1]<-1
jac[3,2]<-1
jac[4,3]<-1
jac[5,4]<-1
jac[6,5]<-1
jac
extract_first_eigen_value(jac) #first eigenvalue is 1

tsteps = 1000
sig_r = 0.1

t0 <- matrix(100,nrow=length(jac[,1]),ncol=1)
tseries <- matrix(NA,nrow=length(jac[,1]),ncol=tsteps)
tseries[,1] <- t0
age1_before_noise <- rep(NA,length=tsteps)
age1_after_noise <- rep(NA,length=tsteps)
set.seed(2)  
for (t in 1:(tsteps-1)){
  step <- jac %*% tseries[,t] #calc the first age vector
  age1_before_noise[t] <- step[1,1]
  age1_after_noise[t] <- step[1,1]*sig_r*exp(rnorm(1,mean=0,sd=1)) #add some noise to age 1
  step[1,1] <- age1_after_noise[t]
  tseries[,t+1] <- step #store age vector
}
plot(colSums(tseries))
plot(age1_after_noise)
plot(age1_before_noise)
plot(exp(sig_r*rnorm(1000,mean=0,sd=1)))




# ****************






# setting 'span' - a vector of odd integers to specify the smoothers
tmp <- ceiling(sqrt(length(1:(tsteps-rm_first_timesteps-1)))) #square root of timeseries length, rounded
if (tmp %% 2 == 0) {m <- tmp+1} else {m <- tmp} #make it odd, if the square root is even
m = m * span.multiplier


out <- spec.pgram(nsize[200:1000],spans=c(m,m),plot = TRUE)
plot(nsize[200:tsteps])
plot(x=out$freq,y=out$spec)
plot(etseries)
