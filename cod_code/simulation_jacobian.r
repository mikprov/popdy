# New simulation method
# Instead of using the BH model, this code simulates
# populations using the Jacobian matrix as if it is 
# a Leslie. We use the Jacobian instead of a Leslie
# because we are evaluating the dynamics about the 
# equilibrium point. 

names(A3dlist) <- codNames

set.seed(32)

tsteps = 1000
sig_r = 2
p=8
t=5

for (p in 1:length(A3dlist)) {
  jac <- A3dlist[[p]][,,8] #choose jac when k=0.9
  jac1 <- jac[1,] / 0.9
  jac[1,] <- jac1
  t0 <- matrix(100,nrow=length(jac[,1]),ncol=1)
  tseries <- matrix(NA,nrow=length(jac[,1]),ncol=tsteps)
  tseries[,1] <- t0
  
  for (t in 1:(tsteps-1)){
    step <- jac %*% tseries[,t] #calc the first age vector
    step[1,1] <- step[1,1]*sig_r*rnorm(1,mean=0,sd=1) #add some noise to age 1
    tseries[,t+1] <- step #store age vector
  }
  plot(colSums(tseries))
  plot(tseries[1,])
  matplot(tseries[,1:10])
}
nsize <- colSums(tseries)



# ****************






# setting 'span' - a vector of odd integers to specify the smoothers
tmp <- ceiling(sqrt(length(1:(tsteps-rm_first_timesteps-1)))) #square root of timeseries length, rounded
if (tmp %% 2 == 0) {m <- tmp+1} else {m <- tmp} #make it odd, if the square root is even
m = m * span.multiplier


out <- spec.pgram(nsize[200:1000],spans=c(m,m),plot = TRUE)
plot(nsize[200:tsteps])
plot(x=out$freq,y=out$spec)
plot(etseries)
