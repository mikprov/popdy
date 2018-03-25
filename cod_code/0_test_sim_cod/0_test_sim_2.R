# test code for Loo

# Plan:
# 1. load leslie matrix for Iceland, test case
# 2. define parms
# 3. set timesteps (this can be adjusted)
# 4. initialize abundance matrix
# 5. fill in timesteps in abundance matrix (for loop)
# 6. plot abundance
 
# ---
# 1. load leslie matrix
# ---
# Leslie matrix for Iceland
#A <- read.table("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/matrix_maxages/Iceland.txt")
#A <- as.matrix(A) #fix format to matrix

A <- read.table("Iceland.txt")
A <- as.matrix(A)
# ---
# 2. define parms
# ---
# set up empty matrix for abundance at age at time
maxage = length(A[,1]) #get the max age, uses dimensions from Leslie matrix
ages = length(seq(1,maxage)) #sequence of ages (1,2,3,... maxage)
N0 = c(100, rep(0,ages-1)) #initial age vector, start with 100 at age 1 and 0's for all other ages
alpha = 100000 #for BH
beta = 1000 #for BH
sig_r=0.2 #Sigma multiplies variability term-- I choose 0.3 to start with


# ---
# 3. set timesteps
# ---
timesteps = 100 #using a short timeseries for now, makes it easier to see matrix output (Nt)
               #increase timesteps to get a longer time series


# ---
# 4. initialize abundance matrix (Nt)
# ---
Nt = matrix(0,ages,timesteps) #Initialize matrix to store numbers at age (rows), for all time steps (columns)
Nt[,1] = N0 #Put in initial age vector (100,0,0,...maxage)

Nt[,2] = A %*% Nt[,1]  # multiply initial age vector with Leslie to get 2nd age vector


# ---
# 5. fill in time steps
# ---
# this for loop fills in the Nt matrix. 
# Since we already have age vectors at time 1 and 2,
# the loop will begin filling in ages at time 3 (aka column 3).
set.seed(1) #Set the seed so that every simulation uses same random sequence
for(t in 1:(timesteps-2)){ #step through each time step
  
  # first, replace number of age 1 individuals (#adults) in current age 
  # vector (column 2 to begin with) by pluggin in this value into BH and multiply by noise
  Nt[1,t+1] = ((alpha*Nt[1,t+1])/(1+beta*Nt[1,t+1]))*exp(sig_r*rnorm(1,mean=0,sd=1))
  
  # second, perform population projection for one time step (i.e., column 3 to start)
  Nt[,t+2] = A %*% Nt[,t+1]
}

# ---
# 6. plot
# ---
#Nt #look at age vector matrix, columns=time steps, rows=age
Nt_total_abundance = colSums(Nt) #total population size at each time step
plot(x=(1:timesteps), y=Nt_total_abundance, type="l")


# save workspace
save.image(file="c:/Users/provo/Documents/GitHub/popdy/cod_code/0_test_sim_3.RData")
