---
title: "Simulation model"
author: "Mikaela Provost"
date: "Sept 19, 2017"
output: html_document
---

Get my packages:
```{r}
library(ggplot2)
library(tidyr)
#install.packages("Reunran","foreach")
#library(Runuran)
library(foreach)
```


###Try Lewis code (I modified it)
```{r}
# load parameters & functions
source("C:/Users/provo/Documents/GitHub/popdy/parameters_pacifichake.r")
source("C:/Users/provo/Documents/GitHub/popdy/functions.r")

# note: in 'functions' this includes - 
# 1) assemble_leslie
# 2) is_mature

# length of simulation in years
time = 5100

# if want recruitment variability, define eta as a normal random variate for larval mortality, # producing lognormal survival deviations
eta = rnorm(n = time, mean = 0, sd = eta_sd) - ((eta_sd ^ 2) / 2)

# if you want M deterministic
M = rep(natural_mort, time)

# WRITE FUNCTION TO SET UP STRUCTURES FOR STORING DATA IN THE CORRECT FORMAT, GIVEN SCENARIO ARGUMENTS
# set up empty matrix for size at age at time
W = matrix(0, a_max, time)
# If deterministic growth, define W as same accross all times (FOR PERFORMANCE, IN FUTURE DON'T DEFINE FOR EACH TIME)
W[,] = W_inf * (1 - exp(-k * ((1:a_max) - t_0))) ^ b
# Else stochastic growth
# if updating W_inf within function, just insert deterministic size at age into first column of empty size at age over time matrix
#W[,1] = W_inf * (1 - exp(-k * ((1:a_max) - t_0))) ^ b
# else calculate random series of W_inf and use to preallocate series of size at age over time
#W_inf = urlnorm(time, meanlog = log(W_inf), sdlog = log(W_inf_sd), lb = 0, ub = W_inf_max - W_inf_min) + W_inf_min
# NEW GROWTH FORM NOT WORKING, USE LENGTH-BASED VERSION BELOW
#for(j in 2:time) {
#  W[1,j] =  W_inf[j-1] * (1 - exp(-k * (1 - t_0))) #^ b #W_0 + (W_inf[j-1] * exp(-b * k * (1-t_0)) * (1 - exp(-b * k)))
#  W[2:a_max,j] = W[1:(a_max-1),j-1] + W_inf[j-1] * exp(-k * ((1:(a_max-1)) - t_0)) * (1 - exp(-k)) 
#}


# set up empty matrix for maturity at age, given size, at time
mat = matrix(0, a_max, time)
# if deterministic growth, fill mat with constant values for maturity at age, based on deterministic growth
# plogis() gives the distribution function for the logistic distribution with parameters 'location' and 'scale'
# the output from plogis() gives the prediction at age and time for maturity. 
# right now it's constant across all ages, this seems wrong (?)
# location parameter = mean, W_mat=size at 50% mat
# scale parameter = variance, standard deviation * (pi^2/3)
# see R help pg: https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Logistic.html

mat[,] <- plogis(W[,1], location = W_mat, scale =  W_mat_sd * (sqrt(3)/pi))

# set up empty matrix for abundance at age at time
n = matrix(NA, a_max, time)

# set up empty vector to store total egg production each year
eggs = rep(NA, time)
# set up empty vector to store catch each year
catch = rep(NA, time)

# keep seed constant for repeatability initially
seed <- 42
set.seed(seed)

# package parms with initial values for some variables, to pass to top-level function
parms_plus = c(eta = eta, M = M, time = time, W = W, mat = mat, n = n, eggs = eggs,
               a_max = a_max, k = k, W_inf = W_inf, b = b, L_W_slope = L_W_slope, catch = catch,
               W_inf = W_inf, t_0 = t_0, W_mat = W_mat, W_mat_sd = W_mat_sd, W_rec = W_rec,
               W_rec_sd = W_rec_sd, W_esc = W_esc, W_esc_sd = W_esc_sd, f_W_slope = f_W_slope,
               f_W_int = f_W_int, alpha = alpha, beta = beta, W_inf_sd = W_inf_sd, 
               W_inf_min = W_inf_min, W_inf_max = W_inf_max, R_0 = R_0, natural_mort = natural_mort)

# write function to simulate, given parms and value for harvest rate, H
# noting that 0.05 to 0.4 is the historical range for Pacific hake
simulate <- function(parms_plus = parms_plus, H = 0) {
  
  # fill in abundance at age for time = 1 given that pop is fished at fishing mortality H at all ages and at the SAD
  n[,1] = R_0 * exp(- (natural_mort + H) * (0:(a_max-1)))
  
  for(j in 2:time){ # step through each time step, now that we have abundance at age in t=1
    
    # if deterministic growth, or stochastic growth with prespecified forcing,
    # assign weight at age as constant by commenting out assignment to W_a_array
    # Else stochastic with update within loop, assign size at age for time j
    #W <- stochastic_growth2(growth_function = '2parmVB', W_a_array = W, year = j, W_inf_sd = W_inf_sd, 
    #                        W_inf_max = Winf_max, W_inf_min = W_inf_min, b = b, k = k, t_0 = t_0, W_0 = W_0)
      
      # if deterministic growth, assign mat at age as constant by commenting out assignment to mat
      # else stochatic, assign mat at age for time j
      # assign maturity at age for time j, given size
      #for(i in 2:a_max) {
      #  mat[i,j] <- is_mature(mat_function = 'logistic', mat_a_array = mat, W_a_array = W, age = i, 
      #                        year = j, location = W_mat, scale = W_mat_sd * (sqrt(3)/pi))
      #}

    # store results from computing Leslie matrix for time j
    leslie_eggs_catch <- assemble_leslie(SR = 'Beverton_Holt', mat_a_array = mat, W_a_array = W, 
                                         year = j, n = n, H = H, W_rec = W_rec, W_rec_sd = W_rec_sd, 
                                         W_esc = W_esc, W_esc_sd = W_esc_sd, M = M,
                                         f_W_slope = f_W_slope, alpha = alpha, beta = beta, eta = eta)
    # extract Leslie matrix for time j
    A <- leslie_eggs_catch$A
    # store egg production for time j
    eggs[j] <- leslie_eggs_catch$E
    # store catch for time j
    catch[j] <- leslie_eggs_catch$C
  
    # perform population projection of one time step
    n[,j] <- A %*% n[,j-1]
  
    # make sure abundances at age that are less than 1 become zero
    for(i in 1:a_max) { if(n[i,j] < 1) n[i,j] = 0 }
    # if whole population goes below 1, stop the simulation
    if(sum(n[,j]) < 1) print("Population crashed!") 
  }
# Return list including numbers at age, egg production, and catch for the whole timeseries, 
# with first 100 time steps trimmed off (to make sure we are around the equilibrium)
# n = a matrix with a_max rows, and time-100 columns
# eggs = a vector, egg production each yr
# catch = a vector, catch each yr
return(list(n=n, eggs=eggs, catch=catch))
}


# simulate for a single harvest rate H, noting that H_crit is close to 3.1
#output = simulate(H=0)

output = foreach(H. = seq(0, 2.4, 0.4)) %do% simulate(H = H.)


```
