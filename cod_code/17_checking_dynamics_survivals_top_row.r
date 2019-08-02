# Checking the dynamics when survival is moved from the
# sub-diagonal to the top row. 

# load parms for cod pop i: L_inf, K (for vonB), TEMP, maxage,B0,B1 (matur)
source(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_pops/',codNames[1], '.r', sep=''))
  
  
Leslieout = assemble_Leslie(maxage=maxage, K=K, L_inf=L_inf, TEMP=TEMP,
                                F.halfmax=0, B0=B0, B1=B1, tknot=0)
# move survivals along subdiagonal to fecundities
Leslieout$A[1,] <- Leslieout$A[1,]*Leslieout$A[2,1]
#I can do it this way bc survival is constant w/age (ie no fishing)

# set suvivals on subdiagonal = 1 (only works if no fishing)
Leslieout$A[Leslieout$A == Leslieout$A[2,1]] <- 1

# transform fecundity-at-age to probability density curve & multiply by k (slope)
Leslieout$A[1,] <- (Leslieout$A[1,]/sum(Leslieout$A[1,]))
    
e1 = extract_first_eigen_value(Leslieout$A)
e2 = extract_second_eigen_value(Leslieout$A)
e12 = e2 / e1

LeslieTEST <- Leslieout$A
NEAR_TEST <- Leslieout$NEAR
rm(L_inf,K,MG,Mp,theta0,theta1,TEMP,A50,S50,maxage,Leslieout,name)

# set params for simulation:
timesteps = 1000 #need this now to create
rm_first_timesteps = 200
betas = 1000
alphas <- 2 #these correspond to 
sig_r = 0.3

output_surv_top = sim_model(A=LeslieTEST, timesteps=timesteps, 
                       alpha=alphas, beta=betas, 
                       sig_r=sig_r, initial_eggs=betas)
plot(output_surv_subdiag$recruits[rm_first_timesteps:(timesteps-2)],type="l")
plot(output_surv_top$recruits[rm_first_timesteps:(timesteps-2)],type="l")


cv_surv_subdiag <- sd(output_surv_subdiag$recruits[rm_first_timesteps:(timesteps-2)])/mean(output_surv_subdiag$recruits[rm_first_timesteps:(timesteps-2)])

cv_surv_top <- sd(output_surv_top$recruits[rm_first_timesteps:(timesteps-2)])/mean(output_surv_top$recruits[rm_first_timesteps:(timesteps-2)])
