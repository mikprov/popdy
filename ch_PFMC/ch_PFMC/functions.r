# PFMC functions
# similar to cod functions

#a) functions for convert age to length (vonB & Schnute) 
vonbertgrowth <- function(maxage,t0,Linf,Kvonb,ages){
  ages <- seq(from=1,to=maxage,by=1)
  L_a = Linf * (1-exp(-Kvonb*(ages-t0)))
  return(L_a)
  rm(ages)
}

schnutegrowth <- function(maxage,t1,t2,L1,L2,Ksch,ages){
  ages <- seq(from=1,to=maxage,by=1)
  L_a = L1 + (L2-L1) * (1-exp(-Ksch*(ages-t1)))/(1-exp(-Ksch*(t2-t1)))
  return(L_a)
  rm(ages)
}


# (b) function for calculate prop mature at age
calc_mat_at_age <- function(maxage,slope,inflection){
  ages <- seq(from=1,to=maxage,by=1)
  propmat_a <- 1/(1+exp(slope*(ages-inflection)))
  return(propmat_a)
  rm(ages)
}


# (c) function for selectivity at age 
vul_a

# (d) function for weight at age
calc_wt_at_age <- function(L_a,cmkga,cmkgb){
  B_a <- cmkga*(L_a)^cmkgb
  return(B_a)
}

# (e) functions for fecundity at age (2 versions)
fecFUN1 <- function(B_a,fecintercept,fecslope){
  eggs_a = B_a*(fecintercept + fecslope*B_a) 
  return(eggs_a)}
fecFUN2 <- function(B_a,fecintercept,fecslope){
  eggs_a = fecintercept + (fecslope*B_a) 
  return(eggs_a)}

# (f) assemble Leslie
assemble_Leslie <- function(maxage,L_a,B_a,propmat_a,eggs_a,vul_a,M,Fval){
  
  # create vector of ages
  Age=1:maxage 
  
  # assemble LTABLE: 
  LTABLE = data.frame(cbind(Age,L_a,B_a,eggs_a,propmat_a,vul_a))  
  LTABLE$M = rep(M,length=length(LTABLE[,1])) #fill in natural mortality
  LTABLE$Fval = rep(Fval,length=length(LTABLE[,1]))#fill in F mortality
  
  # empty Leslie matrix 
  A = matrix(0,length(Age),length(Age)) 
  
  # for each age, get survival at age from F&M
  for(j in 1:length(Age)){ # step through ages
    LTABLE$FISH[j] = LTABLE$vul_a[j]*Fval	# FISH: F rate = vulnerability-at-age * F  
    LTABLE$SURV[j] = exp(-(LTABLE$FISH[j]+LTABLE$M[j])) #SURV:fraction surviving at each age
    
    # set up column for survivorship (amount or fraction present at age)
    LTABLE$Survship = 0 
    LTABLE$Survship[1] = 1}
  
  # for each row in LTABLE (that is, for each age)
  for(k in 1:(nrow(LTABLE)-1)){ # step through ages  
    LTABLE$Survship[k+1] = LTABLE$Survship[k]*LTABLE$SURV[k]} #fraction surv from 1 to a
  
  # LEP at age
  LTABLE$LEP_a <- LTABLE$eggs_a*NEAR$Survship
  
  # insert fecundity (maturity * eggs) in top row
  A[1,] = LTABLE$propmat_a*LTABLE$eggs_a 
  
  # make LEP equal for all spp by adjusting egg (LEP=1.1)
  #A[1,] = A[1,]/(sum(LTABLE$LEP_a)*0.9) 
  
  # insert survival on subdiagonal (has both M&F)
  for(u in 2:length(Age)-1){ 
    A[u+1,u]=LTABLE$SURV[u] }
  
  # -- export Leslie matrix (A) & LTABLE df for one F value
  return(list(A=A,LTABLE=LTABLE))
  
}