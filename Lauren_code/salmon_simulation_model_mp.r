#2a) New Parameters for No Hatchery, Fishing
a=60 
b=0.00017 
tf=2000
N0=c(100,0,0,0,0) 
s=0.28
e=0.1056 
l=0.1056
sx=c(s,s,(s*(1-e)),(s*(l)))
t<-1
ncls = length(N0) #Number of age classes

#2b) Simulation with random variability in recruitment for No Hatchery, Fishing 
AgeStructMatrix_F = function(sx,a,b,tf,N0) { 
  set.seed(1)
  sig_r=0.3
  Nt_F= matrix(0,tf,ncls+2) 
  Nt_F[1,] = c(0,0,N0) 
  for(t in 1:(tf-1)) {
    Pt= (e*Nt_F[t,5])+((1-l)*Nt_F[t,6])+Nt_F[t,7] 
    Nt_F[t+1,1] = Pt 
    Nt_F[t+1,2] = (a*Pt)/(1+(b*Pt)) 
    Nt_F[t+1,3] = (Nt_F[t+1,2])*(exp(sig_r*rnorm(1,mean=0, sd=1))) 
    Nt_F[t+1,4:(ncls+2)] = sx*Nt_F[t,3:(ncls+2-1)] 
  }
  return(Nt_F)
}

#sample run
Nt_F=AgeStructMatrix_F(sx,a,b,tf,N0) 


# output to be fed into wavelet analysis
#Spawning Stock only (age classes 3,4,5)


x_F_S=matrix(0,tf,1)
u<-1
for(u in 1:tf) {
  x_F_S[u,1]=Nt_F[u,1]
}
u<-u+1

Years=c(1901:2000)

x_F_S_wave=c(x_F_S[1901:2000])
Stock=rep('No_Hatchery_and_Fishing',100)
x_F_S_wave=cbind(Stock,Years,x_F_S_wave)
colnames(x_F_S_wave)<-c("Stock", "Years", "Spawners")
x_F_S_wave=as.data.frame(x_F_S_wave)
