# Plot selectivity and maturity vs temp
# By: Mikaela Provost

# Plan
# 1. load selectivity ogive parms from Wang et al
# 2. load maturity ogive parms from Wang et al
# 3. read in eigenvalue table, it has the temp data
# 4. plotting

# load packages
library(ggplot2)

# 1. load maturity ogive parms from Wang et al
popname <- c('Celtic','Irish','Northsea','W_Scotland','W_Baltic','Kat',
       'E_Baltic','Coas','Faroes','Iceland','NE_Arctic','GB',
       'GM','3NO','SGulf','NGulf','3M','2J3KL','3Ps','4X')
mat50 <- c(2.3,1.97,3.73,2.12,2.96,2.28,3.17,
           5.29,2.99,6.76,8.01,2.01,2.44,6.01,
           4.36,5.28,3.59,5.97,6.07,0)
A <- data.frame(popname,mat50)
A

# 2. load selectivity ogive parms from Wang et al
sel50 <- c(1.63,1.41,1.79,3.56,1.81,2.52,3.09,
           6.42,4.7,8,6.15,7.76,3.82,7.31,8.04,
           7.38,4.31,7.47,7.19,0)
A$sel50 <- sel50

ggplot(eigentable, aes(x=cvs_mode,y=dampratio)) +
  geom_point(aes(color=temp), size=3) +  
  geom_text(aes(label=codNames),hjust=-0.1,vjust=-0.2,check_overlap = F) +
  scale_color_gradientn(colors=rev(rainbow(n=10,start=0,end=0.7)))

# 3. read in eigenvalue table, it has the temp data
eigentable <- read.csv('C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLSB/eigentable5.csv',header=T)

