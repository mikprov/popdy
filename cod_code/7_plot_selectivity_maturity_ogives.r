# Plot selectivity and maturity vs temp
# By: Mikaela Provost

# Plan
# 1. load selectivity ogive parms from Wang et al
# 2. load maturity ogive parms from Wang et al
# 3. plotting

# load packages
library(ggplot2)
library(gridExtra)

# 1. load maturity ogive parms from Wang et al
popname <- c('Celtic','Irish','Northsea','W_Scotland','W_Baltic','Kat',
       'E_Baltic','Coas','Faroes','Iceland','NE_Arctic','GB',
       'GM','3NO','SGulf','NGulf','3M','2J3KL','3Ps','4X')
temps <- c(11.9,10.56,8.4,
           9.57,7,6.5,
           5,7.13,7.4,
           5.8,4,8,
           8,1.75,1.75,
           1,3.5,0,
           2.5,6.75)

mat50 <- c(2.3,1.97,3.73,2.12,2.96,2.28,3.17,
           5.29,2.99,6.76,8.01,2.01,2.44,6.01,
           4.36,5.28,3.59,5.97,6.07,NA)
A <- data.frame(popname,temps,mat50)


# 2. load selectivity ogive parms from Wang et al
sel50 <- c(1.63,1.41,1.79,3.56,1.81,2.52,3.09,
           6.42,4.7,8,6.15,7.76,3.82,7.31,8.04,
           7.38,4.31,7.47,7.19,NA)
A$sel50 <- sel50
A

# 3. plotting - make a 3 panel plot

pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/50percent_selectivity_maturity_plots.pdf', width=10, height=8)
a <- ggplot(A, aes(x=temps,y=mat50)) +
  geom_point(aes(color=temps), size=3) +  
  geom_text(aes(label=popname),hjust=-0.1,vjust=-0.2,check_overlap = F, size=2) +
  scale_color_gradientn(colors=rev(rainbow(n=10,start=0,end=0.7))) +
  ylim(0,10) +
  xlim(0,12) +
  ylab("age at 50% maturity") +
  xlab("temp") +
  ggtitle("(a)") +
  theme(plot.title = element_text(hjust=0))
 
b <- ggplot(A, aes(x=temps,y=sel50)) +
  geom_point(aes(color=temps), size=3) +  
  geom_text(aes(label=popname),hjust=-0.1,vjust=-0.2,check_overlap = F, size=2) +
  scale_color_gradientn(colors=rev(rainbow(n=10,start=0,end=0.7))) +
  ylim(0,10) +
  xlim(0,13) +
  ylab("age at 50% selectivity") +
  xlab("temp") +
  ggtitle("(b)") +
  theme(plot.title = element_text(hjust=0))

c <- ggplot(A, aes(x=temps,y=mat50/sel50)) +
  geom_point(aes(color=temps), size=3) +  
  geom_text(aes(label=popname),hjust=-0.1,vjust=-0.2,check_overlap = F, size=2) +
  scale_color_gradientn(colors=rev(rainbow(n=10,start=0,end=0.7))) +
  geom_hline(yintercept = 1, linetype="dashed") +
  ylim(0,3) +
  xlim(0,13) +
  ylab("mat50 / sel50") +
  xlab("temp") +
  ggtitle("(c)") +
  theme(plot.title = element_text(hjust=0))

grid.arrange(a, b, c, ncol=2, top="Age at 50% mature (mat50) and age at 50% selectivity (sel50) with temperature")
dev.off()
  
A$mat50/
