# Plot F vs FLEP (same figure as Fig 3(b) in Botsford et al. 2014)

# Plan
# 1. 

library(gridExtra)
library(reshape2)
# ---
# 1) Define fishing levels 
# ---
F.halfmax = seq(0,100,by=0.01)

# ---
# 2) read in cod data
# ---
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/0_load_cod_data.r")
# this loads 'datalist' which has the data for each population
# combined into one list

# ---
# 3. load functions
# ---
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/2_cod_functions.r")
# assemble_Leslie()
# extract_first_eigen_value()
# calc_LSB_at_age_by_F()

# ---
# 4) Generate LSB for each pop
# ---
tknot = 0 #used in vonB equation 
#k = 1 #this sets no fishing as influenced on Leslie matrix


# loop over pop data in datalist to generate LSB matrices
LSBlist = as.list(rep(NA,length=length(codNames))) # store matrix of LSB at age vs F
names(LSBlist) = codNames # assign names to list
FLEP = matrix(NA,nrow=length(F.halfmax),ncol=length(codNames))
colnames(FLEP) = codNames
for (i in 1:length(codNames)) { # step through each pop in datalist
  # load parms for cod pop i
  source(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_pops/',names(LSBlist)[i], '.r', sep=''))
  # this should load parms: L_inf, K (vonB), TEMP, maxage
  LSBlist[[i]] = calculate_LSB_at_age_by_F(maxage=maxage, 
                                           L_inf=L_inf, 
                                           K=K, 
                                           TEMP=TEMP, 
                                           F.halfmax=F.halfmax,
                                           B0=B0,B1=B1)
  # here I calculate FLEP: LSB at all F levels / LSB at F=0
  FLEP[,i] <- colSums(LSBlist[[i]]) / sum(LSBlist[[i]][,1])
  
}

matplot(x=F.halfmax,y=FLEP,type="l")

# plot F vs FLEP
# store this df: the Fs associated with different FLEP values for each pop
FLEP.F <- as.data.frame(cbind(F.halfmax,FLEP))
write.csv(FLEP.F,file='C:/Users/provo/Documents/GitHub/popdy/cod_code/FLEP.F.csv')
# 
FLEPlong <- melt(FLEP.F,id="F.halfmax")
pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/F_vs_FLEP.pdf', width=7, height=7)
ggplot(data=FLEPlong, aes(x=F.halfmax,y=log(value), 
                          color=variable, linetype=variable, shape=variable)) + 
  geom_line() +
 # scale_linetype_manual(values=c(rep("solid",8),rep(dashed)))
  theme_classic() +
  xlab("F") +
  ylab("ln(FLEP)")

ggplot(data=FLEPlong, aes(x=F.halfmax,y=value, color=variable, linetype=variable, shape=variable)) + 
  geom_line() +
  # scale_linetype_manual(values=c(rep("solid",8),rep(dashed)))
  theme_classic() +
  xlab("F") +
  ylab("FLEP")
dev.off()
# ----
# plot egg production at age for a couple F values
# ----
# the FLEP values chosen correspond to Exploitation Index values (EI)
# FLEP = 1, 0.64, 0.46, 0.28, 0.145
# what are the F values that correspond to these FLEP values?

# second, match up F values
FLEP.F <- as.data.frame(cbind(F.halfmax,FLEP))
# third, find the nearest F values that correspond to the desired FLEP value
target_Fs <- rbind(

  as.data.frame(melt(FLEP.F,id="F.halfmax") %>% 
                  group_by(variable) %>%
                  arrange(abs(value-1)) %>% #desired FLEP=1
                  slice(1) %>%
                  #mutate(value=round(value,2)) %>%
                  mutate(FLEPlevel=rep(1))),
  
  as.data.frame(melt(FLEP.F,id="F.halfmax") %>% 
                  group_by(variable) %>%
                  arrange(abs(value-0.64)) %>% #desired FLEP=0.64
                  slice(1) %>%
                  #mutate(value=round(value,2)) %>%
                  mutate(FLEPlevel=rep(0.64))),

  as.data.frame(melt(FLEP.F,id="F.halfmax") %>% 
                group_by(variable) %>%
                arrange(abs(value-0.46)) %>% #desired FLEP=0.46
                slice(1) %>%
                #mutate(value=round(value,2)) %>%
                mutate(FLEPlevel=rep(0.46))),
  
  as.data.frame(melt(FLEP.F,id="F.halfmax") %>% 
                  group_by(variable) %>%
                  arrange(abs(value-0.28)) %>% #desired FLEP=0.28
                  slice(1) %>%
                  #mutate(value=round(value,2)) %>%
                  mutate(FLEPlevel=rep(0.28))),
  
  as.data.frame(melt(FLEP.F,id="F.halfmax") %>% 
                  group_by(variable) %>%
                  arrange(abs(value-0.145)) %>% #desired FLEP=0.145
                  slice(1) %>%
                  #mutate(value=round(value,2)) %>%
                  mutate(FLEPlevel=rep(0.145)))

)
plot(x=target_Fs$FLEPlevel, y=target_Fs$value)
target_Fs$diffs_in_F = target_Fs$FLEPlevel-target_Fs$value


# for pop i in LSBlist...
p <-list()

for (i in 1:length(LSBlist)) {
  f <- target_Fs[target_Fs$variable == names(LSBlist)[i],]
  f$F.halfmax <- paste("F",f$F.halfmax,sep="_") #converting F values to character, makes easier to grab columns later on
  lsb <- as.data.frame(LSBlist[[i]])
  colnames(lsb) <- paste("F",F.halfmax,sep="_")
  # get data ready to plot 
  forplot <- melt(lsb %>% 
                    select(f$F.halfmax) %>% 
                    mutate(age=seq(from=1,to=length(lsb[,1]),by=1)),
                  id="age")
  # plot egg production at age for different Fs
  p[[i]] <- ggplot(forplot,aes(x=age,y=value,color=variable)) +
    ylab("egg production") +
    xlab("age") +
    geom_line() +
    scale_color_brewer(palette = "Reds") +
    ggtitle(paste(names(LSBlist)[i])) +
    theme_classic()
} 

pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/egg_production_for_diff_Fs.pdf', width=7, height=10)
do.call(grid.arrange,c(p,ncol=2))
dev.off()


