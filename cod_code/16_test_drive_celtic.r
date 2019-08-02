# **************** #
# Test Drive Celtic
# **************** #
# modify the Celtic Leslie matrix so that peak age is 5, 10
#Leslie3d <- Aarray[["Celtic"]]
#cel <- Leslie3d[,,which(1.0==kvals)[[1]]]
#celtop <- cel[1,]
#addages <- c(2,7,22)
# note: peak age of Celtic is 3, with added ages
# new peak ages are: 5(3+2), 10(3+7), 25(3+22)
#celtopplus <- c(rep(0,length=addages[3]),celtop) 

library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
ynar <-  c(0,0,0,0,0,0,0.05,0.25,0.7)
ywide <- c(0,0.01,0.01,0.03,0.1,0.1,0.15,0.2,0.4)
onar <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0.05,0.25,0.7)
owide <- c(0,0,0,0,0,0,0,0,0.01,0.01,0.03,0.1,0.1,0.15,0.2,0.4)
length(ynar)
length(ywide)
length(onar)
length(owide)
vecs <- list(ynar,ywide,onar,owide)
names(vecs) <- c("ynar","ywide","onar","owide")


yn <- matrix(0,nrow=length(ynar),ncol=length(ynar))
yn[1,] <- ynar
yw <- matrix(0,nrow=length(ywide),ncol=length(ywide))
yw[1,] <- ywide
on <- matrix(0,nrow=length(onar),ncol=length(onar))
on[1,] <- onar
ow <- matrix(0,nrow=length(owide),ncol=length(owide))
ow[1,] <- owide

for(i in 1:(length(ynar)-1)){yn[1+i,i] <- 1}
for(i in 1:(length(ywide)-1)){yw[1+i,i] <- 1}
for(i in 1:(length(onar)-1)){on[1+i,i] <- 1}
for(i in 1:(length(owide)-1)){ow[1+i,i] <- 1}

LesliesL <- list(yn,yw,on,ow)
names(LesliesL) <- c("young peak, small Stdev","young peak, wide Stdev","old peak, small Stdev","old peak, wide Stdev")

# set params for simulation:
timesteps = 2000 #need this now to create
rm_first_timesteps = 1000
betas = 1000
alphas <- selectedalphas
sig_r = 0.3
span.multiplier = 1 # adjusting the span in spec.prgm()

outputL <- as.list(rep(NA,length=length(LesliesL)))

for (l in 1:length(LesliesL)){ #step through oynw
  output.matrix <- array(NA,c(timesteps-2,4,length(selectedalphas)))
  
  for (a in 1:length(alphas)) { #for each oynw, step through alphas
    output = sim_model(A=LesliesL[[l]], timesteps=timesteps, 
                     alpha=selectedalphas[a], beta=betas, 
                     sig_r=sig_r, initial_eggs=betas)
    length(output$Nsize) <- length(output$N_t) #trim Nsize ts vector, -2 elements
    output.matrix[,,a] <- do.call(cbind,output) #fill in array for pop i
    colnames(output.matrix) <- names(output)
    print(a)
    }
  outputL[[l]] <- output.matrix
  print(l)
}
names(outputL) <- names(LesliesL)
rm(l,a,output.matrix)

# Format output ts for plotting simulations using output.3d.list
variable_type <- c("Nt","eggs","recruits","Nsize")
var.number <- 3 # recruits

# reformat data to work with ggplot
tsL <- as.list(rep(NA,length=length(outputL)))
names(tsL) <- names(outputL)
for (i in 1:length(outputL)){
  L <- outputL[[i]]
  aa <- as.data.frame(L[,var.number,])
  aa$year <- seq(from=1, to=length(aa[,1]))
  colnames(aa) <- c(selectedalphas,"year")
  aa1 <- aa %>% gather(alphaval,value,1:length(selectedalphas))
  aa1$variable <- rep(names(outputL)[i],length=length(aa1[,1]))
  tsL[[i]] <- aa1
}
rm(L,aa,i,aa1)
names(tsL) <- names(LesliesL)

span.multiplier= 1

# setting 'span' - a vector of odd integers to specify the smoothers
tmp <- ceiling(sqrt(length(1:(timesteps-rm_first_timesteps-1)))) #sq root of timeseries lgth, rounded
if (tmp %% 2 == 0) {m <- tmp+1} else {m <- tmp} #make it odd, if the square root is even
m = m * span.multiplier
  
# spectral analysis on RECRUIT, loop over tsL
spsaveL <- as.list(rep(NA,length=length(tsL)))
for(i in 1:length(tsL)){ #step through yonw
  aa1 <- tsL[[i]]
  spL <- list(rep(NA,length=length(selectedalphas))) #put spec in here
  
  for (b in 1:length(selectedalphas)){
    y = aa1[aa1$alphaval == selectedalphas[b],]$value[rm_first_timesteps:(timesteps-2)]
    yy = y-mean(y)
    sp = spec.pgram(yy,spans=c(m,m),taper=0.1,plot = FALSE)
    spL[[b]] = 2*sp$spec # save matrix of spec values for different FLEP, index by pop i
    print(b)
  }
  spdf <- bind_cols(spL)
  names(spdf) <- c('1.1','1.4','2','3.3','10')
  spdf$variable <- rep(names(tsL)[i],length=spdf[,1])
  spdf$freq <- sp$freq
  spsaveL[[i]] <- spdf
}
rm(i,b,spL,aa1,sp,yy)
spsave <- bind_rows(spsaveL)
spsavelong <- spsave %>% gather(alphavalue, value, 1:length(selectedalphas))
head(spsavelong)

# --- recruits spectra --- #
cs <- scales::seq_gradient_pal("#bdd7e7", "#08519c", "Lab")(seq(0,1,length.out=length(selectedalphas)))
pp <- as.list(rep(NA,length=length(LesliesL)))
str(spsavelong)
spsavelong$alphavalue <- as.numeric(as.character(spsavelong$alphavalue))
spsavelong$kval<- round((1/spsavelong$alphavalue),digits=2)
spsavelong$alphavalue <- factor(spsavelong$alphavalue)
spsavelong$kval <- factor(spsavelong$kval)

for(i in 1:length(LesliesL)){

  df <- spsavelong[spsavelong$variable == names(LesliesL)[i],]
  
  pp[[i]] <- ggplot(data=df, aes(x=freq,y=log(value),group=alphavalue)) + 
    geom_line(aes(color=kval)) + 
    theme_classic() + ylab("") + xlab("") + 
    scale_colour_manual(values=cs) +
    ggtitle(names(LesliesL)[i]) +
    theme(plot.title = element_text(size = 10)) 
}
tiff(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/manuscript/Test_drive_young-old_wide-narrow.tiff', units="in", width=6, height=5.5, res=300) 
do.call(grid.arrange,c(pp,ncol=2))
dev.off()
# -------------- AUC ------------------
freq <- spsavelong[1,1]
#threshold for high/low frequencies is 1/(2T) 
newpeak = 3 + addages[3]
AUC_less_L <- as.list(rep(NA,length=length(selectedalphas)))
AUCperlow_L <- as.list(rep(NA,length=length(selectedalphas)))
AUC_total_L <- as.list(rep(NA,length=length(selectedalphas)))

names(AUC_less_L) <- as.character(selectedalphas)
names(AUC_total_L) <- selectedalphas
names(AUCperlow_L) <- selectedalphas

for (j in 1:length(selectedalphas)){ #for each alpha (ie kval)...
  
  
  AUC_less[j] <- sum(freq*spsavelong[spsavelong$alphavalue == selectedalphas[j] 
                                         & spsavelong$freq <= (1/(25*2)),]$value)
    
    AUC_total[j] <- sum(freq*spsavelong[spsavelong$alphavalue == selectedalphas[j],]$value)
    AUCperlow[j] <- AUC_less[i]/AUC_total[i]
    
  }

rm(i,j)
AUC_less_df <- data.frame(do.call(cbind,AUC_less_L))
AUC_total_df <- data.frame(do.call(cbind,AUC_total_L))
AUCperlow_df <- data.frame(do.call(cbind,AUCperlow_L))

AUC_less_df$codNames <- codNames_ordered_by_peak
AUC_total_df$codNames <- codNames_ordered_by_peak
AUCperlow_df$codNames <- codNames_ordered_by_peak



