# ----------------------------------------------------------
# Re-create spawning biomass distributions from LSB txt files

# prep for plotting
rm(list=ls())
codpop_name <- c("Northsea","Coastal Norway","West Baltic",
                  "Faroe","NE Arctic","Celtic",
                  "Iceland","Kattegat","West Scotland",
                  "N Gulf of St. L","Georges Bank","Gulf of Maine",
                  "3NO","3M","2J3KL",
                  "3Ps")
codpops <- c("Northsea","Coas","W_Baltic",
              "Faroe","NE_Arctic","Celtic",
              "Iceland","Kat","W_Scotland",
              "NGulf","GB","GM",
              "3NO","3M","2J3KL",
              "3Ps")

Age = 1:40  #max age = 40 yrs
F.halfmax = seq(0,3,by=0.01)
Fs <- paste("F",F.halfmax, sep="") 
cod.no <- 1:length(codpops)
SBdfs_list = as.list(rep(NA,length(codpops)))

# multiplot panel of spawning biomass:
pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLSB/SBplot_nofishing4.pdf', width=7, height=10)
par(mfrow=c(5,3))
for (i in seq_along(codpops)) { 
 # datax <- read.table(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_code/huiyu/LSB_varyfishing/'
 #                                 ,codpops1[i], '.txt', sep=''),header=T)
  
  datax <- read.table(file = paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLSB/k1/",
                                   codpops[i], '.txt', sep=''),header=T)
  datax = rbind(Fs,datax) #add F labels to df as first row
  colnames(datax) = datax[1,] #change the first row to column name
  datax = datax[-1,] #remove first row (bc now it's the header)
  dataxx = subset(datax, select=c("F0")) #,"F0.5","F1")) #only plot for these values of F
  SBdfs_list[[i]] = datax #store data frame with fishing rate as header in a list
  matplot(x=Age, y=dataxx, type="l",ylim=c(0,1),main=codpop_name[i],ylab="",xlab="")
  
}
par(mfrow=c(1,1))
dev.off()

