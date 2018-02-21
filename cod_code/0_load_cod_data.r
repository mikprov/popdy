# Load cod data and break up into populations

cod = read.table("C:/Users/provo/Documents/GitHub/popdy/cod_code/huiyu/cod_all_2013.txt",header=T,na.strings="NA")
setwd("C:/Users/provo/Documents/GitHub/popdy/cod_code/huiyu/LSB_varyfishing")

NorthseaD = subset(cod, AREA=="NORTH_SEA")
CoasD = subset(cod, AREA=="COAS")
E_BalticD = subset(cod, AREA=="E_BALTIC")
W_BalticD = subset(cod, AREA=="W_BALTIC")
FaroeD = subset(cod, AREA=="FAROE")
NE_ArcticD = subset(cod, AREA=="NE_ARCTIC")
CelticD = subset(cod, AREA=="CELTIC_SEA")
IcelandD = subset(cod, AREA=="ICELAND")
IrishD = subset(cod, AREA=="IRISH_SEA")
KatD = subset(cod, AREA=="KATTEGAT")
W_ScotlandD = subset(cod, AREA=="W_SCOTLAND")

NGulfD = subset(cod, AREA =="N_GulfSL")
SGulfD = subset(cod, AREA =="S_GulgSL")
GBD = subset(cod, AREA =="GB")
GMD = subset(cod, AREA =="GM")
cod3MD = subset(cod,AREA=="cod3M")
cod3NOD = subset(cod,AREA=="cod3NO")
cod3PsD = subset(cod,AREA=="cod3Ps")
cod4XD = subset(cod,AREA=="cod4X")
cod2J3KLD = subset(cod,AREA=="cod2J3KL")


datalist <- list(NorthseaD,CoasD,W_BalticD,
                 FaroeD,NE_ArcticD,CelticD,
                 IcelandD,KatD,W_ScotlandD,
                 NGulfD,GBD,GMD,
                 cod3NOD,cod3MD,cod2J3KLD,
                 cod3PsD)

codNames <- c("Northsea","Coas","W_Baltic",
              "Faroe","NE_Arctic","Celtic",
              "Iceland","Kat","W_Scotland",
              "NGulf","GB","GM",
              "3NO","3M","2J3KL",
              "3Ps")

names(datalist) <- c("Northsea","Coas","W_Baltic",
                     "Faroe","NE_Arctic","Celtic",
                     "Iceland","Kat","W_Scotland",
                     "NGulf","GB","GM",
                     "3NO","3M","2J3KL",
                     "3Ps")

