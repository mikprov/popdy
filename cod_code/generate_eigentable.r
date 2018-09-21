# Generate 'eigentable'
# This table has key pop information for running different analyses
# temp, lambda1, lambda2, dampratio, mode, cvs_mode, etc

# Plan
# 1. load codNames vector (must run the 0_load_cod_data script)
# 2. in a loop, read in Leslie matricies & calc eigenvalues
# 3. in a loop, read in LSB at age to calculate mode, CV, sd
# 4. manually create temp vector (temps are in Table 1, Wang et al. 2014)
# 5. read in max ages
# 6. combine all information into a single table, export table

# Plan Part 2
# 7. generate new matricies - relative number of spawners at age
# 8. calculate eigenvalues and store in eigentable

# *** Important ***
# before running this script, make sure to have Leslie 
# and LSB-at-age matricies stored locally. These matricies 
# are generated in  1_LSB_Lesliematricies_cod_v3.rmd 

rm(list=ls())

# ---
# 1. load codNames vector (must run the 0_load_cod_data script)
# ---
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/0_load_cod_data.r")
codNames # check, should be 17 pops
# load cod functions (need eigenvalue functions)
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/2_cod_functions.r")

# ---
# 2. in a loop, read in Leslie matricies & calc eigenvalues
# ---
# remove Irish because no maturity coefficients
codNames <- c("Northsea","Coas","W_Baltic",
              "Faroe","NE_Arctic","Celtic",
              "Iceland","Kat","W_Scotland",
              "NGulf","GB","GM",
              "cod3NO","cod3M","cod2J3KL",
              "cod3Ps")

#first = rep(NA, length(codNames))
#second= rep(NA, length(codNames))
#dampratio = rep(NA, length(codNames))



#for (i in 1:length(codNames)) {
  # for each pop i, load the Leslie matrix
#  A = read.table(file=paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/k1_basematrix_maxages/",
#                    codNames[i],".txt",sep=""))
#  first[i] <- extract_first_eigen_value(A) #store lambda1
#  second[i] <- extract_second_eigen_value(A) #store lambda2
#  dampratio[i] <- abs(second[i])/first[i]
#}
#table1 <- data.frame(cbind(codNames, first, second, dampratio))
#rm(i,A,first,second,dampratio)



# ---
# # 3. in a loop, read in LSB at age to calculate mode, CV, sd
# ---
mode_age = rep(NA,length(codNames))
sd_mode = rep(NA,length(codNames))
cvs_mode = rep(NA,length(codNames))

for (i in 1:length(codNames)) {
LSB = read.table(file=paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLSB/LSBmatrix_maxages/",
                            codNames[i],".txt",sep=""))
  p_spawn = LSB[,1] / sum(LSB[,1]) # LSB[,1] is LSB at age with no fishing
  # probability of spawning at age = LSB at age/total LSB
  Age = seq(1,length(LSB[,1]),by=1)
  p_table = data.frame(cbind(Age,p_spawn))
  #keep <- p_table[which(p_table$p_spawn > 0.01),] # remove probabilities less than 0.01
  keep<- p_table
  # using mode
  mode_age[i] = keep$Age[which.max(keep$p_spawn)] # what is the age with highest probability?
  sd_mode[i] = sqrt( sum(keep$p_spawn*(keep$Age-mode_age[i])^2) ) # stdev
  cvs_mode[i] = sd_mode[i]/mode_age[i] # coefficient of variation 
}
table2 <- data.frame(cbind(codNames,mode_age,sd_mode,cvs_mode))
rm(i,LSB,p_spawn,Age,p_table,keep,mode_age,sd_mode,cvs_mode)


# ---
# 4. manually create temp vector (temps are in Table 1, Wang et al. 2014)
# ---
temp = rep(NA, length(codNames))
for (i in 1:length(codNames)) {
  source(file = paste('C:/Users/provo/Documents/GitHub/popdy/cod_pops/',codNames[i], '.r', sep=''))
  temp[i] = TEMP
  rm(TEMP)
}
table3 <- data.frame(cbind(codNames,temp))
rm(i,temp)


# ---
# 5. read in max ages
# ---
max_ages <- rep(NA, length(codNames))
max_ages_occurances <- rep(NA, length(codNames))

for (i in 1:length(codNames)) {
  pop <- datalist[[i]]
  max_ages[i] <- max(pop$AGE)
  max_ages_occurances[i] <- length(pop[pop$AGE == max_ages[i],]$AGE) 
}
table4<- data.frame(cbind(codNames, max_ages))
table4$max_ages <- as.numeric(as.character(table4$max_ages))
rm(i,pop,max_ages,max_ages_occurances)


# ---
# 6. combine all information into a single table, export table
# ---
#e1 <- merge(table2,by="codNames")
e2 <- merge(table2,table3,by="codNames")
eigentable<- merge(e2,table4,by="codNames")
rm(e2)
write.csv(eigentable,file="C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLSB/eigentable.csv")


