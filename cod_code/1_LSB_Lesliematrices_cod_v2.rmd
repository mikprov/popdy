#' ---
#' title: "Create Leslie matrices & LSB matrices from cod data"
#' author: "Mikaela Provost"
#' date: "December 6th, 2017"
#' ---

# ---
# plan:
# 1. read cod data, break in separate populations
# 2. load functions


rm(list=ls())

# ---
# 1. load cod data, break into separate populations
# ---
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/0_load_cod_data.r")

# ---
# 2. load functions
# ---
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/2_cod_functions.r")
  # assemble_Leslie()
  # extract_first_eigen_value()
  # calc_LSB_at_age_by_F()


# ---
# 3. execute functions for each population of cod, export Leslie & LSB matrices
# ---

# prep for looping
Age = 1:40  #max age = 40 yrs
to = 0
F.halfmax = seq(0,3,by=0.01)
k = 1


# loop over pop data in datalist to generate Leslie matrices
Alist = as.list(rep(NA,length(datalist))) # store Leslie matrix
names(Alist) = codNames
for (i in 1:length(datalist)) {
  Alist[[i]]=assemble_Leslie(data=datalist[[i]], codPopname = codNames[i], littlek=k)
  write.table(Alist[[i]],
              file=paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/k1/",
                         codNames[i],".txt",sep=""))
}

# loop over pop data in datalist to generate LSB matrices
LSBlist = as.list(rep(NA,length(datalist))) # store matrix of LSB vs F
names(LSBlist) = codNames
for (i in 1:length(datalist)) { # step through each pop in datalist
  LSBlist[[i]] = calculate_LSB_at_age_by_F(data=datalist[[i]], codPopname=codNames[i], littlek=k)
  # store the output from calculate_LSB.. function in list. output is LSB by age and F
  write.table(LSBlist[[i]],
              file=paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLSB/k1/",
                         codNames[i],".txt",sep=""))
}





