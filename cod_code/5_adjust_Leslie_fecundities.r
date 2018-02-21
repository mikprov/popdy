# ---
# Adjusting the Leslie matrix fecundities so that lambda=1
# Plan: for each population,
# a) test different multipliers to adjust fecundities
# b) you want lambda to equal 1
# c) keep track of which multipliers you use, save one that works
# d) remember to change i in codNames[i] for each test

rm(list=ls())

# ---
# load cod data, break into separate populations
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/0_load_cod_data.r")


# ---
# Load functions:
#  extract_first_eigen_value()
#  extract_second_eigen_value()
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/2_cod_functions.r")


# ---
# Store the multiplier values here
multiplier_values <- rep(NA,length(codNames))


# ---
# Population 1: codNames[1] = Northsea
# Test history:
Multiplier = 1/1.0705509
Multiplier = 1/1.0805509
Multiplier = 1/1.45
A = read.table(paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/k1/",
                     codNames[1],".txt",sep=""))
f = extract_first_eigen_value(Lesliematrix=A)
s = extract_second_eigen_value(Lesliematrix=A)
# This is Loo trying to adjust lambda to equal 1 for each Leslie matrix
AA = A[1,]*Multiplier # multiply fecundities by 1/lambda
A[1,] <- AA # replace the first row, will new fecundities
fadj = extract_first_eigen_value(Lesliematrix=A) # recalc lambda 1
sadj = extract_second_eigen_value(Lesliematrix=A)
c(f,fadj,sadj)

multiplier_values[1] <- 1/1.45


# ---
# Population: codNames[2] = Coas
# Test history:
Multiplier = 1/0.885766581412765
Multiplier = 1/0.7
Multiplier = 1/0.6
Multiplier = 1/0.4
Multiplier = 1/0.39
A = read.table(paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/k1/",
                     codNames[2],".txt",sep=""))
f = extract_first_eigen_value(Lesliematrix=A)
s = extract_second_eigen_value(Lesliematrix=A)
# This is Loo trying to adjust lambda to equal 1 for each Leslie matrix
AA = A[1,]*Multiplier # multiply fecundities by 1/lambda
A[1,] <- AA # replace the first row, will new fecundities
fadj = extract_first_eigen_value(Lesliematrix=A) # recalc lambda 1
sadj = extract_second_eigen_value(Lesliematrix=A)
c(f,fadj,sadj)

multiplier_values[2] <- 1/0.39


# ---
# Population: codNames[3] = W_Baltic
# Test history:
Multiplier = 1/0.944835228663227
Multiplier = 1/0.71

A = read.table(paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/k1/",
                     codNames[3],".txt",sep=""))
f = extract_first_eigen_value(Lesliematrix=A)
s = extract_second_eigen_value(Lesliematrix=A)
# This is Loo trying to adjust lambda to equal 1 for each Leslie matrix
AA = A[1,]*Multiplier # multiply fecundities by 1/lambda
A[1,] <- AA # replace the first row, will new fecundities
fadj = extract_first_eigen_value(Lesliematrix=A) # recalc lambda 1
sadj = extract_second_eigen_value(Lesliematrix=A)
c(f,fadj,sadj)

multiplier_values[3] <- 1/0.71


# ---
# Population: codNames[4] = Faroe
# Test history:
Multiplier = 1/1.07379285729038
Multiplier = 1/1.1
Multiplier = 1/1.3
Multiplier = 1/1.46

A = read.table(paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/k1/",
                     codNames[4],".txt",sep=""))
f = extract_first_eigen_value(Lesliematrix=A)
s = extract_second_eigen_value(Lesliematrix=A)
# This is Loo trying to adjust lambda to equal 1 for each Leslie matrix
AA = A[1,]*Multiplier # multiply fecundities by 1/lambda
A[1,] <- AA # replace the first row, will new fecundities
fadj = extract_first_eigen_value(Lesliematrix=A) # recalc lambda 1
sadj = extract_second_eigen_value(Lesliematrix=A)
c(f,fadj,sadj)

multiplier_values[4] <- 1/1.46


# ---
# Population: codNames[5] = NE_Arctic
# Test history:
Multiplier = 1/1.00446270089909
Multiplier = 1/1.09999
Multiplier = 1/1.1
Multiplier = 1/1.099999999

A = read.table(paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/k1/",
                     codNames[5],".txt",sep=""))
f = extract_first_eigen_value(Lesliematrix=A)
s = extract_second_eigen_value(Lesliematrix=A)
# This is Loo trying to adjust lambda to equal 1 for each Leslie matrix
AA = A[1,]*Multiplier # multiply fecundities by 1/lambda
A[1,] <- AA # replace the first row, will new fecundities
fadj = extract_first_eigen_value(Lesliematrix=A) # recalc lambda 1
sadj = extract_second_eigen_value(Lesliematrix=A)
c(f,fadj,sadj)

multiplier_values[5] <- 1/1.09999999


# ---
# Population: codNames[6] = Celtic 
# Test history:
Multiplier = 1/1.08665642398179
Multiplier = 1/1.08
Multiplier = 1/1.09
Multiplier = 1/1.1
Multiplier = 1/1.2
Multiplier = 1/1.4
Multiplier = 1/1.38


A = read.table(paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/k1/",
                     codNames[6],".txt",sep=""))
f = extract_first_eigen_value(Lesliematrix=A)
s = extract_second_eigen_value(Lesliematrix=A)
# This is Loo trying to adjust lambda to equal 1 for each Leslie matrix
AA = A[1,]*Multiplier # multiply fecundities by 1/lambda
A[1,] <- AA # replace the first row, will new fecundities
fadj = extract_first_eigen_value(Lesliematrix=A) # recalc lambda 1
sadj = extract_second_eigen_value(Lesliematrix=A)
c(f,fadj,sadj)

multiplier_values[6] <- 1/1.38


# ---
# Population: codNames[7] = Iceland 
# Test history:
Multiplier = 1/0.984774076389807
Multiplier = 1/1.1
Multiplier = 1/0.9

A = read.table(paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/k1/",
                     codNames[7],".txt",sep=""))
f = extract_first_eigen_value(Lesliematrix=A)
s = extract_second_eigen_value(Lesliematrix=A)
# This is Loo trying to adjust lambda to equal 1 for each Leslie matrix
AA = A[1,]*Multiplier # multiply fecundities by 1/lambda
A[1,] <- AA # replace the first row, will new fecundities
fadj = extract_first_eigen_value(Lesliematrix=A) # recalc lambda 1
sadj = extract_second_eigen_value(Lesliematrix=A)
c(f,fadj,sadj)

multiplier_values[7] <- 1/0.9


# ---
# Population: codNames[8] = Kat
# Test history:
Multiplier = 1/0.955227297653795
Multiplier = 1/0.8
Multiplier = 1/0.76

A = read.table(paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/k1/",
                     codNames[8],".txt",sep=""))
f = extract_first_eigen_value(Lesliematrix=A)
s = extract_second_eigen_value(Lesliematrix=A)
# This is Loo trying to adjust lambda to equal 1 for each Leslie matrix
AA = A[1,]*Multiplier # multiply fecundities by 1/lambda
A[1,] <- AA # replace the first row, will new fecundities
fadj = extract_first_eigen_value(Lesliematrix=A) # recalc lambda 1
sadj = extract_second_eigen_value(Lesliematrix=A)
c(f,fadj,sadj)

multiplier_values[8] <- 1/0.76


# ---
# Population: codNames[9] = W_Scotland
# Test history:
Multiplier = 1/0.8916947188149
Multiplier = 1/0.7
Multiplier = 1/0.6
Multiplier = 1/0.56
Multiplier = 1/0.57

A = read.table(paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/k1/",
                     codNames[9],".txt",sep=""))
f = extract_first_eigen_value(Lesliematrix=A)
s = extract_second_eigen_value(Lesliematrix=A)
# This is Loo trying to adjust lambda to equal 1 for each Leslie matrix
AA = A[1,]*Multiplier # multiply fecundities by 1/lambda
A[1,] <- AA # replace the first row, will new fecundities
fadj = extract_first_eigen_value(Lesliematrix=A) # recalc lambda 1
sadj = extract_second_eigen_value(Lesliematrix=A)
c(f,fadj,sadj)

multiplier_values[9] <- 1/0.57


# ---
# Population: codNames[10] = NGulf
# Test history:
Multiplier = 1/1.2494941761164
Multiplier = 1/1.6
Multiplier = 1/2
Multiplier = 1/5
Multiplier = 1/8

A = read.table(paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/k1/",
                     codNames[10],".txt",sep=""))
f = extract_first_eigen_value(Lesliematrix=A)
s = extract_second_eigen_value(Lesliematrix=A)
# This is Loo trying to adjust lambda to equal 1 for each Leslie matrix
AA = A[1,]*Multiplier # multiply fecundities by 1/lambda
A[1,] <- AA # replace the first row, will new fecundities
fadj = extract_first_eigen_value(Lesliematrix=A) # recalc lambda 1
sadj = extract_second_eigen_value(Lesliematrix=A)
c(f,fadj,sadj)

multiplier_values[10] <- 1/8


# ---
# Population: codNames[11] = GB
# Test history:
Multiplier = 1/1.1198214438462
Multiplier = 1/1.71

A = read.table(paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/k1/",
                     codNames[11],".txt",sep=""))
f = extract_first_eigen_value(Lesliematrix=A)
s = extract_second_eigen_value(Lesliematrix=A)
# This is Loo trying to adjust lambda to equal 1 for each Leslie matrix
AA = A[1,]*Multiplier # multiply fecundities by 1/lambda
A[1,] <- AA # replace the first row, will new fecundities
fadj = extract_first_eigen_value(Lesliematrix=A) # recalc lambda 1
sadj = extract_second_eigen_value(Lesliematrix=A)
c(f,fadj,sadj)

multiplier_values[11] <- 1/1.71


# ---
# Population: codNames[12] = GM
# Test history:
Multiplier = 1/1.09752032619776
Multiplier = 1/1.65
Multiplier = 1/1.62

A = read.table(paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/k1/",
                     codNames[12],".txt",sep=""))
f = extract_first_eigen_value(Lesliematrix=A)
s = extract_second_eigen_value(Lesliematrix=A)
# This is Loo trying to adjust lambda to equal 1 for each Leslie matrix
AA = A[1,]*Multiplier # multiply fecundities by 1/lambda
A[1,] <- AA # replace the first row, will new fecundities
fadj = extract_first_eigen_value(Lesliematrix=A) # recalc lambda 1
sadj = extract_second_eigen_value(Lesliematrix=A)
c(f,fadj,sadj)

multiplier_values[12] <- 1/1.62


# ---
# Population: codNames[13] = 3NO
# Test history:
Multiplier = 1/1.17731095601338
Multiplier = 1/1.3
Multiplier = 1/2
Multiplier = 1/5
Multiplier = 1/5.1
Multiplier = 1/5.2

A = read.table(paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/k1/",
                     codNames[13],".txt",sep=""))
f = extract_first_eigen_value(Lesliematrix=A)
s = extract_second_eigen_value(Lesliematrix=A)
# This is Loo trying to adjust lambda to equal 1 for each Leslie matrix
AA = A[1,]*Multiplier # multiply fecundities by 1/lambda
A[1,] <- AA # replace the first row, will new fecundities
fadj = extract_first_eigen_value(Lesliematrix=A) # recalc lambda 1
sadj = extract_second_eigen_value(Lesliematrix=A)
c(f,fadj,sadj)

multiplier_values[13] <- 1/5.2


# ---
# Population: codNames[14] = 3M
# Test history:
Multiplier = 1/1.0555877388104
Multiplier = 1/1.1
Multiplier = 1/1.5
Multiplier = 1/1.6

A = read.table(paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/k1/",
                     codNames[14],".txt",sep=""))
f = extract_first_eigen_value(Lesliematrix=A)
s = extract_second_eigen_value(Lesliematrix=A)
# This is Loo trying to adjust lambda to equal 1 for each Leslie matrix
AA = A[1,]*Multiplier # multiply fecundities by 1/lambda
A[1,] <- AA # replace the first row, will new fecundities
fadj = extract_first_eigen_value(Lesliematrix=A) # recalc lambda 1
sadj = extract_second_eigen_value(Lesliematrix=A)
c(f,fadj,sadj)

multiplier_values[14] <- 1/1.6


# ---
# Population: codNames[15] = 
# Test history:
Multiplier = 1/1.22846310973873
Multiplier = 1/1.5
Multiplier = 1/5
Multiplier = 1/6
Multiplier = 1/7
Multiplier = 1/8
Multiplier = 1/8.1

A = read.table(paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/k1/",
                     codNames[15],".txt",sep=""))
f = extract_first_eigen_value(Lesliematrix=A)
s = extract_second_eigen_value(Lesliematrix=A)
# This is Loo trying to adjust lambda to equal 1 for each Leslie matrix
AA = A[1,]*Multiplier # multiply fecundities by 1/lambda
A[1,] <- AA # replace the first row, will new fecundities
fadj = extract_first_eigen_value(Lesliematrix=A) # recalc lambda 1
sadj = extract_second_eigen_value(Lesliematrix=A)
c(f,fadj,sadj)

multiplier_values[15] <- 1/8.1


# ---
# Population: codNames[16] = 3Ps
# Test history:
Multiplier = 1/1.10556880817574
Multiplier = 1/2
Multiplier = 1/2.5

A = read.table(paste("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/k1/",
                     codNames[16],".txt",sep=""))
f = extract_first_eigen_value(Lesliematrix=A)
s = extract_second_eigen_value(Lesliematrix=A)
# This is Loo trying to adjust lambda to equal 1 for each Leslie matrix
AA = A[1,]*Multiplier # multiply fecundities by 1/lambda
A[1,] <- AA # replace the first row, will new fecundities
fadj = extract_first_eigen_value(Lesliematrix=A) # recalc lambda 1
sadj = extract_second_eigen_value(Lesliematrix=A)
c(f,fadj,sadj)

multiplier_values[16] <- 1/2.5

write.csv(multiplier_values,file='C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/multiplier_values.csv')
