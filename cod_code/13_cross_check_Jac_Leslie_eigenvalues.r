# Check if eigenvalues match Leslie and Jacobian matricies
# last updated: 2018 Sept 24

# ****
# Take home: when fecundity in Jacobian is not relative, then eigenvalues match
# ****

# Make sure to run 12_k_vs_CV.rmd to get A3dlist (this has the Jacobian matricies)
# Jacobian matricies in A3dlist should be for when k=1

# Leslie matricies are stored in C:\Users\provo\Documents\GitHub\popdy\cod_code\mikaelaLeslie\matrix_maxages

# load functions to calculate eigenvalues
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/2_cod_functions.r") # load functions

# Northsea
extract_first_eigen_value(Ls3dlist$Northsea[,,1])
A = read.table("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/matrix_maxages/Northsea.txt") # load functions")
extract_first_eigen_value(A)

sum(Ls3dlist$Northsea[1,,1])
# Coas
extract_first_eigen_value(Ls3dlist$Coas[,,1])
A = read.table("C:/Users/provo/Documents/GitHub/popdy/cod_code/mikaelaLeslie/matrix_maxages/Coas.txt") # load functions")
extract_first_eigen_value(A)


