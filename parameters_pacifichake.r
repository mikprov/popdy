# Parms - Pacific Hake
# by: Mikaela Provost
# on: 2017-8-16

# ---
# goals:
# ---
#1. define max age for this species
#2. define growth parms
#3. define size at 50% maturity
#4. define size at first recruitment parms
#5. define natural mortality
#6. define fecundity-weight slope
#7. define BH stock-recruitment parms
#8. define magnitude of recruitment variation
#9. define magnitude of background mortality variation
#10. define magnitude of growth variation

# Parameters from Lewis
# Current values are for Pacific hake, as in Botsford et al. 2014 and from assessments

# ---
# 1. Define max age
# maximum age (yr) from assessment, whereas Botsford et al. use 22
#a_max = 20 #Lewis
a_max = 22

# ---
# 2. Growth Parameters:
# using values for 1990-2011 in 2011 assessment
# von Bertalanffy growth coefficient
k = 0.3720 # Lewis 
#k = 0.3002 # from Botsford et al. 2014

# von Bertalanffy theoretical maximum length (cm)
#L_inf = 51.2346 # Lewis
L_inf = 53.36 # Botsford et al. 2014

# allometric length-weight exponent
b = 2.9624 # same in Botsford et al. 2014

# allometric length-weight scale is 0.006419 in Botsford et al. 2014, 
# but value from assessment is 0.000007 (maybe in kg there?)
L_W_slope = 0.006419 # same in Botsford et al. 2014

# von Bertalanffy theoretical maximum weight (g)
W_inf = L_W_slope * (L_inf ^ b)

# von Bertalanffy theoretical time at L = 0 or W = 0
t_0 = - 0.8704 # because lgth vs wt curve does not go to origin, how did Lewis get this value?

# calculate von Bertalanffy weight at t = 0, W_0, from t_0
W_0 = W_inf * (1 - exp(k * t_0)) ^ b


# ---
# 3. Size at 50% maturity parms:
# size at 50% maturity (cm) from Dorn in 2014 assessment
L_mat = 36.8900
W_mat = L_W_slope * (L_mat ^ b)
# standard deviation of size 50% maturity
L_mat_sd = 3.5395
W_mat_sd = (L_mat_sd / L_mat) * W_mat


# --- 
# 4. Size at recruitment parms:
# mean and sd of size at first recruitment to fishery
# Length-based values from Botsford et al. 2014, but weight-based values
# were tuned to get close to matching selectivity at age from 2014 assessment
L_rec = 26.06 # same as Botsford et al. 2014
L_rec_sd = 9.401 # same as Botsford et al. 2014
W_rec = 3 * (L_W_slope * (L_rec ^ b))
W_rec_sd = 2 * ((L_rec_sd / L_rec) * W_rec)
# if dome-shaped selectivity desired, could specify mean and sd of size at escape from fishing
L_esc = 51.71
L_esc_sd = 1.173
W_esc = L_W_slope * (L_esc ^ b)
W_esc_sd = (L_esc_sd / L_esc) * W_esc


# --- 
# 5. Natural mortality:
# instantaneous rate of natural mortality from 2014 initial value =0.2, was used in first runs,
#   but MLE=0.213 and posterior median=0.222; whereas Botsford et al. 2014 use 0.23, and 0.23 sometimes used elsewhere
# natural_mort = 0.222 # Lewis
natural_mort = 0.23 # Botsford et al. 2014

# ---
# 6. Fecundity-weight slope:
# fecundity-weight slope from Botsford et al. 2014 is 0.0009638, assessment is 1?
# to get it to produce reasonable fecundity estimates, slope should be orders of magnitude higher
# but given that alpha is set for S-R rather than E-R relationship, scaling of S-E-R goes into fecundity
# f_W_slope = 0.0009638 * 1000000
# now value is arbitrary, but gives somewhat reasonable value of F_collapse of 0.815 given old parameters, or 0.46 given new ones
# NEVERMIND, PROBABLY MAKES MORE SENSE TO PUT SCALING INTO ALPHA INSTEAD OF F_W_SLOPE, SO USE ASSESSMENT VALUE HERE
f_W_slope = 1
# fecundity-weight intercept from Botsford et al. 2014, assessment says 1, but I think it is zero?
f_W_int = 0


# ---
# 7. BH stock-recruitment parms:
# Beverton-Holt stock-recruitment parameters, given form R = (E * alpha) / (1 + (beta * E))
# First runs used 2014 assessment initial values for steepness = 0.88 and maximum recruitment ln(R_0) = 15.9
# consider using prior steepness which came from Myers et al. 1999, h = 0.777 (Gadidae) because otherwise too high
# otherwise, try posterior medians: steepness=0.826, M=0.22, R_0=2.72E9
steepness = 0.826
R_0 = 2.72E9
alpha = ((4 * steepness) / (1 - steepness)) * 0.001
beta = alpha / R_0
 

# ---
# 8. Magnitude of recruitment variation:
# magnitude of recruitment variation, in terms of standard deviation, given mean of 0, assessment says 1.4
# also, see Stachura et al. 2014 if want another, more robust, estimate for Alaska stocks
#eta_sd = (39.93565/71.79419) # Lewis
eta_sd = 0.7 # Botsford et al. 2014

# ---
# 9. Magnitude of background mortality variation:
# magnitude of background mortality variation, in terms of standard deviation, given mean of 0
M_sd = (39.93565/71.79419) * natural_mort


# ---
# 10. Magnitude of growth variation
# in terms of standard deviation
# based on observed sd accross fits per 2-year period for female hake between 1975 and 2010 from 2011 assessment
# note mean L_inf from same data was 71.79419 in case want to use this instead of overall L_inf estimate as the mean
#L_inf_sd = 39.93565
L_inf_sd = (39.93565/71.79419) * L_inf
W_inf_sd = (39.93565/71.79419) * W_inf
# include parameters for bounding L_inf deviates here? 2011 assessment data indicates lowest L_inf in the high 40s
L_inf_min = 40
W_inf_min = L_W_slope * (L_inf_min ^ b)
# maximum size reported (for a male) from fishbase -> Cohen et al. 1990
L_inf_max = 91
W_inf_max = L_W_slope * (L_inf_max ^ b)

 
