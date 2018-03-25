# Parms - Pacific Mackerel
# by: Mikaela Provost
# on: 2017-3-1

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

# Stock assessment: 
# Crone, P. R., Hill. K. T. 2015. Pacific mackerel (Scomber japonicus) stock assessment for USA  
# management in the 2015-16 fishing year. Pacific Fishery Management Council, Pacific Fishery 
# Management Council, 7700 NE Ambassador Place, Suite 101, Portland, Oregon 97220, USA. 131 p.

# ---
# 1. Define max age
# maximum age (yr) from table from John Field = 14
# max age from 2015 assessment = 12
a_max = 14

# ---
# 2. Growth Parameters:
# using values for 1983-2014 in 2015 assessment, see pg 16,34
# von Bertalanffy growth coefficient
k = 0.39 # units = time^-1

# von Bertalanffy theoretical maximum length (cm), see pg 16
L_inf = 39.2 

# allometric length-weight exponent
b = 3.4 # pg.16, b = 3.394 in FishBase

# allometric length-weight scale from FishBase is 0.00137, 
# but value from assessment is 0.0000027 (this is in kig)
L_W_slope = 0.00137 

# von Bertalanffy theoretical maximum weight (g)
W_inf = L_W_slope * (L_inf ^ b)

# von Bertalanffy theoretical time at L = 0 or W = 0
# useful resource for calculating t_0: 
#http://derekogle.com/NCNRS349/modules/Growth/BKG.html
# use the Schnute parameterization of the vonB because that's what they use in assessment
L_min_age = 20.5
L_max_age = 39.3
t1 = 0.5
t3 = 12

t_0 = t1 + (1/k)* log((L_max_age - L_min_age)/(L_max_age-L_min_age*exp(-k*(t3-t1)))) 
###t_0 = - 0.8704 # because lgth vs wt curve does not go to origin, how did Lewis get this value?
#t_0 = 0 # assume t_0=0 for now, not sure how to estimate this value

# calculate von Bertalanffy weight at t = 0, W_0, from t_0
#W_0 = (W_inf * (1 - exp(k * t_0))) ^ b

# ---
# 3. Size at 50% maturity parms:
# size at 50% maturity (cm) from Crone in 2015 assessment, p17
L_mat = 27
W_mat = L_W_slope * (L_mat ^ b)
# standard deviation of size 50% maturity
L_mat_sd = 3 # ********** this is a place holder for now
W_mat_sd = (L_mat_sd / L_mat) * W_mat #is_mature


# --- 
# 4. Size at recruitment parms:
# mean and sd of size at first recruitment to fishery
# According to assessment, age at recruitment is age 0 (p36)
L_rec = 20.5 # length-at-age_minimum (p16)
L_rec_sd = 0.10 # SD was fixed to 0.10 based on little info 
W_rec = 3 * (L_W_slope * (L_rec ^ b))
W_rec_sd = 2 * ((L_rec_sd / L_rec) * W_rec)
# if dome-shaped selectivity desired, could specify mean and sd of size at escape from fishing
# assume asymptotic selectivity for now - I set the length of escapement really high
# so that basically call fish do not escape
# assessments says dome shape selectivity used for recreational fishery & dome for commerical
# not sure what to make of that.
L_esc = 100 # assume this for now
L_esc_sd = 1 # assume this for now
W_esc = L_W_slope * (L_esc ^ b)
W_esc_sd = (L_esc_sd / L_esc) * W_esc


# --- 
# 5. Natural mortality:
# p33 in assessment
natural_mort = 0.1 # for all ages, constant over time

# ---
# 6. Fecundity-weight slope:
# a value of 1 assumes fecundity ~ weight
f_W_slope = 1
# fecundity-weight intercept
# for now, assume 0. Look into assessment to see if different
f_W_int = 0


# ---
# 7. BH stock-recruitment parms:
# Beverton-Holt stock-recruitment parameters, given form R = (E * alpha) / (1 + (beta * E))

steepness = 0.48 # p36
R_0 = 0.54E9 # p36   
alpha = ((4 * steepness) / (1 - steepness)) * 0.001
beta = alpha / R_0


# ---
# 8. Magnitude of recruitment variation:
# magnitude of recruitment variation, in terms of standard deviation, 
# table from John says 1, which matches older assessment 
eta_sd = 0.75 # p34

# ---
# 9. Magnitude of background mortality variation:
# magnitude of background mortality variation, in terms of standard deviation, given mean of 0
# Not sure where these numbers came from (this is from Lewis)
###M_sd = (39.93565/71.79419) * natural_mort


# ---
# ********* Have not looked into these parms yet ***********
# 10. Magnitude of growth variation
# in terms of standard deviation
# based on observed sd accross fits per 2-year period for female hake between 1975 and 2010 from 2011 assessment
# note mean L_inf from same data was 71.79419 in case want to use this instead of overall L_inf estimate as the mean
#L_inf_sd = 39.93565
#L_inf_sd = (39.93565/71.79419) * L_inf
#W_inf_sd = (39.93565/71.79419) * W_inf
# include parameters for bounding L_inf deviates here? 2011 assessment data indicates lowest L_inf in the high 40s
#L_inf_min = 40
#W_inf_min = L_W_slope * (L_inf_min ^ b)
# maximum size reported (for a male) from fishbase -> Cohen et al. 1990
#L_inf_max = 91
#W_inf_max = L_W_slope * (L_inf_max ^ b)


