# Find oldest fish in each population & overall oldest
# by: Mikaela Provost

# Motiviation:
# Our spawning biomass distributions have very long tails and
# these tails might not be biologically meaningful. This script
# digs into the data to find the oldest fish observed in each 
# population. This information will be used to truncate the 
# spawning biomass distributions to the age of the oldest 
# observed fish. 

# ---
# load cod data, break into separate populations
source("C:/Users/provo/Documents/GitHub/popdy/cod_code/0_load_cod_data.r")

# packages
library(ggplot2)
library(gridExtra)
# one vector for max ages and one for occurences of max age in each pop
max_ages <- rep(NA, length(codNames))
max_ages_occurances <- rep(NA, length(codNames))

for (i in 1:length(codNames)) {
  pop <- datalist[[i]]
  max_ages[i] <- max(pop$AGE)
  max_ages_occurances[i] <- length(pop[pop$AGE == max_ages[i],]$AGE) 
}

max_ages_table<- data.frame(cbind(codNames, max_ages, max_ages_occurances))
max_ages_table$max_ages <- as.numeric(as.character(max_ages_table$max_ages))
max_ages_table$max_ages_occurances <- as.numeric(as.character(max_ages_table$max_ages_occurances))


# plot abundance at age for every year. 
# x axis --> age, y --> abundance (CANUM), one line for each year
#pdf(file='C:/Users/provo/Documents/GitHub/popdy/cod_figures/abundnace_at_age.pdf', width=20, height=30)
#p <- list()
#for (i in 1:length(codNames)) {
#  pop <- datalist[[12]]
#  pop <- pop[pop$AGE > 5,] # plot abundance at older ages to see contrast better
#  pop <- pop[pop$YEAR %in% c('1965','1970','1975','1980','1985','1990','1995','2000','2005'),]
    # plot only a few years so it isn't as messy
  #pop$CANUM <- as.numeric(as.character(pop$CANUM)) # convert abundances to numeric
  #pop$YEAR <- as.numeric(as.character(pop$YEAR))
  #pop$YEAR <- factor(pop$YEAR)

#  p[[i]] <- ggplot(data=pop, aes(x=AGE, y=CANUM, color=YEAR)) + 
#            geom_line() + 
#            theme(legend.position='none') +
#            ggtitle(codNames[12]) +
#            scale_colour_gradient2()
#}
#do.call(grid.arrange,p)
#rm(p,i)
#dev.off()

# ===> this wasn't working, so trying bubble plots. 
# ===> problem: I want no bubble where there is zero abundance at a particular age
# bubble plot: x axis --> year, y axis -- > age, bubble size corresponds to abundance (CANUS)
#pop$CANUM <- as.numeric(as.character(pop$CANUM)) # convert abundances to numeric
#ggplot(data=pop, aes(YEAR,AGE, size=CANUM)) +
#  geom_point() + scale_radius()
