# popdy
Sensitivities of west coast fish stocks to environmental variability 

`simulation_model.rmd` - this is the main script that generates timeseries of eggs, spawners, and catch. in here it reads two other scripts: parameters_pacifichake.r and functions.r

output is: a list called 'output'. list is made up of lists. each list has timeseries for eggs, spawners, and catch associated with a specific fishing mortality. 

`parameters_pacifichake.r` - parms specific for hake. need to figure out how to write code for other species

`functions.r` - these functions are general, not species specific. but will need to review them if life history is iteroperous/semelparous. 

`plot_timeseries.r` - creates a series of plots. right now it's set up to print plots only associeted with one species. 