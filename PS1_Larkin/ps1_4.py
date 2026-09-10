import numpy as np
import matplotlib.pyplot as plt
from astrotools import constants as c
from astrotools.cloud import collapse
from astrotools import data as d

N2 = 622 #number of class II stellar objects
t2 = 2e6 #lifetime of class 2 star, years

cores = d.load_herschel_cores() #makes cores array
age = 'prestellar' 
filter = cores['core_type'] == age #makes filter for only prestellar cores
p_cores = cores[filter] #makes new array of only prestellar cores
#print(p_cores['core_type']) #check to see if it worked, it did :)
p_cores['n_H2_avd_cm3'] = p_cores['n_H2_avd_cm3'] * 1e6 #changes cm^-3 to m^-3

"""
Need to set up a loop or something starting at 1e10 and running until 1e12, with changeable number of bins. 
Evaluate t_life at each, which involves filtering to find number of prestellar cores with higher density.
Evaluate t_ff at each, only involves sending density to collapse.free_fall_time(). This needs mass density, however, which involves
    sending number density to collapse.number_density_to_mass_density(), and then sending the reciept to get t_ff.
Then plot t_life on y and t_ff on x, with logarithmic scales. Add line of slope 1 as a reference
Fit the plot, slope should give t_life/t_ff.
"""

t_life = [] #empty array of t_life
t_ff = [] #empty array of t_ff

for n in np.linspace(1e10, 1e12, 1000, endpoint=True): #inclusive linspace with customizable number of iterations
    filter = p_cores['n_H2_avd_cm3'] > n 
    p_cores = p_cores[filter] #filters p_cores to only have entries with density greater than n m^-3
    t_life.append((p_cores.size / N2) * t2) # calculates and adds t_life to grouping
    rho = collapse.number_density_to_mass_density(n) #generates mass density based on threshold number density
    t_ff.append(collapse.free_fall_time(rho)) #calculates and adds t_ff to grouping

t_ff = np.array(t_ff)
t_life = np.array(t_life) #online sources suggested this to allow to make the polynomial.

t_life_min = np.min(t_life) #minimum t_life, for start of line with slope of 1
t_ff_min = np.min(t_ff) #minimum t_ff, for start of line of slope 1
t_diff = t_ff_min - t_life_min #difference between t_ff and t_life at minimum times, setting y-intercept of linear line. 


plt.scatter(t_ff, t_life, label = 'inferred lifetime vs free-fall time') #scatterplot of both log axes
plt.xlabel('t_ff (s)') #labels x
plt.ylabel('t_life (s)') #labels y
plt.title('t_life vs t_ff') #title
m, c = np.polyfit(t_ff, t_life, 1) #first order polynomial fit of the data, slope m.
plt.plot(t_ff, m*t_ff + c, color="black", linestyle = '--', label = f"linear fit of t_life vs t_ff, slope = {m}")
plt.axline((t_ff_min,t_life_min), slope = 1, color = 'r', linestyle = '-', label = 'line of slope 1, starting at minimum times')
plt.xscale('log')
plt.yscale('log')
#plt.plot(t_ff, t_ff - t_diff, label = 'slope of 1') #this is a linear line, just looks very funny since it's on a loglog plot.
#the above line for some reason was producing a very strange line, but the one above it using axline does it correctly. The other one was not linear.
plt.legend()
plt.grid(True) #gridlines to help see.
print(f"The slope of the dotted LSRL is {m}.")

mean_ff = np.mean(t_ff)
mean_life = np.mean(t_life)
ratio = mean_life / mean_ff
print(f"The ratio between mean inferred lifetime and freefall time t_life/t_ff = {ratio}.")

plt.tight_layout()
plt.show()