import numpy as np
import matplotlib.pyplot as plt
from astrotools import constants as c
from astrotools.cloud import collapse
from astrotools import data as d

cores = d.load_herschel_cores() #makes cores array
age = 'prestellar' 
filter = cores['core_type'] == age #makes filter for only prestellar cores
p_cores = cores[filter] #makes new array of only prestellar cores
#print(p_cores['core_type']) #check to see if it worked, it did :)

#need to use n_H2_avd_cm3 and T_dust for Mj calcs. n_H2_avd_cm3 is in cm^-3, need to convert to m^-3 then use collapse.number_density_to_mass_density()
#print(p_cores['n_H2_avd_cm3'][1])
p_cores['n_H2_avd_cm3'] = p_cores['n_H2_avd_cm3'] * 1e6 #changes cm^-3 to m^-3
#print(p_cores['n_H2_avd_cm3'][1])
p_cores['n_H2_avd_cm3'] = collapse.number_density_to_mass_density(p_cores['n_H2_avd_cm3']) #changes number density to mass density
#print(p_cores['n_H2_avd_cm3'][1])
mj = collapse.jeans_mass(p_cores['T_dust_K'],p_cores['n_H2_avd_cm3']) #computes Jean's mass for all cores
#print(mj.size) #check to make sure it worked.

#both histograms means one for predicted jean's mass, one for actual mass of the core
#histograms logarithmically spaced, so .1-.9, 1-9, 10-90

mj = mj * 1.989 * 1e-30 #converts jeans masses from kg to solar masses

xmin = 10**-1.5 #minimum boundary for histogram
xmax = 10**1.5 #maximum boundary for histogram
logbins = np.geomspace(xmin, xmax, num=16) #creates logarithmic bins for histograms, lower bound at 10^-1.5 solar masses, upper at 10^1.5 solar masses.
figure, (ax0, ax1) = plt.subplots(1, 2) #creates two plots, one row two columns, called ax0 and ax1

mj_counts, bins, patches = ax0.hist(mj, bins=logbins) #create histogram for jeans mass, logarithmic bins
ax0.set_title('Jeans Mass of Prestellar Cores') 
ax0.set_xlabel('log(mass) [M_sun]')
ax0.set_xscale('log')
ax0.set_ylabel('Quantity')

max_mj = mj_counts.max() #quantity in most frequent bin
peak_mj = np.argmax(mj_counts) #finds index of the bin where Mj peaks
peak_mj_start = bins[peak_mj] #left bound of peak index
peak_mj_end = bins[peak_mj+1] #right bound of peak index
mj_mode = (peak_mj_start + peak_mj_end) / 2

ax0.axhline(max_mj, color = 'black', linestyle = '--') #marks highest quantity with horizontal line
ax0.axvline(mj_mode, color = 'black', linestyle = '--') #marks mode mj mass with vertical line


print(f"Jeans mass mode: {mj_mode} solar masses.")

cores_counts, bins, patches = ax1.hist(p_cores['M_core_msun'], bins=logbins) #creates histogram for core masses
ax1.set_title('Masses of Prestellar Cores')
ax1.set_xlabel('log(mass) [M_sun]')
ax1.set_xscale('log')
ax1.set_ylabel('Quantity')

max_cores = cores_counts.max()
peak_cores = np.argmax(cores_counts) #find index of bin where core mass peaks
peak_cores_start = bins[peak_cores] #left bound of peak index
peak_cores_end = bins[peak_cores+1] #right bound of peak index
cores_mode = (peak_cores_start + peak_cores_end) / 2 #core mode mass (approximated)


ax1.axhline(max_cores, color = 'black', linestyle = '--') #marks most frequent bin with horiontal line
ax1.axvline(cores_mode, color = 'black', linestyle = '--') #marks mode core mass with vertical line

print(f"Core mode: {cores_mode} solar masses.")

mj_median = np.median(mj)
cores_median = np.median(p_cores['M_core_msun'])

median_ratio = (mj_median / cores_median)
mode_ratio = (mj_mode / cores_mode)

print(f"Ratio between histogram peaks Mj/M: {mode_ratio}")
print(f"Ratio between median masses Mj/M: {median_ratio}")

#plt.legend()
plt.tight_layout() #prevents text from getting to close to the plots
plt.show() #displays plots