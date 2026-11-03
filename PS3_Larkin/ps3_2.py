
import numpy as np
import matplotlib.pyplot as plt
from astrotools import constants as c
from astrotools.disk import structure as s
from astrotools.disk import viscous as v

T_star = 4000 #Temperature of star is 4000 K
R = 2*c.R_SUN #radius of star is 2* sun radius [m]

M = c.M_SUN #mass of star is 1 solar mass [kg]
d = 150 * c.PC #distance of 150 parsecs, converted to [m]
cos_i = 1 #inclination is 0, so cos_i = 1. Don't need to send to s.sed

wave = np.linspace(1e-6,0.003,2000) #wavelength array from 1e-6 m to 0.003 m, same number of entries as Temp
wave_cm = wave * 100 #wavelength array in cm, for plotting
nu = c.C_LIGHT / wave #array of frequency from wavelengths

r_in = 0.1 * c.AU #edge of 0.1 AU to [m]
r_out = 100 * c.AU #edge of 100 AU to [m]

edges = s.radial_grid(r_in,r_out) #m, 2001 edges of 2000 annuli
diff_edges = np.diff(edges) #centers of each annulus, redundant? [m]
T_flat = s.temperature(diff_edges,81,-0.75) #temp, flat [K]
T_flared = s.temperature(diff_edges,160,-3/7) #temp, flared [K]


#a

''' Compute nu*F_nu from 1 micrometer to 3mm, for optically thick, 0.1-100 AU. Do both flat disk and flared disk. 
        - nu is set by the wavelength, from 1 to 3000 micrometers. (f = c/lambda)
Plot photosphere nu_F_nu_photo, disk, and the sum on one figure.
Measure alpha_ir between 2.2 and 24 micrometers for each spectrum, compare with power-law value.

'''
thick_flat = s.sed(nu,edges,T_flat,d) #flat sed
thick_flared = s.sed(nu,edges,T_flared,d) #flared sed
F_nu_photo = np.array([])
for v in range(len(nu)):
        F_nu_photo = np.append(F_nu_photo, np.pi * R**2 * s.planck_nu(v,T_star) / d**2) #photosphere sed

nu_flat = nu * thick_flat #nu fnu
nu_flared = nu * thick_flared #nu fnu
nu_photo = nu * F_nu_photo #nu fnu

flat_sum = nu_flat + nu_photo #sum of flat and photosphere
flared_sum = nu_flared + nu_photo #sum of flared and photosphere


plt.plot(wave_cm,nu_flat, color = 'red', label = 'Flat')
#plt.plot(wave_cm,flat_sum, color = 'red', linestyle = '--', label = 'Flat + Photosphere')
plt.plot(wave_cm,nu_flared, color = 'blue', label = 'Flared')
#plt.plot(wave_cm,flared_sum, color = 'blue', linestyle = '--', label = 'Flared + Photosphere')
plt.plot(wave_cm,nu_photo, color = 'black', linestyle = '--', label = 'Photosphere')
plt.xlabel('λ (cm)')
plt.xscale('log')
plt.ylabel(r'$νF_ν$ $(W/m^2)$')
plt.yscale('log')
plt.title(r'$νF_ν$ for Optically Thick Flared and Flat Disks')
plt.legend()

#f_lambda



plt.show()


#b

''' Flared disk with minimum-mass surface density, 0.1-30 AU, dust-to-gas = 0.01, opacity of s.dust_opacity().
- F_nu at 1.3 mm, and dust mass at T_dust = 20 K
- recompute flux as if optically thin, I_nu = tau_nu * B_nu, compare.
'''
'''

r_in = 0.1*c.AU
r_out = 30*c.AU
r = s.radial_grid(r_in,r_out)
sigma = s.surface_density(r)
nu = c.C_LIGHT / (1.3 / 1000) #meters wavelength to hertz
kappa = s.dust_opacity(nu)

#where does dust to gas ratio come in?

F_nu_dust = s.sed(nu,r,T_flared,d,sigma,kappa)

'''