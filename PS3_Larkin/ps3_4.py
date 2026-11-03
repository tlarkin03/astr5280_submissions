import numpy as np
import matplotlib.pyplot as plt
from astrotools import data as d
from astrotools import constants as c
from astrotools.disk import structure as s
from astrotools.disk import viscous as v

lupus = d._load_table("ps3/lupus_disks.csv","Lupus Disks")

#a

lupus_f = lupus[lupus['F890_mJy'] != 0] #fluxes
lupus_mf = lupus_f[lupus_f['log_mdot'] != 0] #accretions
lupus_fl_mf = lupus_mf[lupus_mf['flag'] == ''] #nonflagged disks
lupus_flag = lupus_mf[lupus_mf['flag'] != ''] #all flagged disks
lupus_chrom = lupus_flag[lupus_flag['flag'] == 'chromospheric'] #chromospheric disks
lupus_upper = lupus_flag[lupus_flag['flag'] == 'mdot_upper_limit']
lupus_sublum = lupus_flag[lupus_flag['flag'] == 'subluminous']
lupus_trans = lupus_flag[lupus_flag['flag'] == 'transitional']

kappa = 0.337 #m^2/kg dust
v = c.C_LIGHT / 8.9e-4 #frequency at 890 micrometers, in Hz
T_dust = 20 #K
B_v = s.planck_nu(v,T_dust) #W / m^2-Hz-sr

dist = lupus_fl_mf['distance_pc'] * c.PC #parsecs of each disk to m
dist_flag = lupus_flag['distance_pc'] * c.PC
dist_chrom = lupus_chrom['distance_pc'] * c.PC
dist_upper = lupus_upper['distance_pc'] * c.PC
dist_sublum = lupus_sublum['distance_pc'] * c.PC
dist_trans = lupus_trans['distance_pc'] * c.PC

Flux = lupus_fl_mf['F890_mJy'] * 10**-29 #flux of each disk, mJy, to W/m^2-Hz
Flux_flag = lupus_flag['F890_mJy'] * 10**-29
Flux_chrom = lupus_chrom['F890_mJy'] * 10**-29
Flux_upper = lupus_upper['F890_mJy'] * 10**-29
Flux_sublum = lupus_sublum['F890_mJy'] * 10**-29
Flux_trans = lupus_trans['F890_mJy'] * 10**-29

logMdot = lupus_fl_mf['log_mdot'] #log mdot
Mdot = 10**logMdot #mdot in solar masses/year
Mdot_flag = 10**lupus_flag['log_mdot']
Mdot_chrom = 10**lupus_chrom['log_mdot']
Mdot_upper = 10**lupus_upper['log_mdot']
Mdot_sublum = 10**lupus_sublum['log_mdot']
Mdot_trans = 10**lupus_trans['log_mdot']

M_dust = Flux * dist**2 / (kappa * B_v * T_dust) #according to 2b, kg
M_dust_sun = M_dust / c.M_SUN #kg to M_sun
M_disk = 100 * M_dust

Mdust_flag = (Flux_flag * dist_flag**2 / (kappa * B_v * T_dust)) / c.M_SUN
Mdisk_flag = 100 * Mdust_flag
Mdisk_chrom = 100 * ((Flux_chrom * dist_chrom**2 / (kappa * B_v * T_dust)) / c.M_SUN)
Mdisk_upper = 100 * ((Flux_upper * dist_upper**2 / (kappa * B_v * T_dust)) / c.M_SUN)
Mdisk_sublum = 100 * ((Flux_sublum * dist_sublum**2 / (kappa * B_v * T_dust)) / c.M_SUN)
Mdisk_trans = 100 * ((Flux_trans * dist_trans**2 / (kappa * B_v * T_dust)) / c.M_SUN)

plt.scatter(Mdot, M_disk, color = 'black', label = 'non-flagged')
#plt.scatter(Mdot_chrom, Mdisk_chrom, color = 'blue', label = 'Chromospheric', marker = '*') #chromospheric redundant
#plt.scatter(Mdot_upper, Mdisk_upper, color = 'red', label = 'Upper Limit', marker = '*') #upper limit redundant
plt.scatter(Mdot_sublum, Mdisk_sublum, color = 'yellow', label = 'Subluminous', marker = '*')
plt.scatter(Mdot_trans, Mdisk_trans, color = 'lime', label = 'Transitional', marker = '*')
plt.yscale('log')
plt.ylabel('M_disk, solar masses')
plt.xscale('log')
plt.xlabel('Mdot, solar masses / yr')
plt.legend()
plt.title('Mdot vs M_disk for Lupus Disks')
#plt.show()

#b

t_acc = M_disk / Mdot
median = np.median(t_acc)
print (f"The median of t_acc is {median} years.")