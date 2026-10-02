#Question 4

import numpy as np
import matplotlib.pyplot as plt
from astrotools import constants as c
from astrotools.nbody import twobody as tb
from astrotools.nbody import integrators as ig
from astrotools import data as d

#functions

def rigid(d, m_planet):
    #d = distance from satellite to planet, m
    #m_star = mass of planet, kg

    #returns critical density, kg/m^3

    return (9 * m_planet) / (4* np.pi * d**3)

def fluid(d, R_p, rho_p):
    #d = distance from satellite to planet (m)
    #R_p = Radius of planet (m)
    #rho_p = density of planet (kg/m^3)

    #returns critical density, kg/m^3

    return rho_p * ((d / (2.44 * R_p))**-3)


#a

''' Compute r_h for the 8 planets of our solar system. This requires taking their semimajor axis, their mass, and their eccentricities
and passing them to tb.hill_radius(a, m_planet, m_star = c.M_SUN, e). Check Earth's hill radius once complete.
'''
m_star = c.M_SUN #kg
planets = d.load_planet_elements() #name, a_au, e, mass_kg, radius_m
planets['a_au'] = planets['a_au'] * 1.496e+11 #converts AU to m

r_h = tb.hill_radius(planets['a_au'], planets['mass_kg'], m_star, planets['e']) #generates hill radius for each planet, m
r_h = r_h.reshape(1, -1) #makes r_h a single column array
print(r_h)
print(f" Earth's Hill's Radius = {r_h[0,2]:e} m according to this program.") #This is the same as expected!

#b

''' Show that the period of satillite orbiting at r_h is P_sat = P_planet / sqrt(3) using Kepler's law first. The evaluate P_sat at
1/2 r_h and show in units of P_planet
'''


P_sat = tb.orbital_period(0.5*r_h) #produces satellite periods for each planet at 1/2 r_h, in units of m.
P_planet = tb.orbital_period(planets['a_au']) #produces planet periods, m
P_sat_adj = P_sat/P_planet #satellite periods are now represented as a fraction of their host's period
print(P_sat_adj)

#c

''' Load planetary satellites into array from d.load_satellites(). express satellite
semimajor axis in terms of their host planet's r_H. Separate prograde and retrograde orbits (5th column of the spreadsheet says this).
Report largest a/r_h for each planet. This can be expressed as a table or figure of some sort (over 400 moons, quite a lot).
'''

sat = d.load_satellites() #planet, name, a_m, e, i_deg, direction
#only Jupiter, saturn, uranus, and neptune are in here (rows 4,5,6,7 in planets)


for row in sat:
    name = row['planet']
    if name == 'Jupiter':
        row['a_m'] = row['a_m']/r_h[0,4] #if host is jupiter, divide a by jupiter r_h
    if name == 'Saturn':
        row['a_m'] = row['a_m']/r_h[0,5] #if host is saturn, divide a by saturn r_h
    if name == 'Uranus':
        row['a_m'] = row['a_m']/r_h[0,6] #if host is uranus, divide a by uranus r_h
    if name == 'Neptune':
        row['a_m'] = row['a_m']/r_h[0,7] #if host is neptune, divide a by neptune r_h

retrograde_f = sat['direction'] == 'retrograde'
prograde_f = sat['direction'] == 'prograde'

retrograde = sat[retrograde_f]
prograde = sat[prograde_f]

#filter jupiter retrogrades
J_r = retrograde['planet'] == 'Jupiter'
J_retrograde = retrograde[J_r]
#jupiter progrades
J_p = prograde['planet'] == 'Jupiter'
J_prograde = prograde[J_p]
#saturn retrogrades
S_r = retrograde['planet'] == 'Saturn'
S_retrograde = retrograde[S_r]
#saturn progrades
S_p = prograde['planet'] == 'Saturn'
S_prograde = prograde[S_p]
#uranus retrogrades
U_r = retrograde['planet'] == 'Uranus'
U_retrograde = retrograde[U_r]
#uranus progrades
U_p = prograde['planet'] == 'Uranus'
U_prograde = prograde[U_p]
#neptune retrogrades
N_r = retrograde['planet'] == 'Neptune'
N_retrograde = retrograde[N_r]
#neptune progrades
N_p = prograde['planet'] == 'Neptune'
N_prograde = prograde[N_p]

Jr_max = np.max(J_retrograde['a_m'])
Sr_max = np.max(S_retrograde['a_m'])
Ur_max = np.max(U_retrograde['a_m'])
Nr_max = np.max(N_retrograde['a_m'])

Jp_max = np.max(J_prograde['a_m'])
Sp_max = np.max(S_prograde['a_m'])
Up_max = np.max(U_prograde['a_m'])
Np_max = np.max(N_prograde['a_m'])

rmin = np.min(retrograde['a_m'])
rmax = np.max(retrograde['a_m'])
pmin = np.min(prograde['a_m'])
pmax = np.max(prograde['a_m'])

histmin = 1/2500
histmax = 1/2
logbins = np.geomspace(histmin, histmax, num=16)
plt.hist(retrograde['a_m'], bins=logbins, alpha = 0.5, label = 'Retrograde Orbits', color = 'red')
plt.hist(prograde['a_m'], bins=logbins, alpha = 0.5, label = 'Prograde Orbits', color = 'blue')
plt.xscale('log')
plt.xlabel(f'log(a/$r_h$)')
plt.ylabel('Number of Satellites')
plt.title('Proximity to the Hill Radius for Retrograde and Prograde Satellites')
plt.annotate(text="0.5", xy=(0.5, 0), xytext=(0, -17), textcoords="offset points", ha="center") #marks where the ratio is 1/2
plt.legend()
#plt.show()

print(f"The maximum orbit radius for Jupiter prograde orbit in terms of a/r_h is {Jp_max}")
print(f"the maximum orbit radius for Jupiter retrograde orbit in terms of a/r_h is {Jr_max}")
print(f"The maximum orbit radius for Saturn prograde orbit in terms of a/r_h is {Sp_max}")
print(f"the maximum orbit radius for Saturn retrograde orbit in terms of a/r_h is {Sr_max}")
print(f"The maximum orbit radius for Uranus prograde orbit in terms of a/r_h is {Up_max}")
print(f"the maximum orbit radius for Uranus retrograde orbit in terms of a/r_h is {Ur_max}")
print(f"The maximum orbit radius for Neptune prograde orbit in terms of a/r_h is {Np_max}")
print(f"the maximum orbit radius for Neptune retrograde orbit in terms of a/r_h is {Nr_max}")

#d

''' For Saturn's A ring (1.368e8 m), what mean density does both the rigid criterion and fluid criterion require to be placed 
at the edge of the ring? (Set a to be the radius of the ring for rigid, and d to the radius for the fluid criterion). Which one
is plausible for water ice, and what does it imply about the strength of the ring material?
'''

#saturn is row 6 in planets. mass is column 3, radius column 4

d = 1.368e8 #saturn's A ring radius, m
for row in planets: #for some reason I wasn't able to figure out how to index this other than a loop like this
    name = row['name']
    if name == 'Saturn':
        R_p = row['radius_m']
        m_planet = row['mass_kg']

rho_p = m_planet / ((4/3) * np.pi * (R_p**3)) #density of saturn using radius and mass in spreadsheet


rigid = rigid(d, m_planet) #density of rigid criteria
fluid = fluid(d, R_p, rho_p) #density of fluid criteria

print(f"The critical density described by the rigid criteria is {rigid} kg/m^3")
print(f"The critical density described by the fluid critera is {fluid} kg/m^3")

#only fluid criteria makes much sense. at STP, ice is 1000kg/m^3. 750 isn't too bad, but it does mean that the ice is not densely packed
#at all in the rings


#e

''' Where are the final satellites that we have for each planet relative to r_h, and is it different for prograde vs retrograde moons?
Explain why material near r_h is only partially bound: how many orbits do they complete per orbit of their planet? Why does this
number matter?
'''

#census stops around a ratio of 1/2. Retrograde orbits also seem more tolerant of being closer to the hill radius, although 400 satellites
#of 4 planets is not exactly the largest sample size. 

#a satellite at half their host's hill radius orbit somewhere between 150 and 10000 times per orbit of their parent body. 