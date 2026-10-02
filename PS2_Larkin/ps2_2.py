#question 2

import numpy as np
import matplotlib.pyplot as plt
from astrotools import constants as c
from astrotools.nbody import twobody as tb
from astrotools.nbody import integrators as ig

'''
#Question a
'''

a = c.AU #semimajor axis (m)
e = 0.3 #eccentricity
m_star = c.M_SUN #mass of star (kg)
mu = m_star*c.G #mass of star * gravitational constant (m^3/s^2)

P = tb.orbital_period(a) #period of orbit of small body around big body, (s)
r_p = a*(1-e) #distance at perihelion (m)
r_a = a*(1+e) #distance at aphelion (m)
v_p = np.sqrt(mu*((2/r_p)-(1/a))) #velocity at perihelion, (m/s)
v_a = np.sqrt(mu*((2/r_a)-(1/a))) #velocity at aphelion, (m/s)
E = (-1)*mu/(2*a) #specific orbital energy
j = (2*np.pi*(a**2)*np.sqrt(1-e**2))/P #angular momentum, (m^2/s)

#print(P) #correct, 1 year

perihelion = r_p*v_p  #m^2/s
aphelion = r_a*v_a    #m^2/s

if (perihelion == aphelion == j):
    print("Success! j = v_p*r_p = v_a*r_a")


#print(f"{r_p*v_p:e}") #correct
#print(f"{r_a*v_a:e}") #correct
#print(f"{j:e}") #correct

#print(f"{twobody.hill_radius(a, c.M_EARTH, c.M_SUN, 0.0167):e}") #this checked my hill's radius for Earth, seemed to work
#print(twobody.roche_density(1.855e8, 5.6834e26)) #this checked my roche density equation for saturn at mimas's orbit, seemed to work

'''
Question b
'''

y0 = tb.initial_conditions(a, e) #initial coords
dt = P/1000 #defines dt as period / 1000
t, y = ig.integrate(ig.rk4_step, tb.kepler_derivs, y0, (0.0, P), dt) #time t and position+velocity y
err = np.hypot(*(y[-1, :2] - y0[:2])) / a #defines error as in Jupyter
print(f"Error in returned position vs original position is {err}") #prints error in final coords compared to semimajor axis a

drift = tb.relative_drift(tb.specific_energy(y))
print(f"Error in the specific energy at final time vs original time is {drift[-1]}") #prints energy drift at final time




'''
#everything below borrowed from the Jupyter Notebook and modified as needed
nu = np.linspace(0, 2 * np.pi, 400)
r_exact = a * (1 - e ** 2) / (1 + e * np.cos(nu))
xy_exact = np.array([r_exact * np.cos(nu), r_exact * np.sin(nu)]) / c.AU
#above is exact ellipse with eccentricity e

fig, axes = plt.subplots(1, 4, figsize=(13, 4.5))

for ax, n in zip(axes, [50, 500, 5000, 50000]):
    t, y = ig.integrate(ig.euler_step, tb.kepler_derivs, y0, (0.0, P), P / n)
    ax.plot(*xy_exact, 'k--', lw=1, label='exact')
    ax.plot(y[:, 0] / c.AU, y[:, 1] / c.AU, '-', lw=1.2, label='Euler')
    ax.plot(0, 0, 'y*', ms=14)
    ax.set_aspect('equal')
    ax.set_title(f'$\\Delta t = P/{n}$')
    ax.set_xlabel('$x$ (AU)')

for ax, n in zip(axes, [50, 500, 5000, 50000]):
    t_r, y_r = ig.integrate(ig.rk4_step, tb.kepler_derivs, y0, (0.0, P), P / n)
    ax.plot(y_r[:, 0] / c.AU, y_r[:, 1] / c.AU, '-', lw=1.2, label='rk4')
    ax.set_aspect('equal')
    ax.set_title(f'$\\Delta t = P/{n}$')
    ax.set_xlabel('$x$ (AU)')

axes[0].set_ylabel('$y$ (AU)')
axes[0].legend(frameon=False, loc='upper left')
plt.tight_layout()
#plt.show()

rk4 = np.linalg.norm(y_r)
euler = np.linalg.norm(y)

diff = rk4/euler

print(diff)

'''



#(c) 

print(f"Angular momentum for this orbit is {j:e}")


'''Order goes:
        fastest: fiducal core
                 current orbit
                 break-up value
        This suggests that the core was spinning very fast, slowed down to current levels, and then
        will eventually break up once it slows to the break-up value.
'''