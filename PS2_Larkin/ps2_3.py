#question 3

import numpy as np
import matplotlib.pyplot as plt
from astrotools import constants as c
from astrotools.nbody import twobody as tb
from astrotools.nbody import integrators as ig

#part a,c

a = c.AU
e = 0.3 #set to either 0.3 or 0.9 for part a or c, respectively.
m_star = c.M_SUN
mu = m_star*c.G
P = tb.orbital_period(a)
y0 = tb.initial_conditions(a, e)

#exact elipse as defined in jupyter
nu = np.linspace(0, 2 * np.pi, 400)
r_exact = a * (1 - e ** 2) / (1 + e * np.cos(nu))
xy_exact = np.array([r_exact * np.cos(nu), r_exact * np.sin(nu)]) / c.AU

fig, axes = plt.subplots(2, 4, figsize=(13, 4.5)) #defines plots, 6 of them
axes = axes.flatten() #flattens grid so that computation can be done on them
err_e = np.array([])
err_r = np.array([])

t1 = 230
t2 = 940
t3 = 3750
t4 = 15000

for ax, n in zip(axes[0:4], [t1, t2, t3, t4]): #first 4 subplots, euler, at dt = 25, 100, 400, 1600
    t, y = ig.integrate(ig.euler_step, tb.kepler_derivs, y0, (0.0, P), P / n)
    err = np.hypot(*(y[-1, :2] - y0[:2])) / a
    err_e = np.append(err_e, err)
    ax.plot(*xy_exact, 'k--', lw=1, label='exact')
    ax.plot(y[:, 0] / c.AU, y[:, 1] / c.AU, '-', lw=1.2, label='Euler')
    ax.plot(0, 0, 'y*', ms=14)
    ax.set_aspect('equal')
    ax.set_title(f'$\\Delta t = P/{n}$')
    ax.set_xlabel('$x$ (AU)')
    ax.set_ylabel('$y$ (AU)')
    ax.legend(frameon=False, loc='upper left')

for ax, n in zip(axes[4:], [t1, t2, t3, t4]): #second 4 subplots, euler, at dt = 25, 100, 400, 1600
    tr, yr = ig.integrate(ig.rk4_step, tb.kepler_derivs, y0, (0.0, P), P / n)
    err = np.hypot(*(yr[-1, :2] - y0[:2])) / a
    err_r= np.append(err_r, err)
    ax.plot(*xy_exact, 'k--', lw=1, label='exact')
    ax.plot(yr[:, 0] / c.AU, yr[:, 1] / c.AU, '-', lw=1.2, label='RK4')
    ax.plot(0, 0, 'y*', ms=14)
    ax.set_aspect('equal')
    ax.set_title(f'$\\Delta t = P/{n}$')
    ax.set_xlabel('$x$ (AU)')
    ax.set_ylabel('$y$ (AU)')
    ax.legend(frameon=False, loc='upper left')


n = np.array([t1, t2, t3, t4])
#plt.clf()
plt.figure(2)
plt.yscale('log')
plt.xscale('log')
plt.plot(1/n, err_e, color = 'red', label = 'Euler Error vs $\\Delta t$/P')
plt.plot(1/n, err_r, color = 'blue', label = 'RK4 Error vs $\\Delta t$/P')
plt.axline(xy1 = (10e-5,10e-4), xy2 = (10e-4, 10e-3), color = 'red', linestyle = '--', label = "Slope: -1")
plt.axline(xy1 = (10e-5,10e-10), xy2 = (10e-4,10e-6), color = 'blue', linestyle = '--', label = "Slope: -4" )
plt.title('Comparison of Error in RK4 and Euler')
plt.xlabel('$\\Delta t$/P')
plt.ylabel('Error')
plt.legend()
plt.tight_layout()

m_e, c_e = np.polyfit(np.log10(1/n), np.log10(err_e), 1)
m_r, c_r = np.polyfit(np.log10(1/n), np.log10(err_r), 1)

print(f"The slope of the Euler Error is {m_e}")
print(f"The slope of the RK4 Error is {m_r}")

#plt.show()


#part b

n_per_orbit, n_orbits = 500, 100

te2, ye2 = ig.integrate(ig.euler_step, tb.kepler_derivs, y0, (0.0, n_orbits * P), P / n_per_orbit, store_every=10) #integrate n orbits, with a dt of P/n_per_orbit
tr2, yr2 = ig.integrate(ig.rk4_step, tb.kepler_derivs, y0, (0.0, n_orbits * P), P / n_per_orbit, store_every=10)
drift_e2 = tb.relative_drift(tb.specific_energy(ye2))
drift_r2 = tb.relative_drift(tb.specific_energy(yr2))
plt.figure(3)#figsize=(7, 4))
plt.yscale('log')
#plt.ylim(1e-8, 1)
drift_r2 = np.abs(drift_r2)#negative without abs
plt.plot(te2 / P, drift_e2, label = 'Euler')
plt.plot(tr2 / P, drift_r2, label = 'RK4') 
plt.xlabel('time (orbits)')
plt.ylabel(r'$\Delta E / |E_0|$')
plt.title(r'$\Delta E / |E_0|$ vs Orbit Number for both Euler and RK4')
plt.legend()
plt.tight_layout()


rk4_drift = np.mean(drift_r2)

print(f"The average drift per orbit for RK4 is {rk4_drift}")

plt.show()