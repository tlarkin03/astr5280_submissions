import numpy as np
import matplotlib.pyplot as plt
from astrotools import constants as c
from astrotools.disk import structure as s
from astrotools.disk import viscous as v

r0 = 10 #radius of 10 AU
nu = 1e11 #viscosity of 1e11 m^2/s
r = v.make_grid(1*c.AU,1000*c.AU,2000)#grid 1-1000 AU; n = 2000 is 20 cells/decade, 8000 is 80 cells/decade
m = .01*c.M_SUN #total mass of disk
tau_to_t = (10 * c.AU)**2 / (12 * 1e11) #multiply given tau value by this to get t in seconds
t = 0.01 * tau_to_t #at tau = 0.01, t is this in seconds (from working backwards)
sigma_ring = v.ring_solution(r.r,t,m,r0,nu) #generates ring surface mass density 

#a

''' evolve to tau = 0.1 at resolutions from 20 to 80 cells per decade. plot error as fraction of peak Sigma, against cells/decade 
on log axes. Report order of convergence and resolution for 1% accuracy (what cells/decade gives error of 99%)
'''

dt = v.cfl_timestep(r,nu)
t_end = 0.1 * tau_to_t
time, sigma_e, mass = v.evolve(sigma_ring,r,nu,t_end,dt)



#is the way I've done it correct for how to define t?
#where does order of convergence/error come from?

#loop over time and plot sigma_e to see what is happening to the calculation

#b

'''40 cells/decade, run with dt = 1.5 x the value that v.cfl_timestep returns with saftey = 1. How does the run fail, and 
where does failure start in the grid?
'''



#c

'''Continue 40 cell/decade run to tau = 0.5, find time after tau = 0.1 where error exceeds 1%, and where largest error is. change
inner edge to 0.3, repeat. 
'''