import numpy as np
import matplotlib.pyplot as plt
import astrotools.nbody.drag as d

'''
test_1 = d.initial_state(1.0, 0.1, 0.41)
cap, time, pos = d.encounter(n.rk4_step, d.hill_drag_derivs, test_1, 0.05)
print(cap)
test_2 = d.initial_state(2.0, 0.1, 0.41)
cap, time, pos = d.encounter(n.rk4_step, d.hill_drag_derivs, test_2, 0.05)
print(cap)

#The above works, but it requires that I set st and dv in the hill_drag_derivs function to be None, and then if None default to 0.1 and 0.41.
'''

#a 

dv = 0.0
st = 10**6
r_cap = 0.01
dx0 = 0.02
x0 = np.arange(-3, 3, dx0)

flag = []
times = []
positions = []

for x in x0:
    y_0 = d.initial_state(x, st, dv)
    cap, time, pos = d.encounter(d.rk4_step, d.hill_drag_derivs, y_0, r_cap, st, dv) #encounter modified to now take st and dv for use in stepper and hill_derivatives
    flag.append(cap)
    times.append(time)
    positions.append(pos)

#captured, times, and positions are now lists of booleans, array of time, array of positions for each particle at starting position x0

#find range of x0 where particles are captured. To do this, I'll find the first value where it is captured and the last. might need to check if there are any falses inside the boundary.

captured = [i for i,x in enumerate(flag) if x == True] #all indexes where pebble is captured
x_cap = [x0[i] for i in range(len(flag)) if i in set(captured)] #this shows all x0 values where the pebbles are captured.
xmax = np.max(x_cap)
xmin = np.min(x_cap)
Mdot = 3 * (3 * st)**(2/3) #accretion rate in terms of sigma_peb-omega_k-r_H^2
print(f'The range of captured pebbles is from {xmin} to {xmax} r$_H$. The shear accretion rate Mdot = {Mdot}.')

fig, ax = plt.subplots(figsize = (10, 7))
for x in x0:
    y0 = d.initial_state(x, st, dv)
    cap, time, pos = d.encounter(d.rk4_step, d.hill_drag_derivs, y0, r_cap, st, dv) #encounter modified to now take st and dv for use in stepper and hill_derivatives
    if cap == True:
        ax.plot(pos[:,1], pos[:,0], color = 'red')
    else:
        ax.plot(pos[:,1], pos[:,0], color = 'grey')


ax.plot([],[], color = 'red', label = 'Deflected')
ax.plot([],[], color = 'grey', label = 'Impacted')
ax.set_xlim(-6, 6)
ax.set_ylabel(r'X-Position [r$_H$]')
ax.set_xlabel(r'Y-Position [r$_H$]')
ax.set_ylim(-6, 6)
ax.set_title(f'Paths of pebbles with Stokes number = {st} and r$_c$$_a$$_p$ = {r_cap}')
ax.set_aspect('equal')
legend = ax.legend()
#ax.grid(True, linestyle=':', alpha=0.2, color='white')

plt.tight_layout()
plt.show()