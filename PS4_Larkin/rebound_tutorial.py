import sys
import rebound

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({"figure.dpi": 110, "font.size": 10})

sim = rebound.Simulation()
sim.units = ("yr", "AU", "Msun")
#print(f"G = {sim.G:.4f} (4 pi^2 = {4 * np.pi**2:.4f})")

sim.add(m=1.0) # Sun
sim.add(m=9.54e-4, a=5.20, e=0.048) # Jupiter
sim.add(m=2.86e-4, a=9.58, e=0.054, f=2.0) # Saturn, starting 2 rad ahead in true anomaly

sim.move_to_com()
sim.integrator = "ias15"

''' This prints a, e, P for each planet
for i, p in enumerate(sim.particles[1:], start=1):
    o = p.orbit(primary=sim.particles[0])
    print(f"particle {i}: a = {o.a:.3f} AU, e = {o.e:.3f}, P = {o.P:.2f} yr")
'''

times = np.linspace(0, 60, 3000) # yr
xy = np.empty((len(times), 3, 2))
E0 = sim.energy()

for k, t in enumerate(times):
    sim.integrate(t)
    for i, p in enumerate(sim.particles):
        xy[k, i] = p.x, p.y

#print(f"t = {sim.t:.1f} yr, relative energy error = {abs(sim.energy() - E0) / abs(E0):.1e}")

fig, ax = plt.subplots(figsize=(5, 5))
for i, (name, c) in enumerate([("Sun", "k"), ("Jupiter", "C1"), ("Saturn", "C0")]):
    ax.plot(xy[:, i, 0], xy[:, i, 1], color=c, lw=0.8, label=name)

ax.set(aspect="equal", xlabel="x (AU)", ylabel="y (AU)", title="60 years, center-of-mass frame")
ax.legend(loc="upper right")

def angle_from_star(xy, i, j):
    """Angle of body j relative to body i as seen from body 0, in degrees."""
    d_i = xy[:, i] - xy[:, 0]
    d_j = xy[:, j] - xy[:, 0]
    ang = np.arctan2(d_j[:, 1], d_j[:, 0]) - np.arctan2(d_i[:, 1], d_i[:, 0])
    return np.degrees((ang + np.pi) % (2 * np.pi) - np.pi)

phi = angle_from_star(xy, 1, 2)

# Synodic period: spacing of successive conjunctions, where phi falls through zero

cross = np.where((phi[:-1] > 0) & (phi[1:] <= 0))[0]
P_J, P_S = (sim.particles[i].orbit(primary=sim.particles[0]).P for i in (1, 2))

print(f"measured synodic period {np.mean(np.diff(times[cross])):.1f} yr; " f"predicted {1 / (1 / P_J - 1 / P_S):.1f} yr")

fig1, ax1 = plt.subplots(figsize=(8, 2.8))
ax1.plot(times, phi, lw=0.9)
ax1.set(xlabel="time (yr)", ylabel="Saturn - Jupiter (deg)", ylim=(-180, 180))
#plt.show()

sim2 = rebound.Simulation()
sim2.units = ("yr", "AU", "Msun")
sim2.add(m=1.0)
sim2.add(m=9.54e-4, a=5.20)
r = 3.0
sim2.add(m=0, x=r, y=0, vx=0, vy=np.sqrt(sim2.G / r))
sim2.move_to_com()
sim2.integrator = "ias15"
o = sim2.particles[2].orbit(primary=sim2.particles[0])
print(f"test particle: a = {o.a:.4f} AU, e = {o.e:.1e}")
sim2.integrate(100)
o = sim2.particles[2].orbit(primary=sim2.particles[0])
4
print(f"after 100 yr: a = {o.a:.4f} AU, e = {o.e:.4f} (Jupiter has excited a small eccentricity)")

plt.show()