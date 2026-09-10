#question 2

import numpy as np
from astrotools import constants as c
from astrotools.cloud import collapse

rho_core = collapse.number_density_to_mass_density(c.N_H2_CORE) #makes density of ref core
M_j_ref = collapse.jeans_mass(c.T_CORE, rho_core) #makes mj of ref core
t_ff_ref = collapse.free_fall_time(rho_core) #makes ff time of ref core
L_j_ref = collapse.jeans_length(c.T_CORE, rho_core) #makes Lj of ref core

jeans_m = f"Jean's Mass = {M_j_ref} kg." #makes string
print(jeans_m)

jeans_l = f"Jean's Length = {L_j_ref} m." #makes string
print(jeans_l)

t_ff = f"Free-Fall Time = {t_ff_ref} s."
print(t_ff)