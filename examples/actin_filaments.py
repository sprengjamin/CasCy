"""
We calculate the Casimir binding energy in units of kB*T for two parallel actin filaments with radius 3nm and length
15um for separations from 3 to 20 nm.
"""
import matplotlib.pyplot as plt
import numpy as np
from CasCy.cylinder_cylinder import cylinder_cylinder_system


# set parameters
D = np.linspace(3.7, 20, 100)*1.e-9 # in m
R = 3.e-9    # in m
L = 15.e-6    # in m

# calculate energies
results = [-cylinder_cylinder_system(d, R, R, L).calculate_casimir_energy() for d in D]

# save data to file
# np.savetxt('figure2A.csv', np.vstack((D, results)).T, delimiter=',', header='d (nm), ΔF/kB T')

# plot results
f, ax = plt.subplots()
ax.plot(D*1.e9, results)

ax.set_ylabel('Casimir binding energy $\Delta\mathcal{F}/k_B T$')
ax.set_xlabel('$d$ (nm)')
plt.show()