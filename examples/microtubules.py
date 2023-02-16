"""
We calculate the Casimir binding energy in units of kB*T for two parallel microtubules with radius 12nm and length 50um
for separations from 15 to 60 nm.
"""
import matplotlib.pyplot as plt
import numpy as np
from CasCy.cylinder_cylinder import cylinder_cylinder_system


# set parameters
D = np.linspace(17, 60, 100)*1.e-9 # in m
R = 12.e-9    # in m
L = 50.e-6    # in m

# calculate energies
results = [-cylinder_cylinder_system(d, R, R, L).calculate_casimir_energy() for d in D]

# save data to file
# np.savetxt('figure3.csv', np.vstack((D, results)).T, delimiter=',', header='d (m), DeltaF/(kB T)')

# plot results
f, ax = plt.subplots()
ax.plot(D*1.e9, results)

ax.set_ylabel('Casimir binding energy $\Delta\mathcal{F}/k_B T$')
ax.set_xlabel('$d$ (nm)')
plt.show()