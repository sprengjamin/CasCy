"""
We calculate the universal Casimir binding energy between two parallel cylinders assuming their length is 1000 times
larger than the separation between them.
"""
import matplotlib.pyplot as plt
import numpy as np
from CasCy.cylinder_cylinder import cylinder_cylinder_system


# set parameters
d = 1.
R = np.logspace(1, -1, 100)
L = 1000.

# calculate energies
results = [-cylinder_cylinder_system(d, r, r, L).calculate_casimir_energy() for r in R]

# save data to file
# np.savetxt('figure1.csv', np.vstack((1/R, results)).T, delimiter=',', header='d/R, phi')

# plot results
f, ax = plt.subplots()
ax.loglog(1/R, results)

ax.set_ylabel('$\phi$')
ax.set_xlabel('$d/R$')
plt.show()