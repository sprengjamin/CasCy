"""
Analyze convergence of the Casimir energy with respect to the parameters eta_Nx and Nz.

We consider the geometries of a plane and a cylinder, and of two cylinders with equal radii. The convergence analysis is
performed for varying ratio of cylinder radius and separation.

The full results may take a few minutes to calculate and are saved as png-files.
"""
import matplotlib.pyplot as plt
import numpy as np
from CasCy.plane_cylinder import plane_cylinder_system
from CasCy.cylinder_cylinder import cylinder_cylinder_system

print('Plane-cylinder geometry:')
R = [0.01, 0.1, 1, 10, 100]
labels = ['$R/d = 0.01$', '$0.1$', '$1$','$10$', '$100$']

print('Convergence for eta_Nx')
f, ax = plt.subplots()
ETA = np.arange(1, 7, 1)
for i, r in enumerate(R):
    print('Calculating R/d:', r)
    s = plane_cylinder_system(1., r, 1.)
    Y = np.array([s.calculate_casimir_energy(eta_Nx=eta) for eta in ETA])
    ax.semilogy(ETA[:-1], np.abs(Y[:-1]/Y[-1]-1.), '-', marker='.', label=labels[i])

ax.legend()
ax.set_title('plane-cylinder geometry ')
ax.set_ylabel('relative error')
ax.set_xlabel('$\eta_{N_x}$')
plt.savefig('convergence_plcy_eta_Nx.png')

print()
print('Convergence for Nz')
f, ax = plt.subplots()
NZ = np.arange(5, 65, 5)
for i, r in enumerate(R):
    print('Calculating R/d:', r)
    s = plane_cylinder_system(1., r, 1.)
    Y = np.array([s.calculate_casimir_energy(Nz=Nz) for Nz in NZ])
    ax.semilogy(NZ[:-1], np.abs(Y[:-1] / Y[-1] - 1.), '-', marker='.', label=labels[i])

ax.legend()
ax.set_title('plane-cylinder geometry ')
ax.set_ylabel('relative error')
ax.set_xlabel('$N_z$')
plt.savefig('convergence_plcy_Nz.png')

print()
print('Cylinder-cylinder geometry:')
R = [0.01, 0.1, 1, 10, 100]
labels = ['$R/d = 0.01$', '$0.1$', '$1$','$10$', '$100$']

print('Convergence for eta_Nx')
f, ax = plt.subplots()
ETA = np.arange(2, 14, 2)
for i, r in enumerate(R):
    print('Calculating R/d:', r)
    s = cylinder_cylinder_system(1., r, r, 1.)
    Y = np.array([s.calculate_casimir_energy(eta_Nx=eta) for eta in ETA])
    ax.semilogy(ETA[:-1], np.abs(Y[:-1]/Y[-1]-1.), '-', marker='.', label=labels[i])

ax.legend()
ax.set_title('cylinder-cylinder geometry with equal radii')
ax.set_ylabel('relative error')
ax.set_xlabel('$\eta_{N_x}$')
plt.savefig('convergence_cycy_eta_Nx.png')

print()
print('Convergence for Nz')
f, ax = plt.subplots()
NZ = np.arange(5, 65, 5)
for i, r in enumerate(R):
    print('Calculating R/d:', r)
    s = cylinder_cylinder_system(1., r, r, 1.)
    Y = np.array([s.calculate_casimir_energy(Nz=Nz) for Nz in NZ])
    ax.semilogy(NZ[:-1], np.abs(Y[:-1] / Y[-1] - 1.), '-', marker='.', label=labels[i])

ax.legend()
ax.set_title('cylinder-cylinder geometry ')
ax.set_ylabel('relative error')
ax.set_xlabel('$N_z$')
plt.savefig('convergence_cycy_Nz.png')