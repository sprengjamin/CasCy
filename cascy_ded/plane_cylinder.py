"""
We calculate the Casimir energy within the scattering formalism utilizing plane waves.
"""
import numpy as np
from .quadratures import fcqs_combined, fcqs_semiinfinite
from .cylinder_reflection import pwrc_TMTM, reflection_matrix
from .matrix_operations import logdet1m

class plane_cylinder_system:
    def __init__(self, d, R, L):
        self.d = d
        self.R = R
        self.L = L

        self.x_quad = fcqs_combined
        self.z_quad = fcqs_semiinfinite

    def calculate_casimir_energy(self, eta_Nx=2., Nx=None, Nz=20, eta_mmax=10., mmax=None):
        r"""
            Casimir energy between a cylinder (radius `R`, length `L`) and a plane at separation `d` in units of :math:`k_B T`.

            Parameters
            ----------
            d : float
                Separation.
            R : float
                Cylinder radius.

            Returns
            -------
            energy : float
                Casimir energy in units of :math:`k_B T`.

        """
        rho = max(self.R  / self.d, 50.)
        if Nx == None:
            Nx = int(eta_Nx * np.sqrt(rho))
        if mmax == None:
            mmax = int(eta_mmax * rho)

        # precompute quadrature nodes and weights
        X1, W1 = self.z_quad(Nz)
        Kz = X1 / self.d
        Wz = W1 / self.d

        X2, W2 = self.x_quad(Nx)
        Kx = X2 / self.d
        Wx = W2 / self.d

        energy_per_length = 0.
        for i, kz in enumerate(Kz):
            # cylider reflection
            pwrc = pwrc_TMTM(mmax, kz * self.R)
            R_cy = reflection_matrix(self.R, kz, Kx, Wx, mmax, pwrc)

            # translation
            kappa = np.sqrt(Kx ** 2 + kz ** 2)
            translation = np.exp(-kappa * self.d)

            # round-trip
            M = np.diag(translation) @ R_cy @ np.diag(-1. * translation)

            # energy contribution
            energy_per_length += Wz[i] / 2 / np.pi * logdet1m(M)
        return self.L*energy_per_length


