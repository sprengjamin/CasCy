"""
We calculate the Casimir energy within the scattering formalism utilizing plane waves.
"""
import numpy as np
from .quadratures import fcqs_combined, fcqs_semiinfinite
from .cylinder_reflection import pwrc_TMTM, reflection_matrix
from .matrix_operations import logdet1m

class cylinder_cylinder_system:
    def __init__(self, d, R1, R2, L):
        self.d = d
        self.R1 = R1
        self.R2 = R2
        self.L = L

        self.x_quad = fcqs_combined
        self.z_quad = fcqs_semiinfinite

    def calculate_casimir_energy(self, eta_Nx=4., Nx=None, Nz=20, eta_mmax=10., mmax=None):
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
        rho1 = max(self.R1 / self.L, 50.)
        rho2 = max(self.R2 / self.L, 50.)
        rho = max(rho1, rho2)
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

        result = 0.
        for i, kz in enumerate(Kz):
            # cylider reflection
            pwrc1 = pwrc_TMTM(mmax, kz * self.R1)
            R_cy1 = reflection_matrix(self.R1, kz, Kx, Wx, mmax, pwrc1)
            if self.R1 == self.R2:
                R_cy2 = R_cy1
            else:
                pwrc2 = pwrc_TMTM(mmax, kz * self.R2)
                R_cy2 = reflection_matrix(self.R2, kz, Kx, Wx, mmax, pwrc2)

            # translation
            kappa = np.sqrt(Kx ** 2 + kz ** 2)
            half_translation = np.diag(np.exp(-0.5*kappa * self.d))

            # round-trip
            M1 = half_translation @ R_cy1 @ half_translation
            M2 = half_translation @ R_cy2 @ half_translation
            M = M1 @ M2

            # energy contribution
            result += Wz[i] / 2 / np.pi * logdet1m(M)
        return self.d/self.L*result


