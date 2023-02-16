import numpy as np
from .quadratures import fcqs_combined, fcqs_semiinfinite
from .cylinder_reflection import pwrc_TMTM, reflection_matrix
from .matrix_operations import logdet1m

class plane_cylinder_system:
    r"""
    A class to represent the plane-cylinder geometry.

    Attributes
    ----------
    d : float
        separation between plane and cylinder
    R : float
        cylinder radius
    L : float
        cylinder length
    x_quad : function
        Quadrature scheme for x integration over the interval (-oo, oo).
        The function has the signature `function(int)->(list, list)`, where the return values are the nodes and weights
        of the quadrature scheme, respectively. The default value is`fcqs_combined` (see module `quadratures.py`).
    z_quad : function
        Quadrature scheme for z integration over the interval (0, oo).
        The function has the signature `function(int)->(list, list)`, where the return values are the nodes and weights
        of the quadrature scheme, respectively. The default value is`fcqs_semiinfinite` (see module `quadratures.py`).

    Methods
    -------
    calculate_casimir_energy(eta_Nx=2., Nx=None, Nz=20, eta_mmax=10., mmax=None) -> float
        Calculates the Casimir energy in units of :math:`k_B T` for the defined geometry.
    """
    def __init__(self, d, R, L, x_quad = fcqs_combined, z_quad = fcqs_semiinfinite):
        """
        Constructs all the necessary attributes of the plane-cylinder object.

        Parameters
        ----------
        d : float
            separation between plane and cylinder
        R : float
            cylinder radius
        L : float
            cylinder length
        x_quad : function
            Quadrature scheme for x integration over the interval (-oo, oo).
            The function has the signature `function(int)->(list, list)`, where the return values are the nodes and weights
            of the quadrature scheme, respectively. The default value is`fcqs_combined` (see module `quadratures.py`).
        z_quad : function
            Quadrature scheme for z integration over the interval (0, oo).
            The function has the signature `function(int)->(list, list)`, where the return values are the nodes and weights
            of the quadrature scheme, respectively. The default value is`fcqs_semiinfinite` (see module `quadratures.py`).
        """
        self.d = d
        self.R = R
        self.L = L
        self.x_quad = x_quad
        self.z_quad = z_quad

    def calculate_casimir_energy(self, eta_Nx=2., Nx=None, Nz=20, eta_mmax=10., mmax=None):
        r"""
        Calculates the Casimir energy in units of :math:`k_B T` for the defined geometry.

        Parameters
        ----------
        eta_Nx : float
            Convergence parameter for x integration used to determine the discretization order `Nx` by means of the
            scaling law :math:`N_x = \eta_{N_x} \sqrt{R/L}` which becomes valid when :math:`R \gg L`. The scaling law
            assumes that the quadrature scheme for x integration has not been changed from the default function.
            The larger the value of `eta_Nx`, the more accurate the result. Default value is set to `2.` and corresponds
            to a numerical error of about 1%.
        Nx : int
            Discretization order of the quadrature scheme for the x integration. The larger the value, the more accurate
            the result. Default is `None`. Setting a value here overwrites the value determined by `eta_Nx`.
        Nz : int
            Discretization order of the quadrature scheme for the z integration. The larger the value, the more accurate
            the result. Default is `20` and corresponds to a numerical error of about `10^-4`.
        eta_mmax : float
            Convergence parameter to determine `mmax` by means of the scaling law
            :math:`m_\text{max} = \eta_{m_\text{max}} R/L` which becomes valid when :math:`R \gg L`.
            The larger the value of `eta_mmax`, the more accurate the result. Default value is set to `10.` and
            corresponds to a numerical error of about `10^-8`.
        mmax : int
            Maximum value of the cylindrical multipole index `m` included in the calculation. The larger the value, the
            more accurate the result. Default is `None`. Setting a value here overwrites the value determined by
            `eta_mmax`.

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


