r"""
Reflection of plane waves on a cylinder.

Consider a half-space defined by an axis perpendicular to the cylinder axis (`x`-axis) and the cylinder axis. Here, we
consider scattering on a cylinder for plane waves incident from this half-space and scattered back into the same
half-space.
"""

import numpy as np
from math import exp, sqrt, cosh, acosh
from .bessel import bessel_array
from .scattering_amplitude import scattering_amplitude

def pwrc_TMTM(m_max, kR):
    r"""
    Calculates an array of coefficients for index :math:`m=0` up to `m_max`. The calculated coefficients are an
    exponentially scaled version of the TM-TM partial-wave reflection coefficient (pwrc) evaluated at zero frequency
    defined as

    .. math::
        e^{-2\eta(n, z) } |\braket{m, \text{TM} | \mathcal{R} | m, \text{TM}}|

    where :math:`\eta(n, z)=\sqrt{n^2 + z^2} - n \asinh(n/z)` and

    .. math::
        \braket{m, \text{TM} | \mathcal{R} | m, \text{TM}} = -\frac{i\pi}{2}(-1)^m \frac{I_m'(k R)}{K_m'(k R)}

    wih :math:`I_n` and :math:`K_n` the modified Bessel function of the first and second kind, respectively. The
    exponential scaling here is inherited from the exponentially scaled versions of the modified Bessel functions.

    Parameters
    ----------
    m_max : integer
        Maximum order of index :math:`m`.
    kR : float
        Product of wave vector component along the cylinder axis and the cylinder radius.

    Returns
    -------
    pwrc : nd.array
        List of exponentially scaled partial-wave reflection coefficients.

    """
    Im, Imp, Km, Kmp = bessel_array(m_max, kR)
    return np.pi/2*Imp/Kmp





def pwrme(R, k, Kxi, Kxs, m_max, pwrc):
    r"""
    Exponentially scaled plane-wave reflection matrix elements (pwrme) at zero frequency for TM-TM polarization
    defined as

    .. math::
        e^{-(\kappa^{(s)} + \kappa^{(i)})R} \braket{K_x^{(s)}, \text{TM} | \mathcal{R} | K_x^{(i)}, \text{TM}}

    with

    .. math::
        \braket{K_x^{(s)}, \text{TM} | \mathcal{R} | K_x^{(i)}, \text{TM}} = \frac{2 T}{i \kappa^{(s)}}\,,

    where :math:`T` is the plane-wave scattering amplitude and :math:`\kappa^{(i/s)} = \sqrt{(K_x^{(i/s)})^2 + k^2}`.

    Parameters
    ----------
    R : float
        Cylinder radius.
    k : float
        Wave vector component along the cylinder axis.
    Kxi : float
        Incident wave-vector component perpendicular to the cylinder axis :math:`K_x^{(i)}`.
    Kxs : float
        Scattered wave-vector component perpendicular to the cylinder axis :math:`K_x^{(s)}`.
    m_max : integer
        Maximum order of index :math:`m`.
    pwrc : list
        List of partial-wave reflection coefficients of length `m_max + 1`.

    Returns
    -------
    pwrme : float

    """
    kappai = sqrt(Kxi**2 + k**2)
    kappas = sqrt(Kxs**2 + k**2)
    mcostheta = (Kxi*Kxs+kappai*kappas)/k**2
    if mcostheta <= 1.:
        u = 0.
    else:
        u = acosh(mcostheta)
    exponent = 2*k*R*cosh(u/2) - (kappai+kappas)*R
    if exponent < -37.:
        return 0.
    else:
        e = exp(exponent)
        T = scattering_amplitude(k*R, u, m_max, pwrc)
        return 2/kappas*T*e


def reflection_matrix(R, k, Kx, Wx, m_max, pwrc):
    r"""
    Discretized plane-wave reflection matrix at cylinder. The matrix elements are given by

    .. math::
        M_{ij} = \frac{\sqrt{\mathtt{Wx[i]} \mathtt{Wx[j]}}}{2\pi}\mathtt{pwrme(R, k, Kx[j], Kx[i], m_max, pwrc)}\,.

    Parameters
    ----------
    R : float
        Cylinder radius.
    k : float
        Wave vector component along the cylinder axis.
    Kx, Wx : arrays
        Quadrature nodes and weights of the infinite interval (-inf, inf).
    m_max : integer
        Maximum order of index :math:`m`.
    pwrc : list
        List of partial-wave reflection coefficients of length `m_max + 1`.

    Returns
    -------

    """
    Nx = len(Kx)
    assert (Nx % 2 == 0)
    mat = np.empty((Nx, Nx))
    for i in range(Nx // 2):  # row/s index
        for j in range(i, Nx - i):  # column/i index
            me = pwrme(R, k, Kx[j], Kx[i], m_max, pwrc)
            weight = sqrt(Wx[i] * Wx[j]) / 2 / np.pi
            mat[i, j] = me * weight

            # symmetry w.r.t. center
            mat[Nx - 1 - i, Nx - 1 - j] = mat[i, j]

    # symmetry w.r.t. diagonal
    for i in range(Nx // 2 - 1):  # row/s index
        for j in range(i + 1, Nx - i - 1):  # column/i index
            kappa_s = sqrt(Kx[i] ** 2 + k ** 2)
            kappa_i = sqrt(Kx[j] ** 2 + k ** 2)
            mat[j, i] = kappa_s / kappa_i * mat[i, j]

            mat[Nx - 1 - j, Nx - 1 - i] = kappa_s / kappa_i * mat[Nx - 1 - i, Nx - 1 - j]
    return mat