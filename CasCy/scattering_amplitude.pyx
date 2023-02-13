"""Exponentially scaled plane-wave scattering amplitudes at zero-frequency for TM-TM polarization at zero frequency.
"""
from libc.math cimport exp, sqrt, sinh, cosh, asinh
from cython import cdivision, boundscheck, embedsignature

@cdivision(True)
cdef double expo_diff(int m, double kR, double u, double expo):
    # helper function
    return exp(2*sqrt(m**2 + kR**2) - 2*m*asinh(m/kR) + m*u - expo)


cdef double T_term(int n, double kR, double u, double expo, double pwrc):
    # helper function
    # n>0
    cdef double ed = expo_diff(n, kR, u, expo)
    return (1+exp(-2*n*u))*pwrc*ed

@cdivision(True)
@boundscheck(False)
@embedsignature(True)
cpdef scattering_amplitude(double kR, double u, int m_max, double[:] pwrc):
    r"""
    Exponentially scaled plane-wave scattering amplitudes at zero-frequency for TM-TM polarization at zero frequency.
    The scattering amplitudes are defined as

    .. math::
        T = a_0 + 2 \sum_{m=0}^\infty a_m \cos(m\Theta)

    with :math:`\Theta` being the projection of the scattering angle into the plane perpendicular to the cylinder axis
    and :math:`a_m` the partial-wave reflection coefficients. The parameter `u` is related to :math:`\Theta` by
    the relation :math:`\cosh(u)=-\cos(\Theta)`. The exponentially scaled plane-wave scattering amplitude is then
    defined as
    
    .. math::
        i T\exp(-2 k R \cosh(u/2))\,.

    Parameters
    ----------
    kR : double
        Product of wave vector component along the cylinder axis and the cylinder radius.
    u : double
        Parameter related to the projection of the scattering angle perpendicular to the cylinder axis.
    m_max : integer
        Maximum order of index :math:`m`.
    pwrc : list
        List of partial-wave reflection coefficients of length `m_max + 1`.

    Returns
    -------
    T : double
        Exponentially scaled plane-wave scattering amplitude.

    """
    cdef double T, Tt
    cdef int m

    cdef double expo = 2 * kR * cosh(u / 2)
    cdef int m_init = min(int(kR * sinh(u / 2)), m_max)

    if m_init == 0:
        T = pwrc[0] * expo_diff(0, kR, u, expo)
    else:
        T = T_term(m_init, kR, u, expo, pwrc[m_init])
        T += pwrc[0] * expo_diff(0, kR, u, expo)

    if T == 0.:
        # prevent zero-division error, this happens typically when n_init >> nmax
        return 0.

    # upward summation
    m = m_init + 1
    while (m < m_max + 1):
        Tt = T_term(m, kR, u, expo, pwrc[m])
        T += Tt
        if abs(Tt / T) < 1.e-16:
            break
        m += 1

    # downward summation
    m = m_init - 1
    while (m > 0):
        Tt = T_term(m, kR, u, expo, pwrc[m])
        T += Tt
        if abs(Tt / T) < 1.e-16:
            return T
        m -= 1

    return T