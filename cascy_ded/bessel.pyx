r"""
Exponentially scaled Bessel functions of first and second order.

The exponentially scaled Bessel functions of first and second order for order :math:`n` and argument :math:`z` are
denoted as :math:`\tilde{I}_{n}(z)` and :math:`\tilde{K}_{n}(z)`, respectively. We define those functions as

.. math::
        \tilde{I}_{n}(z) = e^{-\eta(n, z)} I_n(z)

        \tilde{K}_{n}(z) = e^{\eta(n, z)} K_n(z)


where :math:`\eta(n, z)=\sqrt{n^2 + z^2} - n \asinh(n/z)` and :math:`I_n` and :math:`K_n` are the modified Bessel
function of the first and second kind, respectively. This definition of the exponentially scaled versions of the
modified Bessel functions is motivated by the uniform asymptotic expansions of the modified Bessel functions for large
order. In this way, the functions are numerically well-behaved even if the order or the argument becomes very large.
"""
import numpy as np
from libc.math cimport exp, sqrt, log1p, asinh
from scipy.special.cython_special cimport k0e, k1e, i0e, i1e
from cython import cdivision

@cdivision(True)
cdef double fraction(double nu, double x):
    r"""Returns the fraction

    .. math::
        \frac{I_\nu(x)}{I_{\nu+1}(x)}

    where :math:`I_\nu` is the modified Bessel function of the first kind of order :math:`\nu`.

    Parameters
    ----------
    nu : float
        Order.
    x : float
        Argument.

    Returns
    -------
    frac : float
        Fraction of modified Bessel functions.
    """
    cdef double an
    cdef double invx = 1/x

    cdef double a1 = 2*(nu+1)*invx
    cdef double a2 = 2*(nu+2)*invx

    cdef double num = a2+1/a1
    cdef double denom = a2
    cdef double ratio = a1*num/denom
    cdef double ratio_last = 0.
    cdef int l = 3
    while True:
        an = (2*nu+2*l)*invx
        num = an+1/num
        denom = an+1/denom
        ratio *= num/denom
        if(ratio == ratio_last):
            return ratio
        ratio_last = ratio
        l += 1

@cdivision(True)
cdef double delta(int n, double z, int j):
    r"""
    Helper function to calculate the exponentially scaled modified Bessel functions.

    Parameters
    ----------
    n : order
        Index.
    z : float
        Argument.
    j : integer
        Difference of orders.

    Returns
    -------
    delta : float

    """
    cdef double s1 = sqrt(n**2 + z**2)
    cdef double s2 = sqrt((n+j)**2 + z**2)
    cdef double d = (2*n*j + j*j)/(s1 + s2)
    cdef double e = d + n*log1p(-(d+j)/(n+j+s2)) - j*asinh((n+j)/z)
    return exp(e)

@cdivision(True)
cdef double Knp1(double Knm1, double Kn, int n, double z):
    r"""
    Implementation of upward recurrence for exponentially scaled modified Bessel function of the second kind.
    Uses :math:`\tilde{K}_{n-1}(z)` and :math:`\tilde{K}_n(z)` to compute :math:`\tilde{K}_{n+1}(z)`.

    Parameters
    ----------
    Knm1 : float
        Exponentially scaled Bessel function :math:`K_{n-1}(z)`.
    Kn : float
        Exponentially scaled Bessel function :math:`K_{n}(z)`.
    n : integer
        Order
    z : float
        Argument

    Returns
    -------
    Knp1 : float
        Exponentially scaled Bessel function :math:`K_{n+1}(z)`.


    """
    return 2*n/z*Kn*delta(n, z, 1) + Knm1*delta(n-1, z, 2)


cpdef bessel_array(int nmax, double z):
    r"""
    Calculates arrays for the exponentially scaled Bessel functions of first and second kind and their derivatives.
    The arrays contain the Bessel functions for orders from :math:`n=0` to `nmax`.
    The function returns the four arrays as a tuple with the ordering corresponding to
    (:math:`\tilde{I}_{n}(z)`, :math:`\tilde{I}_{n}'(z)`, :math:`\tilde{K}_{n}(z)`, :math:`\tilde{K}_{n}'(z)`).

    Parameters
    ----------
    nmax : integer
        maximum order
    z : float
        argument

    Returns
    -------
    (nd.array, nd.array, nd.array, nd.array)


    """
    cdef double[:] Kn = np.empty(nmax+3)
    cdef double[:] In = np.empty(nmax+2)
    cdef double[:] Knprime = np.empty(nmax+1)
    cdef double[:] Inprime = np.empty(nmax+1)
    Kn[0] = k0e(z)
    In[0] = i0e(z)
    cdef double delta_0z1 = delta(0,z,1)
    Kn[1] = k1e(z)*delta_0z1
    In[1] = i1e(z)/delta_0z1
    Knprime[0] = -Kn[1]/delta_0z1
    Inprime[0] = In[1]*delta_0z1
    for n in range(1, nmax+1):
        Kn[n+1] = Knp1(Kn[n-1], Kn[n], n, z)
        delta_nz1 = delta(n,z,1)
        delta_nm1z1 = delta(n-1,z,1)
        Knprime[n] = -0.5*(Kn[n-1]*delta_nm1z1 + Kn[n+1]/delta_nz1)
    Kn[nmax+2] = Knp1(Kn[nmax], Kn[nmax+1], nmax+1, z)
    for n in range(1, nmax+1):
        In[n+1] = 1/z/(Kn[n+2]/delta(n+1,z,1) + Kn[n+1]/fraction(n+1, z))
        Inprime[n] = 0.5*(In[n-1]/delta(n-1,z,1) + In[n+1]*delta(n,z,1))
    return np.asarray(In[:nmax+1]), np.asarray(Inprime), np.asarray(Kn[:nmax+1]), np.asarray(Knprime)
