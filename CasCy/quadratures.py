"""Quadrature rules for numerical integration.
"""
import numpy as np

def fcqs_semiinfinite(N):
    r"""Fourier-Chebyshev quadrature rule for integration over the semi-infinite interval (0,oo).

    Parameters
    ----------
    N : int
        quadrature order

    Returns
    -------
    points : numpy.ndarray
        quadrature points
    weights : numpy.ndarray
        quadrature weights

    """
    i = np.arange(1, N + 1, 1)
    t = np.pi / (N + 1) * i
    pts = 1. / (np.tan(t / 2)) ** 2
    wts = np.zeros(N)
    for j in i:
        wts += np.sin(j * t) * (1 - np.cos(j * np.pi)) / j
    wts *= 2 * np.sin(t) * (2 / (N + 1)) / (1 - np.cos(t)) ** 2
    return pts, wts


def fcqs_combined(N):
    r"""
    Quadrature rule for integration over the infinite interval (-oo, oo).

    This method applys `fcqs_semiinfinite` of order `N` to the intervals  (-oo, 0) and (0, oo) separately and combines
    the results.

    Parameters
    ----------
    Parameters
    ----------
    N : int
        quadrature order

    Returns
    -------
    points : numpy.ndarray
        quadrature points
    weights : numpy.ndarray
        quadrature weights

    Returns
    -------

    """
    k, w = fcqs_semiinfinite(N)
    return np.hstack((-k,k[::-1])), np.hstack((w,w[::-1]))