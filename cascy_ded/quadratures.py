import numpy as np
from scipy.integrate import quad
from numba import njit

@njit
def fcqs_infinite(N):
    i = np.arange(1, N+1)
    ti = np.pi * i / (N + 1)
    x = 1 / np.tan(ti)
    w = np.pi / (np.sin(ti) ** 2 * (N + 1))
    return x, w

@njit
def weights(N):
    r"""Fourier-Chebyshev quadrature weights.

    """
    i = np.arange(1, N + 1, 1)
    t = np.pi / (N + 1) * i
    wts = np.zeros(N)
    for j in i:
        wts += np.sin(j * t) * (1 - np.cos(j * np.pi)) / j
    wts *= 2 * np.sin(t) * (2 / (N + 1)) / (1 - np.cos(t)) ** 2
    return wts


@njit
def fcqs_semiinfinite(N):
    r"""Fourier-Chebyshev quadrature rule.

    Parameters
    ----------
    N : int
        quadrature order

    Returns
    -------
    points : nd.array
        quadrature points
    weights : nd.array
        quadrature weights

    """
    i = np.arange(1, N + 1, 1)
    t = np.pi / (N + 1) * i
    pts = 1. / (np.tan(t / 2)) ** 2
    wts = weights(N)
    return pts, wts


@njit
def fcqs_combined(N):
    k, w = fcqs_semiinfinite(N)
    return np.hstack((-k,k[::-1])), np.hstack((w,w[::-1]))