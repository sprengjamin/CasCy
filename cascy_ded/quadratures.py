import numpy as np

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

def fcqs_combined(N):
    k, w = fcqs_semiinfinite(N)
    return np.hstack((-k,k[::-1])), np.hstack((w,w[::-1]))