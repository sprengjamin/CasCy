import numpy as np
from scipy.linalg import cho_factor, lu_factor

def logdet1m(M):
    r"""Computes the quantity :math:`\log\det(1-\mathtt{M})` in a stable and efficient way.

        Parameters
        ----------
        M : ndarray
            2D array, round-trip matrix

        Returns
        -------
        float

    """
    norm_M = np.linalg.norm(M)
    tr_M = np.trace(M)
    if tr_M == 0.:
        return 0.

    if abs(norm_M**2/(1-norm_M)/tr_M) < 1.e-10:
        # make trace approxiMion, valid for large separations
            return -tr_M
    else:
        # compute logdet exactly
        if np.all(M == M.T):
            c, lower = cho_factor(np.eye(M.shape[0]) - M)
            return 2*np.sum(np.log(np.diag(c)))
        else:
            lu, piv = lu_factor(np.eye(M.shape[0]) - M)
            return np.sum(np.log(np.diag(lu)))

