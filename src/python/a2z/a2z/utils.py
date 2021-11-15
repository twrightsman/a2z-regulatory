"""
utils module

miscellaneous utility functions
"""

from collections import Counter
import warnings

import numpy as np
import scipy.interpolate


def gc_content(sequence: str) -> float:
    cnt = Counter(sequence)
    return (cnt["G"] + cnt["C"]) / len(sequence)


def sliding_window_starts(start, end, length, stride):
    """ Generator for start of sliding intervals that are fully-enclosed within given interval """
    # skip short intervals
    if end - start >= length:
        yield from range(start, end - length + 1, stride)


def kMedoids(D, k, tmax = 100):
    """ Perform k-medoids clustering given a distance matrix D and k number of medoids

    DOI: 10.13140/2.1.4453.2009
    """
    # determine dimensions of distance matrix D
    m, n = D.shape

    if k > n:
        raise Exception('too many medoids')

    # find a set of valid initial cluster medoid indices since we
    # can't seed different clusters with two points at the same location
    valid_medoid_inds = set(range(n))
    invalid_medoid_inds = set([])
    rs,cs = np.where(D==0)
    # the rows, cols must be shuffled because we will keep the first duplicate below
    index_shuf = list(range(len(rs)))
    np.random.shuffle(index_shuf)
    rs = rs[index_shuf]
    cs = cs[index_shuf]
    for r,c in zip(rs,cs):
        # if there are two points with a distance of 0...
        # keep the first one for cluster init
        if r < c and r not in invalid_medoid_inds:
            invalid_medoid_inds.add(c)
    valid_medoid_inds = list(valid_medoid_inds - invalid_medoid_inds)

    if k > len(valid_medoid_inds):
        raise Exception('too many medoids (after removing {} duplicate points)'.format(
            len(invalid_medoid_inds)))

    # randomly initialize an array of k medoid indices
    M = np.array(valid_medoid_inds)
    np.random.shuffle(M)
    M = np.sort(M[:k])

    # create a copy of the array of medoid indices
    Mnew = np.copy(M)

    # initialize a dictionary to represent clusters
    C = {}
    for t in range(tmax):
        # determine clusters, i. e. arrays of data indices
        J = np.argmin(D[:,M], axis=1)
        for kappa in range(k):
            C[kappa] = np.where(J==kappa)[0]
        # update cluster medoids
        for kappa in range(k):
            J = np.mean(D[np.ix_(C[kappa],C[kappa])],axis=1)
            j = np.argmin(J)
            Mnew[kappa] = C[kappa][j]
        np.sort(Mnew)
        # check for convergence
        if np.array_equal(M, Mnew):
            break
        M = np.copy(Mnew)
    else:
        # final update of cluster memberships
        J = np.argmin(D[:,M], axis=1)
        for kappa in range(k):
            C[kappa] = np.where(J==kappa)[0]

    # return results
    return M, C


def q_values(p_values, robust = False):
    """ Computes q-values from a list of p-values

    Uses method from Remark B, Storey and Tibshirani 2003.
    Implementation also based off of the WGCNA R package.

    Gives the same q-values as WGCNA given the same pi_0 value.
    """
    p_values = np.array(p_values, dtype = np.float)
    sort = np.argsort(p_values)
    m = p_values.shape[0]

    # estimate pi_0 using a natural cubic spline
    lambdas = np.arange(0, 0.90, 0.05)
    pi_0 = min(scipy.interpolate.CubicSpline(
        x = lambdas,
        y = [(p_values >= x).mean() / (1 - x) for x in lambdas],
        bc_type = 'natural'
    )(lambdas.max()), 1)

    def rank(x):
        """ Computes, for each element, the number of list elements less than or equal to it """
        a, inv, cnt = np.unique(x, return_counts = True, return_inverse = True)
        ranks = cnt.cumsum()
        return np.array([ranks[i] for i in inv])

    # compute q-values
    v = rank(p_values)
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', 'divide by zero encountered in true_divide')
        if robust:
            q_values = np.clip((pi_0 * m * p_values) / (v * (1 - ((1 - p_values) ** m))), a_min = 0, a_max = 1)
        else:
            q_values = np.clip((pi_0 * m * p_values) / v, a_min = 0, a_max = 1)

    for i in range(m - 2, -1, -1):
        q_values[sort[i]] = min(
            q_values[sort[i]],
            q_values[sort[i + 1]]
        )

    return q_values
