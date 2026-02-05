# alpha_rra.py
#
# Implementation of the alpha-RRA (α-RRA) enrichment test
# with a permutation-based empirical p-value, inspired by
# the MAGeCK approach, but operating on a single binary vector.
#
# Python 3.8 compatible.

import numpy as np
from scipy.stats import beta as beta_dist


def binary_vector_to_ranks(v):
    """
    Convert a binary vector v (length N) into 1-based ranks of hits.

    Parameters
    ----------
    v : 1D array-like of {0,1}
        Binary vector, where 1 indicates a "hit".

    Returns
    -------
    ranks : 1D numpy array of ints
        Positions (1..N) where v == 1.
    """
    v = np.asarray(v, dtype=int)
    indices = np.where(v == 1)[0]
    return indices + 1  # ranks are 1..N


def alpha_rra_rho_from_ranks(ranks, M, alpha=0.05):
    """
    Compute the alpha-RRA statistic (rho_alpha) for one gene/feature
    given its hit ranks.

    Parameters
    ----------
    ranks : 1D array-like of ints
        1-based ranks (positions 1..M) where this feature has hits.
    M : int
        Total length of the ranked list (N).
    alpha : float
        Truncation parameter (0 < alpha <= 1).
        Only hits with u_i = ranks_i / M <= alpha are considered.

    Returns
    -------
    rho_alpha : float
        Alpha-RRA score in [0,1]. Smaller values indicate stronger
        evidence for enrichment near the top of the ranking.
    """
    ranks = np.asarray(ranks, dtype=float)
    n = ranks.size
    if n == 0:
        return 1.0

    # Convert to percentiles in (0,1]
    u = ranks / float(M)
    u.sort()

    # Keep only hits in the top alpha fraction
    mask = (u <= alpha)
    u_trunc = u[mask]
    j = u_trunc.size

    if j == 0:
        # No hits in the top alpha fraction: no evidence for enrichment
        return 1.0

    # n = total number of hits for this vector (before truncation)
    # Under the RRA null, the k-th order statistic of n Uniform(0,1)
    # follows a Beta(k, n+1-k) distribution.
    pvals = []
    for k in range(1, j + 1):
        x = u_trunc[k - 1]
        p_k = beta_dist.cdf(x, k, n + 1 - k)
        pvals.append(p_k)

    rho_alpha = float(np.min(pvals))
    return rho_alpha


def alpha_rra_test(v, alpha=0.05, n_permutations=1000, random_state=None):
    """
    Full alpha-RRA test on a single binary vector v.

    This function:
      1. Converts v into ranks of hits.
      2. Computes the alpha-RRA statistic rho_alpha.
      3. Uses permutation of hit positions (keeping N and K fixed)
         to obtain an empirical p-value.

    Parameters
    ----------
    v : 1D array-like of {0,1}
        Binary vector of length N, with K ones.
    alpha : float, default 0.05
        Truncation parameter: only hits with percentile <= alpha
        are used in the RRA statistic.
    n_permutations : int, default 1000
        Number of permutations for empirical p-value.
    random_state : int or numpy.random.Generator or None
        Seed or RNG for reproducibility. If None, a new Generator is used.

    Returns
    -------
    rho_alpha : float
        Observed alpha-RRA statistic.
    p_emp : float
        Empirical p-value:
            p_emp = (count_null_rho <= rho_obs + 1) / (n_permutations + 1)

        NOTE: This does NOT include the global floor 1/(G * n_permutations)
        used in MAGeCK for many genes; that should be applied by the caller
        if desired, once G is known.
    """
    v = np.asarray(v, dtype=int)
    N = v.size

    # Prepare RNG
    if isinstance(random_state, np.random.Generator):
        rng = random_state
    else:
        rng = np.random.default_rng(random_state)

    # Observed statistic
    ranks_obs = binary_vector_to_ranks(v)
    rho_obs = alpha_rra_rho_from_ranks(ranks_obs, M=N, alpha=alpha)

    K = ranks_obs.size
    if K == 0 or rho_obs >= 1.0:
        # No hits or no evidence => p-value is 1.0
        return float(rho_obs), 1.0

    # Permutation-based null: reshuffle the hit positions
    count_le = 0
    for _ in range(n_permutations):
        perm_indices = rng.choice(N, size=K, replace=False)
        ranks_perm = perm_indices + 1
        rho_perm = alpha_rra_rho_from_ranks(ranks_perm, M=N, alpha=alpha)
        if rho_perm <= rho_obs:
            count_le += 1

    p_emp = (count_le + 1.0) / (n_permutations + 1.0)
    return float(rho_obs), float(p_emp)
