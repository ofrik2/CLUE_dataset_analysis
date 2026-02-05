from __future__ import annotations
import numpy as np
from typing import Tuple

from xlmhg import xlmhg_test


def run_xlmhg(v: np.ndarray, X: int, L: int) -> Tuple[float, float]:
    """
    Runs XL-mHG on a binary hit vector.

    Parameters
    ----------
    v : np.ndarray
        1D binary array (1 = pathway hit, 0 = no hit), ordered by decreasing relevance.
    X : int
        Minimum number of hits required at a cutoff.
    L : int
        Maximum cutoff position (top L items are considered).

    Returns
    -------
    (statistic, p_value) : Tuple[float, float]
        statistic : XL-mHG test statistic
        p_value   : XL-mHG p-value
    """

    v = np.asarray(v, dtype=int).ravel()
    if v.ndim != 1:
        raise ValueError(f"v must be 1D, got shape {v.shape}")

    N = v.size
    if N == 0:
        raise ValueError("v is empty")

    # Clamp L to N (important when pos/neg lists differ in size)
    L_eff = min(int(L), int(N))
    if L_eff < 1:
        raise ValueError(f"L must be >= 1, got L={L}")

    K = int(v.sum())
    if X < 0:
        raise ValueError(f"X must be >= 0, got X={X}")

    # If total hits < X, no cutoff can satisfy the constraint
    if K < X:
        return 0.0, 1.0

    stat, _cutoff, pval = xlmhg_test(v, X, L_eff)
    return float(stat), float(pval)
