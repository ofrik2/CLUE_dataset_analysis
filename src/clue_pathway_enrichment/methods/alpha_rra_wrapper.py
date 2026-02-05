from __future__ import annotations
import numpy as np
from typing import Tuple

# Use your existing implementation:
# file:/mnt/data/alpha_rra.py
from alpha_rra import alpha_rra_test  # adjust import if you place it under src/

def run_alpha_rra(v: np.ndarray, alpha: float, n_perm: int, seed: int | None) -> Tuple[float, float]:
    """
    Returns (rho, p_value)
    """
    rho, p = alpha_rra_test(v, alpha=alpha, n_permutations=n_perm, random_state=seed)
    return float(rho), float(p)
