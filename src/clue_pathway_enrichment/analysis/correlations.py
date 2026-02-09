from __future__ import annotations

from typing import Mapping

import pandas as pd
from scipy.stats import spearmanr


def spearman_between_ranks(df: pd.DataFrame, rank1: str, rank2: str) -> float:
    rho, _ = spearmanr(df[rank1].values, df[rank2].values)
    return float(rho)


def pearson_correlation(x: Mapping[str, float], y: Mapping[str, float]) -> float:
    """Compute Pearson correlation over the intersection of keys.

    This is a small compatibility helper used by the unit tests.
    """

    keys = [k for k in x.keys() if k in y]
    if len(keys) < 2:
        return float("nan")

    xv = pd.Series([float(x[k]) for k in keys], dtype=float)
    yv = pd.Series([float(y[k]) for k in keys], dtype=float)
    return float(xv.corr(yv, method="pearson"))
