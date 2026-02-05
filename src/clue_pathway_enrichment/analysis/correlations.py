from __future__ import annotations
import pandas as pd
from scipy.stats import spearmanr

def spearman_between_ranks(df: pd.DataFrame, rank1: str, rank2: str) -> float:
    rho, _ = spearmanr(df[rank1].values, df[rank2].values)
    return float(rho)
