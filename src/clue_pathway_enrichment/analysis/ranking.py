from __future__ import annotations
import pandas as pd

def add_ranks(df: pd.DataFrame, p_col: str, out_col: str) -> pd.DataFrame:
    """
    Smaller p-value => better rank (rank=1 is best).
    """
    out = df.copy()
    out[out_col] = out[p_col].rank(method="average", ascending=True)
    return out
