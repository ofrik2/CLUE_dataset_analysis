from __future__ import annotations
import pandas as pd

def split_and_rank_signature(sig: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Returns (pos_ranked, neg_ranked)
    pos_ranked: score > 0 sorted desc
    neg_ranked: score < 0 sorted by most negative first (i.e., score asc)
    """
    df = sig.copy()
    df["gene"] = df["gene"].astype(str).str.strip()

    pos = df[df["score"] > 0].sort_values("score", ascending=False).reset_index(drop=True)
    neg = df[df["score"] < 0].sort_values("score", ascending=True).reset_index(drop=True)  # most negative first
    return pos, neg
