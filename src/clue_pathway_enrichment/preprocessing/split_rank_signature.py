from __future__ import annotations

import pandas as pd


_ALLOWED_SIGNATURE_RANKING = {"signed_split", "abs"}


def rank_signature(sig: pd.DataFrame, mode: str = "signed_split") -> dict[str, pd.DataFrame]:
    """
    mode:
      - "signed_split": returns {"pos": ..., "neg": ...}
      - "abs": returns {"abs": ...} sorted by abs(score) desc
    """
    df = sig.copy()
    df["gene"] = df["gene"].astype(str).str.strip()

    mode = str(mode).strip().lower()
    if mode not in _ALLOWED_SIGNATURE_RANKING:
        raise ValueError(f"Unknown signature ranking mode: {mode!r}. Allowed: {sorted(_ALLOWED_SIGNATURE_RANKING)}")

    if mode == "signed_split":
        pos = df[df["score"] > 0].sort_values("score", ascending=False).reset_index(drop=True)
        neg = df[df["score"] < 0].sort_values("score", ascending=True).reset_index(drop=True)
        return {"pos": pos, "neg": neg}

    # mode == "abs"
    out = df[df["score"] != 0].copy()
    out["__abs__"] = out["score"].abs()
    # Tie-breaker: if abs equal, keep positive ahead of negative by sorting score desc as 2nd key
    out = out.sort_values(["__abs__", "score"], ascending=[False, False]).drop(columns="__abs__")
    out = out.reset_index(drop=True)
    return {"abs": out}


# Backward-compatible wrapper (so any old imports won't break)
def split_and_rank_signature(sig: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    ranked = rank_signature(sig, mode="signed_split")
    return ranked["pos"], ranked["neg"]