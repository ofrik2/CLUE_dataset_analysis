from __future__ import annotations
import pandas as pd
from typing import Optional


def load_signature_table(path: str, score_col: Optional[str] = None) -> pd.DataFrame:
    """
    Loads a signature table into a standardized 2-column DataFrame: gene, score

    Supported:
    - Standard format: columns include 'gene' and 'score'
    - CLUE-style single signature: columns include 'gene_id' and one score column
      (or provide score_col explicitly)

    CSV/TSV detected by extension.
    """
    sep = "\t" if path.lower().endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(path, sep=sep)

    cols = {c.lower().strip(): c for c in df.columns}

    # Case 1: already gene/score
    if "gene" in cols and "score" in cols:
        gene_c = cols["gene"]
        score_c = cols["score"]

    # Case 2: CLUE-style gene_id + signature column
    elif "gene_id" in cols:
        gene_c = cols["gene_id"]

        if score_col is not None:
            if score_col not in df.columns:
                raise ValueError(f"score_col='{score_col}' not found. Found: {list(df.columns)}")
            score_c = score_col
        else:
            other_cols = [c for c in df.columns if c != gene_c]
            if len(other_cols) != 1:
                raise ValueError(
                    "CLUE-style signature expected exactly 2 columns: gene_id + one score column. "
                    f"Found columns: {list(df.columns)}"
                )
            score_c = other_cols[0]
    else:
        raise ValueError(
            "Signature file must have columns: (gene, score) OR (gene_id, <signature_col>). "
            f"Found: {list(df.columns)}"
        )

    out = df[[gene_c, score_c]].rename(columns={gene_c: "gene", score_c: "score"})
    out["gene"] = out["gene"].astype(str).str.strip()
    out["score"] = pd.to_numeric(out["score"], errors="coerce")
    out = out.dropna(subset=["gene", "score"])
    return out

