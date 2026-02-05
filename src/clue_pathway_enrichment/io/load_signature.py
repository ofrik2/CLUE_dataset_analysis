from __future__ import annotations
import pandas as pd

def load_signature_table(path: str) -> pd.DataFrame:
    """
    Loads a 2-column signature table: gene, score
    Supports CSV/TSV by extension.
    """
    sep = "\t" if path.lower().endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(path, sep=sep)

    # Normalize column names
    cols = {c.lower().strip(): c for c in df.columns}
    if "gene" not in cols or "score" not in cols:
        raise ValueError(f"Signature file must have columns: gene, score. Found: {list(df.columns)}")

    df = df[[cols["gene"], cols["score"]]].rename(columns={cols["gene"]: "gene", cols["score"]: "score"})
    df["gene"] = df["gene"].astype(str).str.strip()
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df = df.dropna(subset=["gene", "score"])
    return df
