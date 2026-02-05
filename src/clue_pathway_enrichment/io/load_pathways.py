from __future__ import annotations
import pandas as pd
from typing import Dict, Set

def load_pathway_mapping_csv(path: str) -> Dict[str, Set[str]]:
    """
    Loads a pathway->genes mapping.
    Expected long format: pathway,gene
    """
    df = pd.read_csv(path)
    cols = {c.lower().strip(): c for c in df.columns}
    if "pathway" not in cols or "gene" not in cols:
        raise ValueError(f"Pathway mapping must have columns: pathway, gene. Found: {list(df.columns)}")

    df = df[[cols["pathway"], cols["gene"]]].rename(columns={cols["pathway"]: "pathway", cols["gene"]: "gene"})
    df["pathway"] = df["pathway"].astype(str).str.strip()
    df["gene"] = df["gene"].astype(str).str.strip()

    mapping: Dict[str, Set[str]] = {}
    for p, g in zip(df["pathway"], df["gene"]):
        mapping.setdefault(p, set()).add(g)
    return mapping
