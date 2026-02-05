from __future__ import annotations

import ast
from typing import Dict, Set

import pandas as pd


def _parse_genes_cell(x) -> list[str]:
    """
    Parses a genes cell that may be:
    - a Python list string: "[1, 2, 3]" or "['1','2']"
    - a comma/semicolon/space-separated string: "1,2,3" / "1;2;3" / "1 2 3"
    - already a list
    Returns list of normalized string IDs.
    """
    if x is None or (isinstance(x, float) and pd.isna(x)):
        return []

    if isinstance(x, list):
        return [str(g).strip() for g in x if str(g).strip()]

    s = str(x).strip()
    if not s:
        return []

    # Try literal list parsing first
    if (s.startswith("[") and s.endswith("]")) or (s.startswith("(") and s.endswith(")")):
        try:
            val = ast.literal_eval(s)
            if isinstance(val, (list, tuple)):
                return [str(g).strip() for g in val if str(g).strip()]
        except Exception:
            pass

    # Fallback: split common delimiters
    for delim in [";", ",", "\t", " "]:
        if delim in s:
            parts = [p.strip() for p in s.replace("\t", " ").split(delim)]
            parts = [p for p in parts if p]
            # If we split on space, collapse multiple spaces by re-splitting
            if delim == " ":
                parts = [p for p in s.split() if p.strip()]
            return parts

    # Single token
    return [s]


def load_pathway_mapping_csv(path: str) -> Dict[str, Set[str]]:
    """
    Loads a pathway->genes mapping.

    Supported formats:
    1) Long format: columns include pathway, gene (one row per gene)
    2) Summary format: columns include pathway, and one of: genes / gene_ids
       where that column contains a list or list-like string of gene IDs.

    Returns: dict[pathway] = set(gene_ids_as_strings)
    """
    df = pd.read_csv(path)
    cols = {c.lower().strip(): c for c in df.columns}

    if "pathway" not in cols:
        raise ValueError(f"Pathway mapping must include a 'pathway' column. Found: {list(df.columns)}")

    pathway_c = cols["pathway"]
    mapping: Dict[str, Set[str]] = {}

    # Format 1: long
    if "gene" in cols:
        gene_c = cols["gene"]
        tmp = df[[pathway_c, gene_c]].rename(columns={pathway_c: "pathway", gene_c: "gene"})
        tmp["pathway"] = tmp["pathway"].astype(str).str.strip()
        tmp["gene"] = tmp["gene"].astype(str).str.strip()

        for p, g in zip(tmp["pathway"], tmp["gene"]):
            if p and g:
                mapping.setdefault(p, set()).add(g)
        return mapping

    # Format 2: summary (genes / gene_ids)
    genes_col_key = None
    if "genes" in cols:
        genes_col_key = "genes"
    elif "gene_ids" in cols:
        genes_col_key = "gene_ids"

    if genes_col_key is not None:
        genes_c = cols[genes_col_key]
        tmp = df[[pathway_c, genes_c]].rename(columns={pathway_c: "pathway", genes_c: "genes"})
        tmp["pathway"] = tmp["pathway"].astype(str).str.strip()

        for _, row in tmp.iterrows():
            p = row["pathway"]
            genes = _parse_genes_cell(row["genes"])
            if p:
                mapping.setdefault(p, set()).update(genes)
        return mapping

    raise ValueError(
        "Pathway mapping must have either columns (pathway, gene) or (pathway, genes/gene_ids). "
        f"Found: {list(df.columns)}"
    )
