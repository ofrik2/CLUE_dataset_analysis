# tests/test_io_smoke.py
from __future__ import annotations

import os
from pathlib import Path

import pytest

from clue_pathway_enrichment.io.load_pathways import load_pathway_mapping_csv
from clue_pathway_enrichment.io.load_signature import load_signature_table


def _require_env_path(var_name: str) -> str:
    """
    Require an env var to point to an existing file.
    If not provided, skip (so CI won't fail by default).
    """
    p = os.getenv(var_name)
    if not p:
        pytest.skip(f"{var_name} not set. Provide it to run this smoke test.")
    path = Path(p).expanduser()
    if not path.exists():
        pytest.fail(f"{var_name} points to missing file: {path}")
    if not path.is_file():
        pytest.fail(f"{var_name} is not a file: {path}")
    return str(path)


def test_load_signature_table_smoke():
    sig_path = _require_env_path("CLUE_SIGNATURE_PATH")
    df = load_signature_table(sig_path)

    # Basic sanity
    assert df is not None
    assert len(df) > 0, "Signature loaded but is empty"
    assert set(df.columns) == {"gene", "score"}

    # More sanity
    assert df["gene"].astype(str).str.len().gt(0).all(), "Found empty gene names"
    assert df["score"].notna().all(), "Found NaN scores after loading"


def test_load_pathway_mapping_smoke():
    map_path = _require_env_path("CLUE_PATHWAY_MAPPING_PATH")
    mapping = load_pathway_mapping_csv(map_path)

    # Basic sanity
    assert mapping is not None
    assert isinstance(mapping, dict)
    assert len(mapping) > 0, "Pathway mapping loaded but is empty"

    # More sanity: every pathway has >=1 gene, no empty strings
    for pathway, genes in mapping.items():
        assert str(pathway).strip() != "", "Found empty pathway name"
        assert isinstance(genes, set)
        assert len(genes) > 0, f"Pathway '{pathway}' has no genes"
        assert all(str(g).strip() != "" for g in genes), f"Pathway '{pathway}' contains empty gene"


def test_signature_genes_overlap_with_pathways():
    """
    Optional: checks there is at least some overlap between signature genes and pathway genes.
    This catches cases where one file is in a different gene-id namespace than the other.
    """
    sig_path = _require_env_path("CLUE_SIGNATURE_PATH")
    map_path = _require_env_path("CLUE_PATHWAY_MAPPING_PATH")

    df = load_signature_table(sig_path)
    mapping = load_pathway_mapping_csv(map_path)

    sig_genes = set(df["gene"].astype(str))
    pathway_genes = set().union(*mapping.values()) if mapping else set()

    overlap = sig_genes & pathway_genes
    assert len(overlap) > 0, (
        "No overlap between signature genes and pathway genes. "
        "Possible mismatch (e.g., symbols vs Entrez vs Ensembl)."
    )
