from pathlib import Path
import json
import pytest

# ...existing code...

def _assert_file_readable(p: Path):
    assert p.exists(), f"File does not exist: {p}"
    assert p.is_file(), f"Not a file: {p}"
    with open(p, "rb") as fh:
        chunk = fh.read(1024)
    assert len(chunk) > 0, f"Unable to read bytes from file: {p}"

def test_load_pathway_mapping_files(pathway_files):
    """
    For each provided pathway file, call load_pathway_mapping_csv and assert structure.
    Prints a short summary so uploaded results can be inspected in pytest output.
    """
    from clue_pathway_enrichment.io.load_pathways import load_pathway_mapping_csv

    for p in pathway_files:
        _assert_file_readable(p)
        mapping = load_pathway_mapping_csv(str(p))
        assert isinstance(mapping, dict), f"Expected dict mapping, got {type(mapping)} for file {p}"
        total_pathways = len(mapping)
        total_genes = sum(len(s) for s in mapping.values())
        sample_keys = list(mapping.keys())[:5]
        print(f"Loaded pathways from {p}: {total_pathways} pathways, {total_genes} total genes, sample pathways: {sample_keys}")

def test_load_signature_files(signature_files):
    """
    For each provided signature file, call load_signature_table and assert DataFrame contains gene and score.
    Prints head() so uploaded results can be inspected.
    """
    import pandas as pd
    from clue_pathway_enrichment.io.load_signature import load_signature_table

    for p in signature_files:
        _assert_file_readable(p)
        df = load_signature_table(str(p))
        assert isinstance(df, pd.DataFrame), f"Expected DataFrame, got {type(df)} for file {p}"
        assert "gene" in df.columns and "score" in df.columns, f"Missing required columns in signature file {p}"
        # print a concise preview
        preview = df.head(10).to_dict(orient="records")
        print(f"Loaded signature from {p}: {len(df)} rows, preview: {preview}")
