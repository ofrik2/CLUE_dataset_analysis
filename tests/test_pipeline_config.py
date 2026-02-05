import os
import json
from pathlib import Path

import pytest

from clue_pathway_enrichment.pipeline.config import load_pipeline_config
from clue_pathway_enrichment.pipeline.run_pipeline import run_from_config


def test_load_pipeline_config_and_run_from_config_smoke_json(tmp_path: Path):
    # Minimal signature + pathways files compatible with the package loaders.
    sig_path = tmp_path / "sig.tsv"
    sig_path.write_text("gene\tscore\nA\t1\nB\t-1\n")

    pw_path = tmp_path / "pathways.csv"
    pw_path.write_text("pathway,gene\np1,A\np1,X\np2,B\n")

    cfg_path = tmp_path / "config.json"
    cfg_path.write_text(
        json.dumps(
            {
                "pipeline": {
                    "signature_path": "sig.tsv",
                    "pathway_csv": "pathways.csv",
                    "direction": "both",
                    "alpha": 0.2,
                    "n_perm": 0,
                    "seed": 0,
                    "X": 1,
                    "L": 100,
                }
            }
        )
    )

    # Use relative paths by setting cwd to tmp_path for this test
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        cfg = load_pipeline_config(cfg_path)
        assert cfg.signature_path.endswith("sig.tsv")
        assert cfg.pathway_csv.endswith("pathways.csv")

        df, corr = run_from_config(str(cfg_path), n_perm=0)
        assert hasattr(df, "shape")
        assert isinstance(corr, dict)
    finally:
        os.chdir(cwd)


def test_load_pipeline_config_invalid_direction_json(tmp_path: Path):
    sig_path = tmp_path / "sig.tsv"
    sig_path.write_text("gene\tscore\nA\t1\n")

    pw_path = tmp_path / "pathways.csv"
    pw_path.write_text("pathway,gene\np1,A\n")

    cfg_path = tmp_path / "config.json"
    cfg_path.write_text(
        json.dumps(
            {
                "pipeline": {
                    "signature_path": "sig.tsv",
                    "pathway_csv": "pathways.csv",
                    "direction": "sideways",
                }
            }
        )
    )

    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        with pytest.raises(ValueError):
            load_pipeline_config(cfg_path)
    finally:
        os.chdir(cwd)


def test_load_pipeline_config_toml_backward_compatible(tmp_path: Path):
    sig_path = tmp_path / "sig.tsv"
    sig_path.write_text("gene\tscore\nA\t1\n")

    pw_path = tmp_path / "pathways.csv"
    pw_path.write_text("pathway,gene\np1,A\n")

    cfg_path = tmp_path / "config.toml"
    cfg_path.write_text(
        """
[pipeline]
signature_path = "sig.tsv"
pathway_csv = "pathways.csv"
direction = "both"
""".strip()
    )

    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        cfg = load_pipeline_config(cfg_path)
        assert cfg.direction == "both"
    finally:
        os.chdir(cwd)
