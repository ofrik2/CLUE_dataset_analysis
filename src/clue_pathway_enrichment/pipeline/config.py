from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional, Union

_ALLOWED_DIRECTIONS = {"pos", "neg", "both"}

'''
NOTE: this config loader still handles a single signature run. 
the pipeline parameters and run parameters are separated just for the option to use batch mode
it does not have anything to do with this config loader

'''


# -------------------------
# Dataclasses
# -------------------------


@dataclass(frozen=True)
class PipelineConfig:
    """
    Algorithm + shared parameters for `clue_pathway_enrichment.pipeline.run_pipeline.run`.

    New JSON format:

    {
      "pipeline": {
        "pathway_csv": "...",
        "direction": "both",
        "alpha": 0.2,
        "n_perm": 200,
        "seed": 0,
        "X": 1,
        "L": 100,
        "spearman_plot_zoom_top_fraction": 0.1,
        "show_progress": true
      },
      "single_run": {
        "signature_path": "...",
        "output_csv": "...",
        "output_spearman_json": "...",
        "output_spearman_plot": "...",
        "output_spearman_plot_zoom": "..."
      },
      "batch": { ... }
    }

    TOML is supported as a backward-compatible fallback.
    """

    pathway_csv: str

    direction: str = "both"
    alpha: float = 0.2
    n_perm: int = 200
    seed: Optional[int] = 0
    X: int = 1
    L: int = 100

    # zoomed rank agreement plot
    spearman_plot_zoom_top_fraction: float = 0.1

    # UX / logging
    show_progress: bool = True

    def validate(self) -> None:
        if self.direction not in _ALLOWED_DIRECTIONS:
            raise ValueError(f"direction must be one of {_ALLOWED_DIRECTIONS}, got: {self.direction}")
        if self.alpha <= 0 or self.alpha >= 1:
            raise ValueError(f"alpha must be in (0, 1), got: {self.alpha}")
        if self.n_perm < 0:
            raise ValueError(f"n_perm must be >= 0, got: {self.n_perm}")
        if self.X < 1:
            raise ValueError(f"X must be >= 1, got: {self.X}")
        if self.L < 1:
            raise ValueError(f"L must be >= 1, got: {self.L}")
        if not (0 < float(self.spearman_plot_zoom_top_fraction) <= 1):
            raise ValueError(
                f"spearman_plot_zoom_top_fraction must be in (0, 1], got: {self.spearman_plot_zoom_top_fraction}"
            )

        pc = Path(self.pathway_csv)
        if not pc.exists():
            raise FileNotFoundError(f"pathway_csv does not exist: {pc}")


@dataclass(frozen=True)
class SingleRunConfig:
    """
    Single-run I/O for `pipeline/cli.py`.

    In the new JSON format, this lives under "single_run".
    For backward compatibility, if "single_run" is missing we will also accept
    these fields under "pipeline".
    """

    signature_path: str

    output_csv: Optional[str] = None
    output_spearman_json: Optional[str] = None
    output_spearman_plot: Optional[str] = None
    output_spearman_plot_zoom: Optional[str] = None

    def validate(self) -> None:
        sp = Path(self.signature_path)
        if not sp.exists():
            raise FileNotFoundError(f"signature_path does not exist: {sp}")


# -------------------------
# Helpers
# -------------------------


def _coerce_int_or_none(v: Any) -> Optional[int]:
    if v is None:
        return None
    if isinstance(v, bool):
        raise TypeError("seed must be int or null")
    return int(v)


def _load_toml(path: Path) -> Dict[str, Any]:
    """Load TOML with stdlib tomllib (3.11+) or tomli fallback."""
    import importlib

    try:
        toml_loader = importlib.import_module("tomllib")
    except ModuleNotFoundError:  # pragma: no cover
        toml_loader = importlib.import_module("tomli")

    with path.open("rb") as fh:
        return toml_loader.load(fh)


def _load_json(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def _load_raw(path: Union[str, Path]) -> Dict[str, Any]:
    """Load config file content as a dict (JSON preferred; TOML supported)."""
    p = Path(path)
    suffix = p.suffix.lower()

    if suffix == ".json":
        return _load_json(p)
    if suffix == ".toml":
        return _load_toml(p)

    # best-effort: JSON first, then TOML
    try:
        return _load_json(p)
    except Exception:
        return _load_toml(p)


def _extract_table(raw: Dict[str, Any], key: str) -> Dict[str, Any]:
    table = raw.get(key, {})
    if not isinstance(table, dict):
        raise TypeError(f"{key} must be an object/table")
    return table


def _extract_single_run_table(raw: Dict[str, Any]) -> Dict[str, Any]:
    """
    New format: raw["single_run"].
    Backward compatibility: if missing, fall back to raw["pipeline"] for signature/output keys.
    """
    if "single_run" in raw:
        table = raw.get("single_run", {})
        if not isinstance(table, dict):
            raise TypeError("single_run must be an object/table")
        return table

    # legacy fallback
    pipeline_table = raw.get("pipeline", {})
    if not isinstance(pipeline_table, dict):
        raise TypeError("pipeline must be an object/table")
    return pipeline_table


# -------------------------
# Public loaders
# -------------------------


def load_pipeline_config(path: Union[str, Path]) -> PipelineConfig:
    """
    Load only the 'pipeline' section (algorithm/shared params).

    This is used by both:
    - pipeline/cli.py (single run)
    - batch/cli.py (batch run)
    """
    raw = _load_raw(path)
    table = _extract_table(raw, "pipeline")

    cfg = PipelineConfig(
        pathway_csv=str(table["pathway_csv"]),
        direction=str(table.get("direction", "both")),
        alpha=float(table.get("alpha", 0.2)),
        n_perm=int(table.get("n_perm", 200)),
        seed=_coerce_int_or_none(table.get("seed", 0)),
        X=int(table.get("X", 1)),
        L=int(table.get("L", 100)),
        spearman_plot_zoom_top_fraction=float(table.get("spearman_plot_zoom_top_fraction", 0.1)),
        show_progress=bool(table.get("show_progress", True)),
    )
    cfg.validate()
    return cfg


def load_single_run_config(path: Union[str, Path]) -> SingleRunConfig:
    """
    Load single-run I/O section.

    New format: "single_run".
    Legacy fallback: read from "pipeline".
    """
    raw = _load_raw(path)
    table = _extract_single_run_table(raw)

    cfg = SingleRunConfig(
        signature_path=str(table["signature_path"]),
        output_csv=(str(table["output_csv"]) if table.get("output_csv") is not None else None),
        output_spearman_json=(
            str(table["output_spearman_json"]) if table.get("output_spearman_json") is not None else None
        ),
        output_spearman_plot=(
            str(table["output_spearman_plot"]) if table.get("output_spearman_plot") is not None else None
        ),
        output_spearman_plot_zoom=(
            str(table["output_spearman_plot_zoom"]) if table.get("output_spearman_plot_zoom") is not None else None
        ),
    )
    cfg.validate()
    return cfg
