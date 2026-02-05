from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional, Union


_ALLOWED_DIRECTIONS = {"pos", "neg", "both"}


@dataclass(frozen=True)
class PipelineConfig:
    """Configuration for `clue_pathway_enrichment.pipeline.run_pipeline.run`.

    Preferred format is JSON:

    ```json
    {
      "pipeline": {
        "signature_path": "signature.tsv",
        "pathway_csv": "pathways.csv",
        "direction": "both",
        "alpha": 0.2,
        "n_perm": 200,
        "seed": 0,
        "X": 1,
        "L": 100,
        "output_csv": "results.csv",
        "output_spearman_json": "spearman.json"
      }
    }
    ```

    TOML is supported as a backward-compatible fallback.
    """

    signature_path: str
    pathway_csv: str

    direction: str = "both"
    alpha: float = 0.2
    n_perm: int = 200
    seed: Optional[int] = 0
    X: int = 1
    L: int = 100

    output_csv: Optional[str] = None
    output_spearman_json: Optional[str] = None

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

        sp = Path(self.signature_path)
        pc = Path(self.pathway_csv)
        if not sp.exists():
            raise FileNotFoundError(f"signature_path does not exist: {sp}")
        if not pc.exists():
            raise FileNotFoundError(f"pathway_csv does not exist: {pc}")


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


def _extract_pipeline_table(raw: Dict[str, Any]) -> Dict[str, Any]:
    table = raw.get("pipeline", {})
    if not isinstance(table, dict):
        raise TypeError("pipeline must be an object/table")
    return table


def load_pipeline_config(path: Union[str, Path]) -> PipelineConfig:
    """Load pipeline config from JSON (preferred) or TOML.

    Selection rules:
    - `.json` => JSON
    - `.toml` => TOML
    - otherwise: try JSON first, then TOML
    """

    p = Path(path)
    suffix = p.suffix.lower()

    if suffix == ".json":
        raw = _load_json(p)
    elif suffix == ".toml":
        raw = _load_toml(p)
    else:
        # best-effort: JSON first, then TOML
        try:
            raw = _load_json(p)
        except Exception:
            raw = _load_toml(p)

    table = _extract_pipeline_table(raw)

    cfg = PipelineConfig(
        signature_path=str(table["signature_path"]),
        pathway_csv=str(table["pathway_csv"]),
        direction=str(table.get("direction", "both")),
        alpha=float(table.get("alpha", 0.2)),
        n_perm=int(table.get("n_perm", 200)),
        seed=_coerce_int_or_none(table.get("seed", 0)),
        X=int(table.get("X", 1)),
        L=int(table.get("L", 100)),
        output_csv=(str(table["output_csv"]) if table.get("output_csv") is not None else None),
        output_spearman_json=(
            str(table["output_spearman_json"]) if table.get("output_spearman_json") is not None else None
        ),
    )
    cfg.validate()
    return cfg
