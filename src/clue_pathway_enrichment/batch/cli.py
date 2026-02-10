from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Optional, Any

from clue_pathway_enrichment.pipeline.config import load_pipeline_config, PipelineConfig
from clue_pathway_enrichment.batch.scan_signatures import scan_signatures


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Batch scan CLUE signatures and flag low Spearman agreement")
    p.add_argument("--config", required=True, help="Path to config (.json preferred)")
    p.add_argument("--max-signatures", type=int, default=None, help="Optional cap for debugging (overrides config)")
    p.add_argument("--no-progress", action="store_true", help="Disable the batch progress bar / ETA")
    return p


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _read_sig_ids(path: str) -> list[str]:
    p = Path(path)
    lines = p.read_text(encoding="utf-8").splitlines()
    sigs = [ln.strip() for ln in lines if ln.strip()]

    # preserve order, remove duplicates
    seen: set[str] = set()
    out: list[str] = []
    for s in sigs:
        if s not in seen:
            seen.add(s)
            out.append(s)
    return out


def _get_batch_table(raw: dict[str, Any]) -> dict[str, Any]:
    batch = raw.get("batch")
    if batch is None:
        raise KeyError("Missing required top-level key: 'batch'")
    if not isinstance(batch, dict):
        raise TypeError("'batch' must be an object")
    return batch


def _get_signature_store_path(batch: dict) -> str:
    store = batch.get("signature_store")
    if not isinstance(store, dict):
        raise TypeError("'batch.signature_store' must be an object")
    store_type = str(store.get("type", "")).lower()
    if store_type != "zarr":
        raise ValueError(f"Only signature_store.type='zarr' supported. Got: {store_type!r}")
    path = store.get("path")
    if not path:
        raise KeyError("Missing 'batch.signature_store.path'")
    return str(path)



def _cli_loading_progress(*, enabled: bool, total: int):
    """Finite progress indicator for batch CLI setup steps."""

    if not enabled:
        return None

    try:
        from tqdm import tqdm  # type: ignore
        import sys
    except ModuleNotFoundError:  # pragma: no cover
        return None

    return tqdm(
        total=total,
        desc="batch cli",
        unit="step",
        mininterval=0.1,
        dynamic_ncols=True,
        file=sys.stderr,
        leave=True,
    )

# src/clue_pathway_enrichment/batch/cli.py

# src/clue_pathway_enrichment/batch/cli.py

def main(argv: Optional[list[str]] = None) -> int:
    args = _build_parser().parse_args(argv)
    cfg_path = Path(args.config)

    raw = _read_json(cfg_path)
    batch = _get_batch_table(raw)

    execution = batch.get("execution", {})
    if not isinstance(execution, dict):
        raise TypeError("'batch.execution' must be an object")

    show_progress_cfg = execution.get("show_progress", None)
    show_progress_cli = bool(show_progress_cfg) if show_progress_cfg is not None else True
    show_progress_cli = show_progress_cli and (not args.no_progress)

    # Keep setup steps aligned with updates below.
    _SETUP_STEPS = 6
    pbar = _cli_loading_progress(enabled=show_progress_cli, total=_SETUP_STEPS)

    def _pbar_postfix(msg: str) -> None:
        if pbar is not None:
            pbar.set_postfix_str(msg)

    def _pbar_step() -> None:
        if pbar is not None:
            pbar.update(1)

    _pbar_postfix("load pipeline config")
    pipeline_cfg: PipelineConfig = load_pipeline_config(cfg_path)
    _pbar_step()

    _pbar_postfix("get zarr path")
    zarr_path = _get_signature_store_path(batch)
    _pbar_step()

    _pbar_postfix("read sig ids")
    sig_ids_path = str(batch["sig_ids_path"])
    sig_ids = _read_sig_ids(sig_ids_path)
    _pbar_step()

    _pbar_postfix("parse thresholding")
    threshold = float(batch.get("threshold", 0.2))
    threshold_metric = str(batch.get("threshold_metric", "min_all"))
    top_fractions = batch.get("top_fractions", [0.1])
    if top_fractions is None:
        top_fractions = []
    if not isinstance(top_fractions, list):
        raise TypeError("'batch.top_fractions' must be a list (e.g., [0.1])")
    _pbar_step()

    _pbar_postfix("parse outputs")
    outputs = batch.get("outputs", {})
    if not isinstance(outputs, dict):
        raise TypeError("'batch.outputs' must be an object")
    summary_path = str(outputs["summary_path"])
    hits_dir = str(outputs.get("hits_dir", "")) or None
    save_artifacts_for_hits = bool(outputs.get("save_artifacts_for_hits", True))
    _pbar_step()

    _pbar_postfix("parse execution")
    resume = bool(execution.get("resume", True))
    n_workers = int(execution.get("n_workers", 1))
    chunk_size = int(execution.get("chunk_size", 50))
    max_signatures_cfg = execution.get("max_signatures", None)

    if show_progress_cfg is None:
        show_progress = bool(getattr(pipeline_cfg, "show_progress", True))
    else:
        show_progress = bool(show_progress_cfg)
    show_progress = show_progress and (not args.no_progress)
    _pbar_step()

    if pbar is not None:
        pbar.close()

    print("[batch] zarr_path:", zarr_path, flush=True)
    print("[batch] sig_ids_path:", sig_ids_path, flush=True)

    if args.max_signatures is not None:
        sig_ids = sig_ids[: args.max_signatures]
    elif max_signatures_cfg is not None:
        sig_ids = sig_ids[: int(max_signatures_cfg)]

    scan_signatures(
        zarr_path=zarr_path,
        sig_ids=sig_ids,
        pipeline_cfg=pipeline_cfg,
        threshold=threshold,
        threshold_metric=threshold_metric,
        top_fractions=[float(x) for x in top_fractions],
        summary_path=summary_path,
        hits_dir=hits_dir,
        save_artifacts_for_hits=save_artifacts_for_hits,
        resume=resume,
        n_workers=n_workers,
        chunk_size=chunk_size,
        show_progress=show_progress,
    )

    return 0




if __name__ == "__main__":
    raise SystemExit(main())
