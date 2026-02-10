from __future__ import annotations

import json
import time
from dataclasses import asdict
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd
import zarr

from clue_pathway_enrichment.pipeline.config import PipelineConfig
from clue_pathway_enrichment.pipeline.run_pipeline import run


def _iter_with_progress(it, *, total: int, desc: str, enabled: bool):
    """Yield from an iterable, optionally wrapped in a tqdm progress bar (with ETA)."""

    if not enabled:
        yield from it
        return

    try:
        from tqdm import tqdm  # type: ignore
        import sys
    except ModuleNotFoundError:  # pragma: no cover
        yield from it
        return

    # tqdm is much more reliable when writing to stderr (stdout may be buffered)
    yield from tqdm(
        it,
        total=total,
        desc=desc,
        unit="sig",
        mininterval=0.5,
        dynamic_ncols=True,
        file=sys.stderr,
        leave=True,
    )


def _zarr_init_progress(*, enabled: bool):
    """Small helper to show quick progress during Zarr open + metadata reads."""

    if not enabled:
        return None

    try:
        from tqdm import tqdm  # type: ignore
        import sys
    except ModuleNotFoundError:  # pragma: no cover
        return None

    return tqdm(
        total=4,
        desc="Loading zarr matrix",
        unit="step",
        mininterval=0.2,
        dynamic_ncols=True,
        file=sys.stderr,
        leave=True,
    )


# -------------------------
# Public entrypoint
# -------------------------


def scan_signatures(
    *,
    zarr_path: str,
    sig_ids: list[str],
    pipeline_cfg: PipelineConfig,
    threshold: float,
    threshold_metric: str,
    top_fractions: list[float],
    summary_path: str,
    hits_dir: Optional[str],
    save_artifacts_for_hits: bool,
    resume: bool,
    n_workers: int,
    chunk_size: int,
    show_progress: bool = True,
) -> None:
    """
    Scan many signatures:
      - load each signature vector from a Zarr store
      - run the pathway enrichment pipeline on it
      - compute spearman metrics and flag "hits" below a threshold
      - write a summary table, optionally saving artifacts for hits

    Zarr format assumptions:
      root["scores"] : 2D array shape (n_genes, n_sigs)
      root["gene_ids"]: 1D array len n_genes (strings)
      root["sig_ids"] : 1D array len n_sigs (strings)
    """
    if n_workers != 1:
        raise NotImplementedError("Parallel batch scan is not implemented yet. Set batch.execution.n_workers=1.")

    if chunk_size < 1:
        raise ValueError(f"chunk_size must be >= 1. Got: {chunk_size}")

    metric = str(threshold_metric).strip().lower()
    supported_metrics = {"pos_all", "neg_all", "min_all"}
    if metric not in supported_metrics:
        raise ValueError(
            f"Unsupported threshold_metric={threshold_metric!r}. "
            f"Supported: {sorted(supported_metrics)}. "
            f"(Zoom/top-fraction metrics can be added next.)"
        )

    # Prepare summary writing strategy
    summary_path_p = Path(summary_path)
    summary_is_parquet = summary_path_p.suffix.lower() == ".parquet"
    summary_csv_path = (
        summary_path_p.with_suffix(".csv") if summary_is_parquet else summary_path_p
    )  # we always append to CSV
    summary_csv_path.parent.mkdir(parents=True, exist_ok=True)

    # Resume: determine processed sig_ids
    processed = set()
    if resume:
        processed = _load_processed_sig_ids(summary_csv_path, summary_path_p if summary_is_parquet else None)

    # Open Zarr store + read metadata (this can take time on large stores)
    zbar = _zarr_init_progress(enabled=show_progress)
    if zbar is not None:
        zbar.set_postfix_str("open_group")
    root = zarr.open_group(zarr_path, mode="r")
    if zbar is not None:
        zbar.update(1)

    if "scores" not in root:
        if zbar is not None:
            zbar.close()
        raise KeyError(f"Zarr store missing dataset 'scores': {zarr_path}")
    if "gene_ids" not in root:
        if zbar is not None:
            zbar.close()
        raise KeyError(f"Zarr store missing dataset 'gene_ids': {zarr_path}")
    if "sig_ids" not in root:
        if zbar is not None:
            zbar.close()
        raise KeyError(f"Zarr store missing dataset 'sig_ids': {zarr_path}")

    scores = root["scores"]

    if zbar is not None:
        zbar.set_postfix_str("read gene_ids")
    gene_ids = _as_str_list(root["gene_ids"][:])
    if zbar is not None:
        zbar.update(1)

    if zbar is not None:
        zbar.set_postfix_str("read sig_ids")
    store_sig_ids = _as_str_list(root["sig_ids"][:])
    if zbar is not None:
        zbar.update(1)

    if zbar is not None:
        zbar.set_postfix_str("index sig_ids")
    # Map signature id -> column index
    sig_to_col = {s: i for i, s in enumerate(store_sig_ids)}
    if zbar is not None:
        zbar.update(1)
        zbar.close()

    # Prepare hits dir
    hits_dir_p: Optional[Path] = None
    if hits_dir:
        hits_dir_p = Path(hits_dir)
        hits_dir_p.mkdir(parents=True, exist_ok=True)

    # Summary CSV header management
    write_header = not summary_csv_path.exists()

    # Main loop
    n_total = len(sig_ids)
    n_done = 0
    n_skipped = 0
    n_errors = 0
    n_hits = 0

    buffer_rows: list[dict[str, Any]] = []

    iter_ids = _iter_with_progress(
        enumerate(sig_ids),
        total=n_total,
        desc="batch scan",
        enabled=show_progress,
    )

    for idx, sig_id in iter_ids:
        if sig_id in processed:
            n_skipped += 1
            continue

        t0 = time.time()
        row: dict[str, Any] = {
            "sig_id": sig_id,
            "status": "ok",
            "error_msg": None,
            "runtime_sec": None,
            "hit": False,
            "threshold_metric": metric,
            "threshold": float(threshold),
        }

        try:
            col = sig_to_col.get(sig_id)
            if col is None:
                raise KeyError(f"sig_id not found in zarr sig_ids: {sig_id}")

            # Load vector (one column)
            vec = scores[:, col]
            # Ensure we have a real numpy array
            vec_np = np.asarray(vec, dtype=np.float32)

            # Build signature DF expected by pipeline
            sig_df = pd.DataFrame({"gene": gene_ids, "score": vec_np})

            # Run pipeline
            res_df, spearman_by_dir = run(
                signature_path=sig_df,  # your run() now supports DF
                pathway_csv=pipeline_cfg.pathway_csv,
                direction=pipeline_cfg.direction,
                alpha=pipeline_cfg.alpha,
                n_perm=pipeline_cfg.n_perm,
                seed=pipeline_cfg.seed,
                X=pipeline_cfg.X,
                L=pipeline_cfg.L,
                output_spearman_plot=None,
                output_spearman_plot_zoom=None,
                spearman_plot_zoom_top_fraction=float(getattr(pipeline_cfg, "spearman_plot_zoom_top_fraction", 0.1)),
                show_progress=False,  # don't nest progress bars (batch controls the outer loop)
            )

            # Flatten spearman metrics
            sp_pos = _safe_float(spearman_by_dir.get("pos"))
            sp_neg = _safe_float(spearman_by_dir.get("neg"))
            row["spearman_pos"] = sp_pos
            row["spearman_neg"] = sp_neg

            # pathway counts (if direction column exists)
            if "direction" in res_df.columns:
                row["n_pathways_pos"] = int((res_df["direction"] == "pos").sum())
                row["n_pathways_neg"] = int((res_df["direction"] == "neg").sum())
            else:
                # single direction run
                row["n_pathways_pos"] = int(len(res_df)) if pipeline_cfg.direction == "pos" else 0
                row["n_pathways_neg"] = int(len(res_df)) if pipeline_cfg.direction == "neg" else 0

            # Hit decision
            score_for_threshold = _compute_threshold_score(metric, sp_pos, sp_neg)
            row["threshold_score"] = score_for_threshold
            row["hit"] = (score_for_threshold is not None) and (score_for_threshold < threshold)

            # Save artifacts if hit
            if row["hit"]:
                n_hits += 1
                if save_artifacts_for_hits and hits_dir_p is not None:
                    _save_hit_artifacts(
                        hits_dir=hits_dir_p,
                        sig_id=sig_id,
                        res_df=res_df,
                        spearman_by_dir=spearman_by_dir,
                        pipeline_cfg=pipeline_cfg,
                        top_fractions=top_fractions,
                    )

        except Exception as e:
            n_errors += 1
            row["status"] = "error"
            row["error_msg"] = f"{type(e).__name__}: {e}"

        row["runtime_sec"] = round(time.time() - t0, 4)
        buffer_rows.append(row)
        n_done += 1

        # Flush periodically
        if len(buffer_rows) >= chunk_size:
            _append_summary_rows_csv(summary_csv_path, buffer_rows, write_header=write_header)
            write_header = False
            buffer_rows.clear()

            # lightweight progress print (helpful when stderr isn't interactive)
            done_plus_skipped = n_done + n_skipped
            print(
                f"[batch] processed={done_plus_skipped}/{n_total} "
                f"(done={n_done}, skipped={n_skipped}, errors={n_errors}, hits={n_hits})"
            )

    # Flush remainder
    if buffer_rows:
        _append_summary_rows_csv(summary_csv_path, buffer_rows, write_header=write_header)
        buffer_rows.clear()

    # Convert to parquet if requested
    if summary_is_parquet:
        _convert_csv_to_parquet(summary_csv_path, summary_path_p)

    print(
        f"[batch] finished. total={n_total} done={n_done} skipped={n_skipped} "
        f"errors={n_errors} hits={n_hits} summary={summary_path}"
    )


# -------------------------
# Internal helpers
# -------------------------


def _as_str_list(arr: Any) -> list[str]:
    # zarr may return numpy array of dtype object/bytes/unicode
    out: list[str] = []
    for x in arr:
        if isinstance(x, bytes):
            out.append(x.decode("utf-8"))
        else:
            out.append(str(x))
    return out


def _safe_float(x: Any) -> Optional[float]:
    if x is None:
        return None
    try:
        v = float(x)
    except Exception:
        return None
    if np.isnan(v):
        return None
    return v


def _compute_threshold_score(metric: str, sp_pos: Optional[float], sp_neg: Optional[float]) -> Optional[float]:
    if metric == "pos_all":
        return sp_pos
    if metric == "neg_all":
        return sp_neg
    if metric == "min_all":
        vals = [v for v in (sp_pos, sp_neg) if v is not None]
        return min(vals) if vals else None
    raise ValueError(f"Unknown metric: {metric}")


def _append_summary_rows_csv(path: Path, rows: list[dict[str, Any]], *, write_header: bool) -> None:
    df = pd.DataFrame(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, mode="a", header=write_header, index=False)


def _load_processed_sig_ids(summary_csv_path: Path, summary_parquet_path: Optional[Path]) -> set[str]:
    """
    Resume strategy:
      - If CSV exists, read it and take sig_id column.
      - Else if parquet exists, read it, write a CSV version for appending, and take sig_id.
      - Else return empty set.
    """
    if summary_csv_path.exists():
        try:
            df = pd.read_csv(summary_csv_path, usecols=["sig_id"])
            return set(df["sig_id"].dropna().astype(str))
        except Exception:
            # if file exists but is corrupted, don't silently skip everything
            raise

    if summary_parquet_path is not None and summary_parquet_path.exists():
        df = pd.read_parquet(summary_parquet_path, columns=["sig_id"])
        # create a CSV to resume appends
        df.to_csv(summary_csv_path, index=False)
        return set(df["sig_id"].dropna().astype(str))

    return set()


def _convert_csv_to_parquet(csv_path: Path, parquet_path: Path) -> None:
    df = pd.read_csv(csv_path)
    parquet_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(parquet_path, index=False)


def _save_hit_artifacts(
    *,
    hits_dir: Path,
    sig_id: str,
    res_df: pd.DataFrame,
    spearman_by_dir: dict[str, Any],
    pipeline_cfg: PipelineConfig,
    top_fractions: list[float],
) -> None:
    out_dir = hits_dir / sig_id
    out_dir.mkdir(parents=True, exist_ok=True)

    # save result table
    res_df.to_csv(out_dir / "results.csv", index=False)

    # save spearman
    (out_dir / "spearman.json").write_text(json.dumps(spearman_by_dir, indent=2, sort_keys=True), encoding="utf-8")

    # record pipeline params used
    try:
        (out_dir / "pipeline_params.json").write_text(
            json.dumps(asdict(pipeline_cfg), indent=2, sort_keys=True), encoding="utf-8"
        )
    except Exception:
        # PipelineConfig might include non-serializable fields in future; ignore.
        pass

    # NOTE:
    # We are NOT generating plots here yet because:
    # - batch runs should be lightweight
    # - your pipeline's run() currently doesn't return zoom spearman metrics
    # We can add rank-agreement plots for hits in the next step easily
    # by importing save_rank_agreement_plot and writing per-direction images.
    _ = top_fractions  # reserved for future use
