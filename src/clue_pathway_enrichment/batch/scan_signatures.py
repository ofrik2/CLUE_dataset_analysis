from __future__ import annotations

import json
import time
from dataclasses import asdict
from pathlib import Path
from typing import List, Dict, Tuple, Any, Optional, Iterable

import numpy as np
import pandas as pd
import zarr
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import multiprocessing as mp


from clue_pathway_enrichment.pipeline.config import PipelineConfig
from clue_pathway_enrichment.pipeline.run_pipeline import run

from concurrent.futures import ProcessPoolExecutor, as_completed


# -----------------------
# Multiprocessing helpers
# -----------------------





def _chunks(xs: list[str], chunk_size: int) -> Iterable[list[str]]:
    for i in range(0, len(xs), chunk_size):
        yield xs[i:i + chunk_size]


def _pipeline_kwargs_from_cfg(pipeline_cfg) -> dict[str, Any]:
    """
    Convert PipelineConfig to kwargs for run().
    """
    return dict(
        pathway_csv=pipeline_cfg.pathway_csv,
        direction=pipeline_cfg.direction,
        signature_ranking=pipeline_cfg.signature_ranking,
        alpha=pipeline_cfg.alpha,
        n_perm=pipeline_cfg.n_perm,
        seed=pipeline_cfg.seed,
        X=pipeline_cfg.X,
        L=pipeline_cfg.L,
        spearman_plot_zoom_top_fraction=float(
            getattr(pipeline_cfg, "spearman_plot_zoom_top_fraction", 0.1)
        ),
        show_progress=False,
    )


def _process_block_worker(
    *,
    zarr_path: str,
    sig_ids_block: list[str],
    pipeline_kwargs: dict[str, Any],
    metric: str,
    threshold: float,
) -> tuple[list[dict[str, Any]], list[str], int]:
    """
    Worker computes summary rows for a block of sig_ids.
    Returns: (rows, hit_sig_ids, n_errors)
    """
    # local import inside worker is ok; avoids some spawn quirks
    from clue_pathway_enrichment.pipeline.run_pipeline import run

    root = zarr.open_group(zarr_path, mode="r")

    scores = root["scores"]
    gene_ids = _as_str_list(root["gene_ids"][:])
    store_sig_ids = _as_str_list(root["sig_ids"][:])
    sig_to_col = {s: i for i, s in enumerate(store_sig_ids)}

    # Determine the minimal contiguous slice that covers this block.
    # Even if there are gaps (due to resume), reading a slightly wider slice is OK.
    cols = []
    for s in sig_ids_block:
        c = sig_to_col.get(s)
        if c is None:
            cols.append(None)
        else:
            cols.append(c)

    valid_cols = [c for c in cols if c is not None]
    if not valid_cols:
        # nothing valid in this block
        rows = []
        for s in sig_ids_block:
            rows.append({
                "sig_id": s,
                "status": "error",
                "error_msg": "KeyError: sig_id not found in zarr sig_ids",
                "runtime_sec": 0.0,
                "hit": False,
                "threshold_metric": metric,
                "threshold": float(threshold),
                "spearman_pos": None,
                "spearman_neg": None,
                "threshold_score": None,
            })
        return rows, [], len(sig_ids_block)

    j0, j1 = min(valid_cols), max(valid_cols) + 1
    mat = np.asarray(scores[:, j0:j1], dtype=np.float32)  # (genes, width)

    # Reuse one DF to avoid rebuilding index every signature (big speed win)
    sig_df = pd.DataFrame({"gene": gene_ids, "score": np.zeros(len(gene_ids), dtype=np.float32)})

    rows: list[dict[str, Any]] = []
    hit_ids: list[str] = []
    n_errors = 0

    for s, c in zip(sig_ids_block, cols):
        t0 = time.time()
        row: dict[str, Any] = {
            "sig_id": s,
            "status": "ok",
            "error_msg": None,
            "runtime_sec": None,
            "hit": False,
            "threshold_metric": metric,
            "threshold": float(threshold),
        }

        try:
            if c is None:
                raise KeyError(f"sig_id not found in zarr sig_ids: {s}")

            vec = mat[:, c - j0]
            sig_df["score"].values[:] = vec

            res_df, spearman_by_dir = run(signature_path=sig_df, **pipeline_kwargs)

            sp_pos = _safe_float(spearman_by_dir.get("pos"))
            sp_neg = _safe_float(spearman_by_dir.get("neg"))
            sp_abs = _safe_float(spearman_by_dir.get("abs"))

            row["spearman_pos"] = sp_pos
            row["spearman_neg"] = sp_neg
            row["spearman_abs"] = sp_abs

            if isinstance(res_df, pd.DataFrame) and "direction" in res_df.columns:
                row["n_pathways_pos"] = int((res_df["direction"] == "pos").sum())
                row["n_pathways_neg"] = int((res_df["direction"] == "neg").sum())
                row["n_pathways_abs"] = int((res_df["direction"] == "abs").sum())
            else:
                row["n_pathways_pos"] = None
                row["n_pathways_neg"] = None
                row["n_pathways_abs"] = None

            score_for_threshold = _compute_threshold_score(metric, sp_pos, sp_neg, sp_abs)

            row["threshold_score"] = score_for_threshold
            row["hit"] = (score_for_threshold is not None) and (score_for_threshold < threshold)

            if row["hit"]:
                hit_ids.append(s)

        except Exception as e:
            n_errors += 1
            row["status"] = "error"
            row["error_msg"] = f"{type(e).__name__}: {e}"
            row["spearman_pos"] = None
            row["spearman_neg"] = None
            row["spearman_abs"] = None
            row["threshold_score"] = None
            row["hit"] = False

        row["runtime_sec"] = round(time.time() - t0, 4)
        rows.append(row)

    return rows, hit_ids, n_errors



# -----------------------
# Progress bars helpers
# -----------------------

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
      root["gene_ids"]: 1D array len n_genes (strings/bytes)
      root["sig_ids"] : 1D array len n_sigs (strings/bytes)
    """
    if chunk_size < 1:
        raise ValueError(f"chunk_size must be >= 1. Got: {chunk_size}")

    metric = str(threshold_metric).strip().lower()
    supported_metrics = {"pos_all", "neg_all", "min_all", "abs_all"}

    if metric not in supported_metrics:
        raise ValueError(
            f"Unsupported threshold_metric={threshold_metric!r}. "
            f"Supported: {sorted(supported_metrics)}."
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
        processed = _load_processed_sig_ids(
            summary_csv_path,
            summary_path_p if summary_is_parquet else None,
        )

    # Open Zarr store + read metadata
    zbar = _zarr_init_progress(enabled=show_progress)
    if zbar is not None:
        zbar.set_postfix_str("open_group")
    root = zarr.open_group(zarr_path, mode="r")
    if zbar is not None:
        zbar.update(1)

    for key in ("scores", "gene_ids", "sig_ids"):
        if key not in root:
            if zbar is not None:
                zbar.close()
            raise KeyError(f"Zarr store missing dataset {key!r}: {zarr_path}")

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

    # Filter to do list (resume)
    todo_sig_ids = [s for s in sig_ids if s not in processed]
    n_total = len(sig_ids)
    n_todo = len(todo_sig_ids)

    n_done = 0
    n_skipped = n_total - n_todo
    n_errors = 0
    n_hits = 0

    # ---------- SINGLE PROCESS PATH ----------
    if n_workers <= 1:
        buffer_rows: list[dict[str, Any]] = []

        iter_ids = _iter_with_progress(
            enumerate(sig_ids),
            total=n_total,
            desc="batch scan",
            enabled=show_progress,
        )

        for idx, sig_id in iter_ids:
            if sig_id in processed:
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

                vec_np = np.asarray(scores[:, col], dtype=np.float32)

                sig_df = pd.DataFrame({"gene": gene_ids, "score": vec_np})

                res_df, spearman_by_dir = run(
                    signature_path=sig_df,
                    pathway_csv=pipeline_cfg.pathway_csv,
                    direction=pipeline_cfg.direction,
                    signature_ranking=pipeline_cfg.signature_ranking,
                    alpha=pipeline_cfg.alpha,
                    n_perm=pipeline_cfg.n_perm,
                    seed=pipeline_cfg.seed,
                    X=pipeline_cfg.X,
                    L=pipeline_cfg.L,
                    spearman_plot_zoom_top_fraction=float(
                        getattr(pipeline_cfg, "spearman_plot_zoom_top_fraction", 0.1)),
                    show_progress=False,
                )

                sp_pos = _safe_float(spearman_by_dir.get("pos"))
                sp_neg = _safe_float(spearman_by_dir.get("neg"))
                sp_abs = _safe_float(spearman_by_dir.get("abs"))

                row["spearman_pos"] = sp_pos
                row["spearman_neg"] = sp_neg
                row["spearman_abs"] = sp_abs

                if "direction" in res_df.columns:
                    row["n_pathways_pos"] = int((res_df["direction"] == "pos").sum())
                    row["n_pathways_neg"] = int((res_df["direction"] == "neg").sum())
                    row["n_pathways_abs"] = int((res_df["direction"] == "abs").sum())
                else:
                    row["n_pathways_pos"] = 0
                    row["n_pathways_neg"] = 0
                    row["n_pathways_abs"] = int(len(res_df)) if pipeline_cfg.signature_ranking == "abs" else 0

                score_for_threshold = _compute_threshold_score(metric, sp_pos, sp_neg, sp_abs)

                row["threshold_score"] = score_for_threshold
                row["hit"] = (score_for_threshold is not None) and (score_for_threshold < threshold)

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

            if len(buffer_rows) >= chunk_size:
                _append_summary_rows_csv(summary_csv_path, buffer_rows, write_header=write_header)
                write_header = False
                buffer_rows.clear()

                done_plus_skipped = n_done + n_skipped
                print(
                    f"[batch] processed={done_plus_skipped}/{n_total} "
                    f"(done={n_done}, skipped={n_skipped}, errors={n_errors}, hits={n_hits})"
                )

        if buffer_rows:
            _append_summary_rows_csv(summary_csv_path, buffer_rows, write_header=write_header)

        if summary_is_parquet:
            _convert_csv_to_parquet(summary_csv_path, summary_path_p)

        print(
            f"[batch] finished. total={n_total} done={n_done} skipped={n_skipped} "
            f"errors={n_errors} hits={n_hits} summary={summary_path}"
        )
        return

    # ---------- MULTIPROCESS PATH ----------
    pipeline_kwargs = _pipeline_kwargs_from_cfg(pipeline_cfg)

    # Progress over blocks (not per sig) to avoid interleaved stdout from workers
    blocks = list(_chunks(todo_sig_ids, chunk_size))


    # Parent will buffer rows and flush periodically
    buffer_rows: list[dict[str, Any]] = []
    flush_every_blocks = 1  # flush after each completed block; safe for long runs

    # submit all jobs first (fast)
    ctx = mp.get_context("spawn")  # safest on macOS
    with ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as ex:
        futures = [
            ex.submit(
                _process_block_worker,
                zarr_path=zarr_path,
                sig_ids_block=block,
                pipeline_kwargs=pipeline_kwargs,
                metric=metric,
                threshold=float(threshold),
            )
            for block in blocks
        ]

        # NOW show progress as blocks complete (this reflects real work)
        # for fut in tqdm(as_completed(futures), total=len(futures), desc="batch scan (done blocks)")
        for fut in as_completed(futures):
            rows, hit_ids, block_errors = fut.result()

            n_done += len(rows)
            n_errors += int(block_errors)
            n_hits += len(hit_ids)

            buffer_rows.extend(rows)

            # flush each completed block (safe for long runs)
            if buffer_rows:
                _append_summary_rows_csv(summary_csv_path, buffer_rows, write_header=write_header)
                write_header = False
                buffer_rows.clear()

            print(
                f"[batch] processed={n_done + n_skipped}/{n_total} "
                f"(done={n_done}, skipped={n_skipped}, errors={n_errors}, hits={n_hits})"
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


def _compute_threshold_score(
    metric: str,
    sp_pos: Optional[float],
    sp_neg: Optional[float],
    sp_abs: Optional[float],
) -> Optional[float]:
    metric = str(metric).strip().lower()

    if metric == "pos_all":
        return sp_pos

    if metric == "neg_all":
        return sp_neg

    if metric == "abs_all":
        return sp_abs

    if metric == "min_all":
        vals = [x for x in (sp_pos, sp_neg, sp_abs) if x is not None]
        return min(vals) if vals else None

    raise ValueError(f"Unsupported threshold metric: {metric}")


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
