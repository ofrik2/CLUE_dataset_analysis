from __future__ import annotations
from typing import Optional
import pandas as pd

from clue_pathway_enrichment.io.load_signature import load_signature_table
from clue_pathway_enrichment.io.load_signature import standardize_signature_df
from clue_pathway_enrichment.io.load_pathways import load_pathway_mapping_csv
from clue_pathway_enrichment.preprocessing.create_ranked_list import pathway_to_binary_vector
from clue_pathway_enrichment.preprocessing.split_rank_signature import split_and_rank_signature
from clue_pathway_enrichment.methods.alpha_rra_wrapper import run_alpha_rra
from clue_pathway_enrichment.methods.xlmhg_wrapper import run_xlmhg
from clue_pathway_enrichment.analysis.ranking import add_ranks
from clue_pathway_enrichment.analysis.correlations import spearman_between_ranks
from clue_pathway_enrichment.analysis.visualization import save_rank_agreement_plot


def _iter_with_progress(it, *, total: int, desc: str, enabled: bool):
    """Yield from an iterable, optionally wrapped in a tqdm progress bar."""

    if not enabled:
        yield from it
        return

    try:
        from tqdm import tqdm  # type: ignore
    except ModuleNotFoundError:  # pragma: no cover
        # Fall back silently if tqdm isn't installed.
        yield from it
        return

    yield from tqdm(it, total=total, desc=desc, unit="pw", mininterval=0.5)



def run(
    signature_path: str | pd.DataFrame,
    pathway_csv: str,
    direction: str = "both",   # "pos" | "neg" | "both"
    alpha: float = 0.2,
    n_perm: int = 200,
    seed: Optional[int] = 0,
    X: int = 1,
    L: int = 100,
    # output_spearman_plot: Optional[str] = None,  # e.g. "out/rank_agreement_{direction}.png"
    output_spearman_plot_zoom: Optional[str] = None,  # e.g. "out/rank_agreement_top10pct_{direction}.png"
    spearman_plot_zoom_top_fraction: float = 0.1,
    show_progress: bool = True,
) -> tuple[pd.DataFrame, dict[str, float] | dict[str, dict[str, float]]]:
    if isinstance(signature_path, str):
        sig = load_signature_table(signature_path)
    elif isinstance(signature_path, pd.DataFrame):
        sig = standardize_signature_df(signature_path)
    else:
        raise TypeError(f"signature_path must be a str path or a pandas DataFrame, got: {type(signature_path)!r}")

    pathways = load_pathway_mapping_csv(pathway_csv)

    pos_ranked, neg_ranked = split_and_rank_signature(sig)
    ranked_by_dir = {"pos": pos_ranked, "neg": neg_ranked}

    enabled = ["pos", "neg"] if direction == "both" else [direction]
    for d in enabled:
        if d not in ranked_by_dir:
            raise ValueError(f"direction must be pos|neg|both, got: {direction}")

    rows = []
    pathway_items = list(pathways.items())
    total_pathways = len(pathway_items)

    for d in enabled:
        ranked = ranked_by_dir[d]
        if ranked.empty:
            continue

        desc = f"scoring pathways ({d})"
        # iterating pathways
        for pname, geneset in _iter_with_progress(
            pathway_items,
            total=total_pathways,
            desc=desc,
            enabled=show_progress,
        ):
            v = pathway_to_binary_vector(ranked, geneset)
            k = int(v.sum())
            if k == 0:
                continue

            rra_stat, rra_p = run_alpha_rra(v, alpha=alpha, n_perm=n_perm, seed=seed)
            xlmhg_stat, xlmhg_p = run_xlmhg(v, X=X, L=L)

            rows.append({
                "direction": d,
                "pathway": pname,
                "K_hits": k,
                "N": len(v),
                "alpha_rra_stat": rra_stat,
                "alpha_rra_p": rra_p,
                "xlmhg_stat": xlmhg_stat,
                "xlmhg_p": xlmhg_p,
            })

    res = pd.DataFrame(rows)
    if res.empty:
        return res, {}

    # rank within each direction separately
    out = []
    spearman_by_dir: dict[str, float] | dict[str, dict[str, float]] = {}

    for d in enabled:
        sub = res[res["direction"] == d].copy()
        if sub.empty:
            continue
        sub = add_ranks(sub, "alpha_rra_p", "rank_alpha_rra")
        sub = add_ranks(sub, "xlmhg_p", "rank_xlmhg")
        full_rho = spearman_between_ranks(sub, "rank_alpha_rra", "rank_xlmhg")

        # If we end up computing zoom rho, we’ll promote the structure to nested dicts.
        zoom_rho: Optional[float] = None

        # if output_spearman_plot:
        #     plot_path = output_spearman_plot.format(direction=d)
        #     save_rank_agreement_plot(
        #         sub,
        #         direction=d,
        #         spearman_rho=full_rho,
        #         out_path=plot_path,
        #     )
        #
        if spearman_plot_zoom_top_fraction:
            # plot_path = output_spearman_plot_zoom.format(direction=d)

            # Recompute Spearman on the same subset that will be shown in the zoom plot:
            zf = float(spearman_plot_zoom_top_fraction)
            if not (0.0 < zf <= 1.0):
                raise ValueError(f"spearman_plot_zoom_top_fraction must be in (0, 1], got: {zf}")

            n = int(len(sub))
            if n >= 2:
                zoom_n = max(1, int((n * zf + 0.999999999)))  # ceil
                tmp = sub[["rank_alpha_rra", "rank_xlmhg"]].copy()
                tmp["__min_rank__"] = tmp[["rank_alpha_rra", "rank_xlmhg"]].min(axis=1)
                zoom_sub = sub.loc[tmp["__min_rank__"] <= zoom_n]
                zoom_rho = (
                    spearman_between_ranks(zoom_sub, "rank_alpha_rra", "rank_xlmhg")
                    if len(zoom_sub) >= 2
                    else float("nan")
                )
            else:
                zoom_rho = float("nan")

            # save_rank_agreement_plot(
            #     sub,
            #     direction=d,
            #     spearman_rho=zoom_rho,
            #     out_path=plot_path,
            #     zoom_top_fraction=spearman_plot_zoom_top_fraction,
            # )

        # Store return values
        if zoom_rho is None:
            # keep old behavior: dict[str, float]
            if isinstance(spearman_by_dir, dict) and (not spearman_by_dir or isinstance(next(iter(spearman_by_dir.values()), 0.0), float)):
                spearman_by_dir[d] = full_rho  # type: ignore[index]
            else:
                # already promoted
                spearman_by_dir[d] = {"full": full_rho}  # type: ignore[index]
        else:
            # promote to nested dict if needed
            if spearman_by_dir and isinstance(next(iter(spearman_by_dir.values())), float):  # type: ignore[arg-type]
                spearman_by_dir = {k: {"full": v} for k, v in spearman_by_dir.items()}  # type: ignore[assignment]
            spearman_by_dir[d] = {"full": full_rho, "top_k": zoom_rho}  # type: ignore[index]

        out.append(sub)

    res2 = pd.concat(out, ignore_index=True) if out else pd.DataFrame()
    return res2, spearman_by_dir



def run_from_config(
    config_path: str,
    *,
    direction: Optional[str] = None,
    alpha: Optional[float] = None,
    n_perm: Optional[int] = None,
    seed: Optional[int] = None,
    X: Optional[int] = None,
    L: Optional[int] = None,
    output_spearman_plot: Optional[str] = None,
    output_spearman_plot_zoom: Optional[str] = None,
    show_progress: Optional[bool] = None,
) -> tuple[pd.DataFrame, dict[str, float]]:
    """Load a config (JSON preferred; TOML supported) and run the pipeline.

    This is a thin helper around :func:`run`.
    """

    from clue_pathway_enrichment.pipeline.config import load_pipeline_config

    cfg = load_pipeline_config(config_path)

    return run(
        signature_path=cfg.signature_path,
        pathway_csv=cfg.pathway_csv,
        direction=direction if direction is not None else cfg.direction,
        alpha=alpha if alpha is not None else cfg.alpha,
        n_perm=n_perm if n_perm is not None else cfg.n_perm,
        seed=seed if seed is not None else cfg.seed,
        X=X if X is not None else cfg.X,
        L=L if L is not None else cfg.L,
        output_spearman_plot=output_spearman_plot if output_spearman_plot is not None else getattr(cfg, "output_spearman_plot", None),
        output_spearman_plot_zoom=output_spearman_plot_zoom if output_spearman_plot_zoom is not None else getattr(cfg, "output_spearman_plot_zoom", None),
        spearman_plot_zoom_top_fraction=float(getattr(cfg, "spearman_plot_zoom_top_fraction", 0.1)),
        show_progress=(show_progress if show_progress is not None else bool(getattr(cfg, "show_progress", True))),
    )
