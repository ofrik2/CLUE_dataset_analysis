from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd


def save_rank_agreement_plot(
    df: pd.DataFrame,
    *,
    direction: str,
    rank1: str = "rank_alpha_rra",
    rank2: str = "rank_xlmhg",
    spearman_rho: Optional[float] = None,
    out_path: str,
    title: Optional[str] = None,
    max_points: int = 5000,
    point_size: float = 12.0,
    zoom_top_fraction: Optional[float] = None,
) -> None:
    """Save a rank-vs-rank scatter plot for a single direction.

    Parameters
    ----------
    df:
        Results DataFrame containing rank columns.
    direction:
        Direction label used for filtering and in the default title.
    rank1, rank2:
        Column names to plot on x/y axes.
    spearman_rho:
        If provided, included in the title.
    out_path:
        Where to save the plot. Extension decides format (e.g. .png, .pdf).
    title:
        Plot title override.
    max_points:
        Upper bound on plotted points (plots the top rows deterministically).
    zoom_top_fraction:
        If set (e.g. 0.1), plot only the top fraction of pathways by rank.

        The subset is defined as rows where *either* method places the pathway in
        the top fraction, i.e. `min(rank1, rank2) <= ceil(n * zoom_top_fraction)`.
        This focuses the plot on the most relevant hits while still reflecting
        disagreement between methods.
    """

    try:
        import matplotlib

        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt  # type: ignore
    except ModuleNotFoundError as e:  # pragma: no cover
        raise ModuleNotFoundError(
            "Plotting requires matplotlib. Install it (e.g. `pip install matplotlib`) or disable output_spearman_plot."
        ) from e

    sub = df[df["direction"] == direction].copy() if "direction" in df.columns else df.copy()
    if sub.empty:
        # Still create an empty plot so pipelines can rely on the file existing.
        sub = pd.DataFrame({rank1: [], rank2: []})

    if rank1 not in sub.columns or rank2 not in sub.columns:
        raise KeyError(f"Missing rank columns for plotting: {rank1!r}, {rank2!r}. Available: {list(sub.columns)}")

    sub = sub[[rank1, rank2]].dropna()

    zoom_n: Optional[int] = None
    if zoom_top_fraction is not None:
        if zoom_top_fraction <= 0 or zoom_top_fraction > 1:
            raise ValueError(f"zoom_top_fraction must be in (0, 1], got: {zoom_top_fraction}")
        if len(sub):
            # ceil, with a minimum of 1 point.
            zoom_n = max(1, int((len(sub) * float(zoom_top_fraction) + 0.999999999)))
            tmp = sub.copy()
            tmp["__min_rank__"] = tmp[[rank1, rank2]].min(axis=1)
            sub = tmp[tmp["__min_rank__"] <= zoom_n].drop(columns=["__min_rank__"])

    if len(sub) > max_points:
        sub = sub.iloc[:max_points]

    x = sub[rank1].astype(float).to_numpy()
    y = sub[rank2].astype(float).to_numpy()

    if title is None:
        bits = [f"Rank agreement ({direction})"]
        if zoom_top_fraction is not None:
            bits.append(
                f"zoom top={zoom_top_fraction:.0%}"
                if zoom_n is None
                else f"zoom top={zoom_top_fraction:.0%} (≤{zoom_n})"
            )

        # --- FIX: allow spearman_rho to be float OR dict keyed by direction ---
        rho_val = spearman_rho
        if isinstance(rho_val, dict):
            rho_val = rho_val.get(direction)

        if rho_val is not None:
            try:
                bits.append(f"Spearman ρ={float(rho_val):.3f}")
            except (TypeError, ValueError):
                # If it's something weird (e.g., list/object), skip formatting
                pass
        # --- end fix ---

        bits.append(f"n={len(sub)}")
        title = " | ".join(bits)


    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(x, y, s=point_size, alpha=0.7, edgecolors="none")
    ax.set_title(title)
    ax.set_xlabel(rank1)
    ax.set_ylabel(rank2)

    # Make it easier to see agreement.
    if len(x) and len(y):
        lo = float(min(x.min(), y.min()))
        hi = float(max(x.max(), y.max()))
        ax.plot([lo, hi], [lo, hi], linestyle="--", linewidth=1, color="black", alpha=0.4)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)

    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
