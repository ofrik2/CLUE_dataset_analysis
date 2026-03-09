from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Optional

from clue_pathway_enrichment.pipeline.config import load_pipeline_config, load_single_run_config
from clue_pathway_enrichment.pipeline.run_pipeline import run
from clue_pathway_enrichment.analysis.visualization import save_rank_agreement_plot



def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run CLUE pathway enrichment from a JSON config file")
    p.add_argument("--config", required=True, help="Path to config (.json preferred; .toml supported)")

    # Optional overrides
    p.add_argument("--direction", choices=["pos", "neg", "both"], default=None)
    p.add_argument("--alpha", type=float, default=None)
    p.add_argument("--n-perm", type=int, default=None)
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--X", type=int, default=None)
    p.add_argument("--L", type=int, default=None)

    p.add_argument("--output-csv", default=None, help="Write result table to CSV (overrides config)")
    p.add_argument("--output-spearman-json", default=None, help="Write Spearman dict to JSON (overrides config)")
    p.add_argument(
        "--output-spearman-plot",
        default=None,
        help=(
            "Write a rank-agreement plot (rank_alpha_rra vs rank_xlmhg) as an image (overrides config). "
            "If direction=both, writes one file per direction by appending '.pos' / '.neg' before the extension."
        ),
    )
    p.add_argument(
        "--output-spearman-plot-zoom",
        default=None,
        help=(
            "Write an additional zoomed-in rank-agreement plot focusing on the top fraction of ranks (overrides config). "
            "If direction=both, writes one file per direction by appending '.pos' / '.neg' before the extension."
        ),
    )
    p.add_argument(
        "--spearman-plot-zoom-top-fraction",
        type=float,
        default=None,
        help="Top fraction for zoomed plot (0 < f <= 1). Default: 0.1 (top 10%).",
    )

    p.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable the per-pathway progress bar / ETA.",
    )
    return p


def _with_dir_suffix(path: str, direction: str) -> str:
    p = Path(path)
    if p.suffix:
        return str(p.with_name(f"{p.stem}.{direction}{p.suffix}"))
    return str(p.with_name(f"{p.name}.{direction}"))


def main(argv: Optional[list[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv) #parses config

    pipe_cfg = load_pipeline_config(args.config)
    single_cfg = load_single_run_config(args.config)

    # merge overrides
    direction = args.direction if args.direction is not None else pipe_cfg.direction
    alpha = args.alpha if args.alpha is not None else pipe_cfg.alpha
    n_perm = args.n_perm if args.n_perm is not None else pipe_cfg.n_perm
    seed = args.seed if args.seed is not None else pipe_cfg.seed
    X = args.X if args.X is not None else pipe_cfg.X
    L = args.L if args.L is not None else pipe_cfg.L

    output_csv = args.output_csv if args.output_csv is not None else single_cfg.output_csv
    output_spearman = args.output_spearman_json if args.output_spearman_json is not None else single_cfg.output_spearman_json
    output_spearman_plot = (
        args.output_spearman_plot if args.output_spearman_plot is not None else single_cfg.output_spearman_plot
    )
    output_spearman_plot_zoom = (
        args.output_spearman_plot_zoom
        if args.output_spearman_plot_zoom is not None
        else getattr(single_cfg, "output_spearman_plot_zoom", None)
    )
    spearman_plot_zoom_top_fraction = (
        args.spearman_plot_zoom_top_fraction
        if args.spearman_plot_zoom_top_fraction is not None
        else float(getattr(single_cfg, "spearman_plot_zoom_top_fraction", 0.1))
    )

    signature_ranking = pipe_cfg.signature_ranking

    show_progress = (not args.no_progress) and bool(getattr(pipe_cfg, "show_progress", True))

    res_df, spearman_by_dir = run(
        signature_path=single_cfg.signature_path,
        pathway_csv=pipe_cfg.pathway_csv,
        direction=direction,
        signature_ranking=signature_ranking,
        alpha=alpha,
        n_perm=n_perm,
        seed=seed,
        X=X,
        L=L,
        # output_spearman_plot=output_spearman_plot,  # <-- wire through
        # output_spearman_plot_zoom=output_spearman_plot_zoom,
        spearman_plot_zoom_top_fraction=spearman_plot_zoom_top_fraction,
        show_progress=show_progress,
    )

    if output_csv:
        Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
        res_df.to_csv(output_csv, index=False)

    if output_spearman:
        Path(output_spearman).parent.mkdir(parents=True, exist_ok=True)
        Path(output_spearman).write_text(json.dumps(spearman_by_dir, indent=2, sort_keys=True))

    # (Re-)create plots from the final result DF so CLI can suffix per-direction filenames.
    # (Re-)create plots: create one file per direction suffix (.pos/.neg)
    if output_spearman_plot or output_spearman_plot_zoom:
        # Use the actual directions present in the result dataframe.
        # This supports both:
        # - signed_split mode: ["pos", "neg"]
        # - abs mode: ["abs"]
        dirs = sorted(res_df["direction"].dropna().unique().tolist())

        for d in dirs:
            ddf = res_df[res_df["direction"] == d].copy()

            # Skip completely empty subsets just in case
            if ddf.empty:
                continue

            # spearman_by_dir[d] can be either:
            # - float (only full), or
            # - {"full": float, "top_k": float} when zoom is enabled
            val = spearman_by_dir.get(d)
            if isinstance(val, dict):
                rho_full = val.get("full")
                rho_zoom = val.get("top_k")
            else:
                rho_full = val
                rho_zoom = None

            # Full plot
            if output_spearman_plot:
                save_rank_agreement_plot(
                    ddf,
                    direction=d,
                    spearman_rho=rho_full,
                    out_path=_with_dir_suffix(output_spearman_plot, d),
                    title=f"Rank agreement ({d})"
                          + (f" | Spearman ρ={rho_full:.3f}" if isinstance(rho_full, float) else ""),
                )

            # Zoom plot
            if output_spearman_plot_zoom and spearman_plot_zoom_top_fraction:
                save_rank_agreement_plot(
                    ddf,
                    direction=d,
                    spearman_rho=rho_zoom,
                    out_path=_with_dir_suffix(output_spearman_plot_zoom, d),
                    title=f"Rank agreement ({d}) — top {int(spearman_plot_zoom_top_fraction * 100)}%"
                          + (f" | Spearman ρ={rho_zoom:.3f}" if isinstance(rho_zoom, float) else ""),
                    zoom_top_fraction=spearman_plot_zoom_top_fraction,
                )

    # also print a tiny summary for interactive use
    print(f"rows={len(res_df)}")
    if spearman_by_dir:
        print(f"spearman_by_dir={spearman_by_dir}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
