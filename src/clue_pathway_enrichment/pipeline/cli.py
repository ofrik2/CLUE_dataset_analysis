from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Optional

from clue_pathway_enrichment.pipeline.config import load_pipeline_config
from clue_pathway_enrichment.pipeline.run_pipeline import run


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
    return p


def main(argv: Optional[list[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    cfg = load_pipeline_config(args.config)

    # merge overrides
    direction = args.direction if args.direction is not None else cfg.direction
    alpha = args.alpha if args.alpha is not None else cfg.alpha
    n_perm = args.n_perm if args.n_perm is not None else cfg.n_perm
    seed = args.seed if args.seed is not None else cfg.seed
    X = args.X if args.X is not None else cfg.X
    L = args.L if args.L is not None else cfg.L

    output_csv = args.output_csv if args.output_csv is not None else cfg.output_csv
    output_spearman = args.output_spearman_json if args.output_spearman_json is not None else cfg.output_spearman_json

    res_df, spearman_by_dir = run(
        signature_path=cfg.signature_path,
        pathway_csv=cfg.pathway_csv,
        direction=direction,
        alpha=alpha,
        n_perm=n_perm,
        seed=seed,
        X=X,
        L=L,
    )

    if output_csv:
        Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
        res_df.to_csv(output_csv, index=False)

    if output_spearman:
        Path(output_spearman).parent.mkdir(parents=True, exist_ok=True)
        Path(output_spearman).write_text(json.dumps(spearman_by_dir, indent=2, sort_keys=True))

    # also print a tiny summary for interactive use
    print(f"rows={len(res_df)}")
    if spearman_by_dir:
        print(f"spearman_by_dir={spearman_by_dir}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
