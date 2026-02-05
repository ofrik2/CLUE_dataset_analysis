from __future__ import annotations
import pandas as pd

from clue_pathway_enrichment.io.load_signature import load_signature_table
from clue_pathway_enrichment.io.load_pathways import load_pathway_mapping_csv
from clue_pathway_enrichment.preprocessing.create_ranked_list import pathway_to_binary_vector
from clue_pathway_enrichment.preprocessing.split_rank_signature import split_and_rank_signature
from clue_pathway_enrichment.methods.alpha_rra_wrapper import run_alpha_rra
from clue_pathway_enrichment.methods.xlmhg_wrapper import run_xlmhg
from clue_pathway_enrichment.analysis.ranking import add_ranks
from clue_pathway_enrichment.analysis.correlations import spearman_between_ranks

def run(
    signature_path: str,
    pathway_csv: str,
    direction: str = "both",   # "pos" | "neg" | "both"
    alpha: float = 0.2,
    n_perm: int = 200,
    seed: int | None = 0,
    X: int = 1,
    L: int = 100,
) -> tuple[pd.DataFrame, dict[str, float]]:
    sig = load_signature_table(signature_path)
    pathways = load_pathway_mapping_csv(pathway_csv)

    pos_ranked, neg_ranked = split_and_rank_signature(sig)

    ranked_by_dir = {"pos": pos_ranked, "neg": neg_ranked}
    enabled = ["pos", "neg"] if direction == "both" else [direction]
    for d in enabled:
        if d not in ranked_by_dir:
            raise ValueError(f"direction must be pos|neg|both, got: {direction}")

    rows = []
    for d in enabled:
        ranked = ranked_by_dir[d]
        if ranked.empty:
            continue

        for pname, geneset in pathways.items():
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
    spearman_by_dir: dict[str, float] = {}
    for d in enabled:
        sub = res[res["direction"] == d].copy()
        if sub.empty:
            continue
        sub = add_ranks(sub, "alpha_rra_p", "rank_alpha_rra")
        sub = add_ranks(sub, "xlmhg_p", "rank_xlmhg")
        spearman_by_dir[d] = spearman_between_ranks(sub, "rank_alpha_rra", "rank_xlmhg")
        out.append(sub)

    res2 = pd.concat(out, ignore_index=True) if out else pd.DataFrame()
    return res2, spearman_by_dir
