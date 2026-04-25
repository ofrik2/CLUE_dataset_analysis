import pandas as pd

from clue_pathway_enrichment.pipeline.run_pipeline import run


def test_run_pipeline_forwards_alpha_rra_mp_params_and_uses_pathway_specific_seed(monkeypatch):
    calls = []

    # Fake signature loader / standardizer
    sig_df = pd.DataFrame({
        "gene": ["g1", "g2", "g3", "g4"],
        "score": [4.0, 3.0, 2.0, 1.0],
    })

    # ranked_by_dir expected by run()
    ranked_pos = pd.DataFrame({"gene": ["g1", "g2", "g3", "g4"], "score": [4.0, 3.0, 2.0, 1.0]})
    ranked_neg = pd.DataFrame({"gene": ["g4", "g3", "g2", "g1"], "score": [-4.0, -3.0, -2.0, -1.0]})

    pathways = {
        "p1": {"g1", "g2"},
        "p2": {"g3"},
    }

    def fake_load_signature_table(path):
        return sig_df

    def fake_load_pathway_mapping_csv(path):
        return pathways

    def fake_rank_signature(sig, mode):
        return {"pos": ranked_pos, "neg": ranked_neg}

    def fake_pathway_to_binary_vector(ranked, geneset):
        genes = list(ranked["gene"])
        return pd.Series([1 if g in geneset else 0 for g in genes]).to_numpy()

    def fake_run_alpha_rra(v, alpha, n_perm, seed, n_workers, mp_chunk_size):
        calls.append({
            "seed": seed,
            "n_workers": n_workers,
            "mp_chunk_size": mp_chunk_size,
            "n_perm": n_perm,
            "alpha": alpha,
            "k": int(v.sum()),
            "N": len(v),
        })
        return 0.1, 0.2

    def fake_run_xlmhg(v, X, L):
        return 0.3, 0.4

    def fake_add_ranks(df, p_col, rank_col):
        df = df.copy()
        df[rank_col] = range(1, len(df) + 1)
        return df

    def fake_spearman_between_ranks(df, c1, c2):
        return 1.0

    monkeypatch.setattr("clue_pathway_enrichment.pipeline.run_pipeline.load_signature_table", fake_load_signature_table)
    monkeypatch.setattr("clue_pathway_enrichment.pipeline.run_pipeline.load_pathway_mapping_csv", fake_load_pathway_mapping_csv)
    monkeypatch.setattr("clue_pathway_enrichment.pipeline.run_pipeline.rank_signature", fake_rank_signature)
    monkeypatch.setattr("clue_pathway_enrichment.pipeline.run_pipeline.pathway_to_binary_vector", fake_pathway_to_binary_vector)
    monkeypatch.setattr("clue_pathway_enrichment.pipeline.run_pipeline.run_alpha_rra", fake_run_alpha_rra)
    monkeypatch.setattr("clue_pathway_enrichment.pipeline.run_pipeline.run_xlmhg", fake_run_xlmhg)
    monkeypatch.setattr("clue_pathway_enrichment.pipeline.run_pipeline.add_ranks", fake_add_ranks)
    monkeypatch.setattr("clue_pathway_enrichment.pipeline.run_pipeline.spearman_between_ranks", fake_spearman_between_ranks)

    res_df, spearman_by_dir = run(
        signature_path="dummy_sig.csv",
        pathway_csv="dummy_pathways.csv",
        direction="both",
        signature_ranking="signed_split",
        alpha=0.2,
        n_perm=5000,
        seed=100,
        X=1,
        L=None,
        show_progress=False,
        alpha_rra_n_workers=4,
        alpha_rra_mp_chunk_size=1000,
    )

    # two pathways x two directions = four calls
    assert len(calls) == 4

    # New params are forwarded
    assert all(c["n_workers"] == 4 for c in calls)
    assert all(c["mp_chunk_size"] == 1000 for c in calls)
    assert all(c["n_perm"] == 5000 for c in calls)

    # Seeds should be pathway-specific and direction-specific
    # pos: 100, 101 ; neg: 102, 103
    assert [c["seed"] for c in calls] == [100, 101, 102, 103]

    assert not res_df.empty
    assert set(res_df["direction"]) == {"pos", "neg"}
    assert set(spearman_by_dir.keys()) == {"pos", "neg"}