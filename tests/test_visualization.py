from __future__ import annotations

from pathlib import Path

import pandas as pd

from clue_pathway_enrichment.analysis.visualization import save_rank_agreement_plot


def test_save_rank_agreement_plot_writes_file(tmp_path: Path):
    df = pd.DataFrame(
        {
            "direction": ["pos", "pos", "pos"],
            "rank_alpha_rra": [1, 2, 3],
            "rank_xlmhg": [1, 3, 2],
        }
    )

    out = tmp_path / "spearman.png"
    save_rank_agreement_plot(df, direction="pos", spearman_rho=0.5, out_path=str(out))

    assert out.exists()
    assert out.stat().st_size > 0


def test_save_rank_agreement_plot_zoom_writes_file(tmp_path: Path):
    df = pd.DataFrame(
        {
            "direction": ["pos"] * 20,
            "rank_alpha_rra": list(range(1, 21)),
            "rank_xlmhg": list(range(20, 0, -1)),
        }
    )

    out = tmp_path / "spearman.zoom.png"
    save_rank_agreement_plot(df, direction="pos", out_path=str(out), zoom_top_fraction=0.1)

    assert out.exists()
    assert out.stat().st_size > 0
