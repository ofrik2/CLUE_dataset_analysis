from types import SimpleNamespace

from clue_pathway_enrichment.pipeline import cli


def test_cli_uses_pipe_cfg_for_alpha_rra_when_no_cli_override(monkeypatch):
    captured = {}

    fake_pipe_cfg = SimpleNamespace(
        pathway_csv="pathways.csv",
        direction="both",
        alpha=0.2,
        n_perm=5000,
        seed=11,
        X=1,
        L=None,
        signature_ranking="signed_split",
        show_progress=True,
        alpha_rra_n_workers=4,
        alpha_rra_mp_chunk_size=1000,
    )

    fake_single_cfg = SimpleNamespace(
        signature_path="sig.csv",
        output_csv=None,
        output_spearman_json=None,
        output_spearman_plot=None,
        output_spearman_plot_zoom=None,
    )

    def fake_load_pipeline_config(path):
        return fake_pipe_cfg

    def fake_load_single_run_config(path):
        return fake_single_cfg

    def fake_run(**kwargs):
        captured.update(kwargs)
        return __import__("pandas").DataFrame({"direction": ["pos"]}), {"pos": 1.0}

    monkeypatch.setattr(cli, "load_pipeline_config", fake_load_pipeline_config)
    monkeypatch.setattr(cli, "load_single_run_config", fake_load_single_run_config)
    monkeypatch.setattr(cli, "run", fake_run)

    rc = cli.main(["--config", "dummy.json"])

    assert rc == 0
    assert captured["alpha_rra_n_workers"] == 4
    assert captured["alpha_rra_mp_chunk_size"] == 1000


def test_cli_cli_override_wins_for_alpha_rra(monkeypatch):
    captured = {}

    fake_pipe_cfg = SimpleNamespace(
        pathway_csv="pathways.csv",
        direction="both",
        alpha=0.2,
        n_perm=5000,
        seed=11,
        X=1,
        L=None,
        signature_ranking="signed_split",
        show_progress=True,
        alpha_rra_n_workers=4,
        alpha_rra_mp_chunk_size=1000,
    )

    fake_single_cfg = SimpleNamespace(
        signature_path="sig.csv",
        output_csv=None,
        output_spearman_json=None,
        output_spearman_plot=None,
        output_spearman_plot_zoom=None,
    )

    def fake_load_pipeline_config(path):
        return fake_pipe_cfg

    def fake_load_single_run_config(path):
        return fake_single_cfg

    def fake_run(**kwargs):
        captured.update(kwargs)
        return __import__("pandas").DataFrame({"direction": ["pos"]}), {"pos": 1.0}

    monkeypatch.setattr(cli, "load_pipeline_config", fake_load_pipeline_config)
    monkeypatch.setattr(cli, "load_single_run_config", fake_load_single_run_config)
    monkeypatch.setattr(cli, "run", fake_run)

    rc = cli.main([
        "--config", "dummy.json",
        "--alpha-rra-n-workers", "2",
        "--alpha-rra-mp-chunk-size", "250",
    ])

    assert rc == 0
    assert captured["alpha_rra_n_workers"] == 2
    assert captured["alpha_rra_mp_chunk_size"] == 250