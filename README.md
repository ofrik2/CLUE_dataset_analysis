# CLUE Pathway Enrichment

This repository contains a small Python package, `clue_pathway_enrichment`, that runs a simple **pathway enrichment** analysis on a CLUE-style gene signature.

At a high level it:

1. Loads a **signature** table: gene IDs + a numeric score.
2. Splits it into **positive** and **negative** directions (or both).
3. For each pathway gene set, builds a **binary hit vector** over the ranked genes.
4. Runs two enrichment-style methods:
   - **alpha-RRA** (permutation p-value)
   - **XL-mHG**
5. Ranks pathways by p-value and reports rank correlations.

The code is intentionally lightweight and meant as a starting point you can extend.

---

## Project layout

- `config.json`: editable run configuration (committed to the repo)
- `src/clue_pathway_enrichment/`: package code
  - `io/`: reading signature/pathway files
  - `preprocessing/`: ranking / splitting / hit vectors
  - `methods/`: alpha-RRA and XL-mHG wrappers
  - `analysis/`: ranking utilities and correlations (+ plotting)
  - `pipeline/`: orchestration (`run_pipeline.py`), config/CLI
- `tests/`: unit/integration tests

---

## Quick start (single run)

From the repo root:

```bash
PYTHONPATH=$PWD/src python3 -m clue_pathway_enrichment.pipeline.cli --config config.json
```

If you install the package (e.g. `python -m pip install -e .`), you also get a console script:

```bash
clue-pathway-enrichment --config config.json
```

You can also run the module after an editable install:

```bash
python -m clue_pathway_enrichment.pipeline.cli --config config.json
```

---

## Configuration (`config.json`)

### Supported formats

- **JSON** is preferred.
- **TOML** is supported as a backward-compatible fallback.

Config selection rules used by the loaders in `pipeline/config.py`:
- `*.json` => JSON
- `*.toml` => TOML
- otherwise: it tries JSON first, then TOML

### Top-level structure

The repo’s config is structured into three top-level sections:

- `pipeline`: shared algorithm parameters (used by both single-run and batch).
- `single_run`: single-signature inputs/outputs (used by `pipeline/cli.py`).
- `batch`: batch execution inputs/outputs (used by the batch runner).

---

## Single run

### `pipeline` section (shared params)

Fields:

- `pathway_csv` (string, **required**)
  Path to your pathway mapping file (CSV).

- `direction` (string, default: `"both"`)
  One of: `"pos"`, `"neg"`, `"both"`.

- `alpha` (float, default: `0.2`)
  alpha-RRA parameter. Must be in `(0, 1)`.

- `n_perm` (int, default: `200`)
  Number of permutations for alpha-RRA p-value estimation. Use `0` for a fast “no permutations” run.

- `seed` (int or null, default: `0`)
  Random seed for alpha-RRA permutations. Use `null` for non-deterministic.

- `X` (int, default: `1`)
  XL-mHG parameter.

- `L` (int, default: `100`)
  XL-mHG parameter.

- `spearman_plot_zoom_top_fraction` (float, default: `0.1`)
  Top fraction of ranks to include in the zoom plot (`0 < f <= 1`).

  The zoom subset is defined as pathways where:
  `min(rank_alpha_rra, rank_xlmhg) <= ceil(n * f)`.

- `show_progress` (bool, default: `true`)
  If `true`, show a tqdm progress bar with ETA during the per-pathway scoring loop.

### `single_run` section (inputs/outputs)

Fields:

- `signature_path` (string, **required**)
  Path to your signature file (CSV or TSV).

- `output_dir` (string, optional)
  Optional convenience field (not required by the core pipeline). If you use it,
  you can put outputs under this folder in your own wrapper code.

- `output_csv` (string or null)
  If set, the CLI writes the main results table to CSV.

- `output_spearman_json` (string or null)
  If set, the CLI writes Spearman correlations as JSON.

- `output_spearman_plot` (string or null)
  If set, writes a rank-agreement scatter plot (`rank_alpha_rra` vs `rank_xlmhg`).

- `output_spearman_plot_zoom` (string or null)
  If set, writes an **additional** rank-agreement plot zoomed into the top fraction.

  Notes on filenames:
  - If `direction="both"`, the CLI writes one file per direction by inserting `.pos` / `.neg`
    before the file extension.

### Spearman output structure

The return value `spearman_by_dir` from `pipeline/run_pipeline.run()` can have two shapes:

- If only the full-list Spearman is computed:
  - `{ "pos": <rho>, "neg": <rho> }`

- If the zoom plot is requested and the zoom Spearman is computed on the zoom subset:
  - `{ "pos": {"full": <rho>, "top_k": <rho_topk>}, "neg": {"full": <rho>, "top_k": <rho_topk>} }`

`top_k` is the Spearman correlation recomputed on the same subset that is shown in the zoom plot.

---

## Input file formats

> Important: the pipeline supports gene **names** (symbols) as well as **numeric gene IDs**.
> The key requirement is **consistency**: the identifier format in the signature must match the identifier format used in `pathway_csv`.

### Signature file (`single_run.signature_path`)

Loaded by `clue_pathway_enrichment.io.load_signature.load_signature_table`.

Supported formats:

1) **Standard**: must include columns named `gene` and `score` (case-insensitive)

2) **CLUE-style**: must include column `gene_id` and exactly one other score column, e.g.:

```tsv
gene_id	my_signature
2547	1.2
1956	-0.7
```

The loader:
- trims gene strings
- coerces score to numeric
- drops rows with missing/invalid gene or score

### Pathway mapping file (`pipeline.pathway_csv`)

Loaded by `clue_pathway_enrichment.io.load_pathways.load_pathway_mapping_csv`.

Supported formats:

1) **Long format**: columns `pathway` and `gene` (one row per member gene)

2) **Summary format**: columns `pathway` and one of `genes` or `gene_ids`.

---

## Running

### CLI (single run)

The CLI is implemented in `src/clue_pathway_enrichment/pipeline/cli.py`.

```bash
PYTHONPATH=$PWD/src python3 -m clue_pathway_enrichment.pipeline.cli --config config.json
```

Useful CLI flags:
- `--output-csv`, `--output-spearman-json`, `--output-spearman-plot`, `--output-spearman-plot-zoom`
- `--spearman-plot-zoom-top-fraction`
- `--no-progress`

### Python API

The main entry point is `clue_pathway_enrichment.pipeline.run_pipeline.run`.

```python
import pandas as pd
from clue_pathway_enrichment.pipeline.run_pipeline import run

# You can pass a path...
results_df, spearman = run(
    signature_path="data/signature.tsv",
    pathway_csv="data/pathways.csv",
    direction="both",
)

# ...or an in-memory DataFrame (must match the signature formats described above):
# - columns (gene, score) OR
# - columns (gene_id, <one score col>)
#
# The pipeline will standardize it internally.
sig_df = pd.DataFrame({"gene": ["A", "B"], "score": [1.0, -1.0]})
results_df2, spearman2 = run(signature_path=sig_df, pathway_csv="data/pathways.csv")
```
### CLI (batch run)

Batch execution is implemented in `src/clue_pathway_enrichment/batch/cli.py`.

From the repo root, run without installing:

```bash
PYTHONPATH=$PWD/src python3 -m clue_pathway_enrichment.batch.cli --config config.json
```

After an editable install (`python -m pip install -e .`), run:

```bash
python -m clue_pathway_enrichment.batch.cli --config config.json
```

> Note: there currently isn’t a dedicated `clue-pathway-enrichment-batch` console script in `pyproject.toml`.
> If you want one, we can add it under `[project.scripts]`.

---

## Batch mode

Batch mode is configured under the top-level `batch` key in `config.json`.

What it typically does:

1. Loads signatures from `batch.signature_store` (e.g. a zarr store).
2. Iterates `sig_ids_path`.
3. Runs the same per-signature pathway enrichment as the single-run pipeline.
4. Writes a summary table to `batch.outputs.summary_path`.
5. Optionally writes per-hit artifacts under `batch.outputs.hits_dir` (if enabled).

### (Optional) Resume / limits

Use the config fields:

- `batch.execution.resume` to resume from existing outputs when supported
- `batch.execution.max_signatures` to cap the run during debugging
- `batch.execution.n_workers` and `batch.execution.chunk_size` to control throughput

### Progress bar / ETA

Batch runs show a progress bar with ETA by default (requires `tqdm`).

Control it via:
- `batch.execution.show_progress` (bool) in `config.json`
- `--no-progress` flag in the batch CLI

---

## Output

### Results table (`single_run.output_csv`)

The pipeline returns a pandas `DataFrame` (and optionally writes it to CSV) with one row per tested `(direction, pathway)`.

Columns typically include:
- `direction`: `pos` or `neg`
- `pathway`: pathway name
- `K_hits`: number of hits for that pathway in the ranked list
- `N`: ranked list size
- `alpha_rra_stat`, `alpha_rra_p`
- `xlmhg_stat`, `xlmhg_p`
- `rank_alpha_rra`: rank of `alpha_rra_p` (lower p-value = better rank)
- `rank_xlmhg`: rank of `xlmhg_p`

---

## Testing

```bash
python -m pip install -r requirements.txt
python -m pytest -q
```

---

## Troubleshooting

### Import errors (package not found)

```bash
PYTHONPATH=$PWD/src ...
```

### Config says files don’t exist

`load_pipeline_config()` validates `pipeline.pathway_csv` exists.
`load_single_run_config()` validates `single_run.signature_path` exists (or falls back to `pipeline.signature_path` for legacy configs).

---

## Notes / next steps

Common extensions:
- multiple-testing correction (FDR/q-values)
- richer batch analytics and reporting
- support selecting a score column for CLUE-style signatures with multiple score columns
