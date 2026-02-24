# `config_example.json` — field-by-field guide

`config_example.json` is a *template* for the repo’s `config.json`.

The real config is split into three top-level sections:

- `pipeline`: shared algorithm params (used by both single-run and batch)
- `single_run`: inputs/outputs for running one signature
- `batch`: inputs/outputs for scanning many signatures

> Paths can be absolute or relative.
> Relative paths are interpreted relative to your current working directory.

---

## 1) `pipeline` (shared algorithm parameters)

### `pipeline.pathway_csv` (string, required)
Path to the pathway mapping CSV.

Supported formats:
- **Long format**: columns `pathway`, `gene` (one row per gene)
- **Summary format**: columns `pathway` and one of `genes` / `gene_ids` (list-like strings or lists)

### `pipeline.direction` (string, default: `"both"`)
One of:
- `"pos"` — use only positive scores
- `"neg"` — use only negative scores
- `"both"` — run both and return separate results per direction

### `pipeline.alpha` (float, default: `0.2`)
alpha parameter for alpha-RRA.

Constraints: `0 < alpha < 1`.

### `pipeline.n_perm` (int, default: `200`)
Number of permutations used by alpha-RRA to estimate p-values.

- `0` is allowed and is useful for fast debugging runs.

### `pipeline.seed` (int or null, default: `0`)
Random seed for alpha-RRA permutations.

- Use `null` to disable seeding.

### `pipeline.X` (int, default: `1`)
XL-mHG parameter: minimum number of hits required at a cutoff.

### `pipeline.L` (int, default: `100`)
XL-mHG parameter: maximum cutoff position (considers top `L` items).

### `pipeline.spearman_plot_zoom_top_fraction` (float, default: `0.1`)
Used for the *zoomed* rank agreement plot (top fraction shown).

Constraints: `0 < f <= 1`.

### `pipeline.show_progress` (bool, default: `true`)
Controls the per-pathway progress bar in single-run pipeline execution.

---

## 2) `single_run` (one signature)

### `single_run.signature_path` (string, required)
Path to a signature file (CSV/TSV) containing gene IDs + a numeric score.

Accepted column layouts:

1) Standard:
- `gene`, `score` (case-insensitive)

2) CLUE-style:
- `gene_id`, plus exactly one score column

> Numeric gene IDs are supported (they’re treated as strings internally).

### Output paths
These are **file paths**, not directories. Set any to `null` to disable writing.

- `single_run.output_csv`: results CSV path
- `single_run.output_spearman_json`: Spearman dict JSON path
- `single_run.output_spearman_plot`: PNG path for the full rank-agreement plot
- `single_run.output_spearman_plot_zoom`: PNG path for the zoomed plot

### `single_run.output_dir` (string, optional)
Not required by the core CLI. Purely a convenience key you can use in your own wrappers.

---

## 3) `batch` (scan many signatures)

Batch mode reads many signatures from a store, runs the pipeline, and writes a summary.

### `batch.signature_store` (object, required)
Currently supported format:

- `type`: must be `"zarr"`
- `path`: path to the `.zarr` directory

Expected datasets in the Zarr root group:
- `scores`: 2D array, shape `(n_genes, n_signatures)`
- `gene_ids`: 1D array, length `n_genes`
- `sig_ids`: 1D array, length `n_signatures`

### `batch.sig_ids_path` (string, required)
Path to a text file containing signature IDs to process.

- One ID per line is typical.
- Duplicates are removed while preserving order.

### Thresholding

After running per-signature enrichment, the batch summary stores Spearman correlations for each direction.
A signature is a **hit** if `threshold_score < threshold`.

- `batch.threshold` (float)
- `batch.threshold_metric` (string)
  Supported values:
  - `"pos_all"`: threshold score = Spearman(pos)
  - `"neg_all"`: threshold score = Spearman(neg)
  - `"min_all"`: threshold score = min(Spearman(pos), Spearman(neg)) ignoring missing

- `batch.top_fractions` (array of float)
  Reserved for future per-hit artifacts/metrics.

### `batch.outputs` (object)

- `summary_path` (string, required)
  Where to write the batch summary.
  If it ends with `.parquet`, the runner will also maintain an adjacent `.csv` for appending.

- `hits_dir` (string or null)
  Directory for per-hit artifacts (one subfolder per signature).

- `save_artifacts_for_hits` (bool)
  If `true`, saves artifacts (CSV results + JSON spearman + params) for signatures marked as hits.

### `batch.execution` (object)

- `resume` (bool)
  If `true`, signatures already present in the existing summary are skipped.

- `max_signatures` (int or null)
  Optional cap for debugging.

- `n_workers` (int)
  Number of worker processes.

  Notes:
  - Some implementations only support `1` worker; others support multiprocessing.

- `chunk_size` (int)
  Used to flush summary rows to disk periodically, and/or split blocks in multiprocessing.

- `show_progress` (bool)
  Controls the batch progress bars.

---

## Common pitfalls

1) **Your current `config.json` uses directories for `single_run.output_csv` / `output_spearman_json` / plots**

Those should be **full file paths**, like:

- `.../results.csv`
- `.../spearman.json`
- `.../rank_agreement.png`

2) **Batch run finishes immediately with `skipped=...`**

That usually means `resume=true` and those sig IDs already exist in the summary file.

3) **No progress bar appears**

- If `batch.execution.show_progress=false` or `--no-progress` is set, tqdm is disabled.
- If `tqdm` isn’t installed, progress bars fall back silently.
