# Pipeline guide

How the analysis is organised, what each step reads and writes, and how to run
it. Run order lives in [`pipeline.yaml`](../pipeline.yaml); paths live in
[`src/config.py`](../src/config.py). If this page and those files disagree,
the files are right.

## Datasets

| Key        | Notebooks            | Day | Condition                 |
|------------|----------------------|-----|---------------------------|
| `D2_DZ`    | `notebooks/d2-dz`    | 2   | Control                   |
| `D2_Lapa`  | `notebooks/d2-lapa`  | 2   | Lapatinib (EGFR/HER2)     |
| `D4_DZ`    | `notebooks/d4-dz`    | 4   | Control                   |
| `D4_Lapa`  | `notebooks/d4-lapa`  | 4   | Lapatinib                 |
| `D4_AS`    | `notebooks/d4-as`    | 4   | AS1842856 (FOXO1)         |
| `G6`       | `notebooks/g6`       | 10  | Control for D10           |
| `D10_Lapa` | `notebooks/d10-lapa` | 10  | Lapatinib                 |

"DZ" in file names is the control condition. Each library pools organoids from
several donors, demultiplexed into `obs["participant"]`; cells called
`Doublet` or `Negative` there are dropped before pseudobulk DGE.

## Data flow

Every dataset runs the same stages. Paths are under `data/data-objects/`.

```
raw/egfDuod_<KEY>.h5ad                      input (see "Getting the data")
 │
 ├─ 0  Quality control      → qc/doublet_egfDuod_<KEY>.h5ad
 │     MT/ribo annotation, scanpy QC metrics, Scrublet doublet scores.
 │     Cells are not filtered here; scores are carried forward for inspection.
 │
 ├─ 1  Clustering           → clustered/clustered_egfDuod_<KEY>.h5ad
 │     HVGs → scale → regress S/G2M scores → PCA → kNN graph → Leiden → UMAP
 │
 ├─ 2  Doublet plotting       per-cluster doublet-score UMAPs (figures only;
 │                            D4_AS also writes the scores into the clustered object)
 ├─ 3  Marker plotting        dotplots / violins / UMAPs of src.markers genes
 │
 ├─ 4  CellAssign           → cellassign/<key>_predictions.h5ad
 │     scvi-tools CellAssign with a marker matrix from resources/cellassign/
 │
 ├─ 5  Manual labelling     → labelled/<key>_manual_labels.h5ad
 │     cluster → cell-type calls, refined with marker thresholds (obs["manual_label"])
 │
 └─ 6  Decoupler              CollecTRI TF activity (ULM) → figures/decoupler/
```

Then, across datasets:

```
labelled objects ─▶ dge/pseudobulk_deseq2.ipynb ─▶ data/analysis-files/dge/pydeseq-output/
                    (one run per comparison in dge/dge_comparisons.yaml)
DGE results      ─▶ gsea/gsea_prerank.ipynb     ─▶ data/analysis-files/gsea/gseapy-output/
```

Datasets differ in places: not every dataset has every step, and numbering does
not always match run order (`d2-lapa/5.1` runs after `5`, and D10 uses
`3-unsupervised_assignment` for CellAssign and `4-eec_subanalysis` for the EEC
subtypes). `python scripts/run_pipeline.py --list` prints the exact plan.

### D10 EEC sub-analysis

`d10-lapa/4-eec_subanalysis.ipynb` produces the EEC subtype figure in the
README. It subsets CellAssign-called EECs, reclusters them, gates peptide
hormones (GHRL, SST, CCK, GIP) by expression, holds out NEUROG3-high cells as
Early EECs, and resolves ARX+ cells without a clear hormone by kNN label
transfer. The UMAP layout is restored from a saved embedding
(`resources/embeddings/d10_eec_canonical_embedding.npz`) so the published
panel is reproduced exactly. Final labels are in `obs["granular_EEC_label_v5"]`.

## Running it

Run order and the environment for each notebook come from `pipeline.yaml`.

```bash
python scripts/run_pipeline.py --list              # the plan, and which inputs are present
python scripts/run_pipeline.py D2_DZ               # one dataset
python scripts/run_pipeline.py D2_DZ --from 3_Marker_plotting.ipynb   # resume
python scripts/run_pipeline.py --all               # everything, then DGE and GSEA
python scripts/run_pipeline.py --step dge --comparison d2_lapa_iscs
```

Executed copies (with outputs and any traceback) go to
`data/executed-notebooks/`; the tracked notebooks are not modified unless you
pass `--inplace`. Each notebook runs with the Jupyter kernel named after its
conda env. If you use a single environment for everything, add
`--kernel <name>`.

Notebooks can also be run by hand in Jupyter, in the order above. They find
`src/` from any working directory.

## Environments

| Env               | File                                  | Used for                         |
|-------------------|---------------------------------------|----------------------------------|
| `single-cell-env` | `conda/single-cell-env-minimal.yaml`  | steps 0–5, plotting, papermill   |
| `decoupler-env`   | `conda/decoupler-env.yaml`            | step 6 (Decoupler) and GSEA      |
| `pydeseq-env`     | `conda/pydeseq-env.yaml`              | pseudobulk DGE                   |

`conda/single-cell-env.yaml` is the full pinned export of the environment the
published results were produced with; use it for an exact rebuild on macOS ARM.

```bash
conda env create -f conda/single-cell-env-minimal.yaml
conda activate single-cell-env && python -m ipykernel install --user --name single-cell-env
# repeat for decoupler-env and pydeseq-env
```

## Checks

```bash
python scripts/check_notebooks.py    # static: imports, paths, step order. No data needed.
pytest tests/                        # src/ on synthetic data, plus the static checks
```

`check_notebooks.py` catches the ways notebooks break on a fresh clone: an
import of something no longer in `src/`, a file one step reads that no earlier
step writes, a path relative to the working directory, or a variable that only
existed in an old kernel session. Run it after renaming or moving anything.

## Conventions

### AnnData slots

| Slot                        | Contents                                                        |
|-----------------------------|-----------------------------------------------------------------|
| `layers["counts"]`          | raw integer counts. Used by CellAssign and pseudobulk DGE.      |
| `.X` in raw/ and qc/        | log-normalised (`log1p` of counts scaled to 10,000). Exception: D10's qc/ object holds raw counts in `.X` (log values in `layers["logcounts"]`), so clustering reads raw/ and takes only the QC obs columns. |
| `.X` in clustered/ onward   | **z-scored** (`sc.pp.scale`, clipped at 10). Use `layers["counts"]` or re-normalise for expression plots. |
| `obs["participant"]`        | donor (from freemuxlet); `Doublet` / `Negative` are not donors  |
| `obs["Condition"]`          | `Dz`, `Lapa` or `AS`                                            |
| `obs["leiden"]`             | clusters from step 1                                            |
| `obs["doublet_score"]`      | Scrublet score (step 0)                                         |
| `obs["initial_cellassign_prediction"]` | CellAssign argmax; `cell_assign_prob_<type>` holds probabilities |
| `obs["manual_label"]`       | final cell type, used by DGE                                    |

### Cell-type names

Some labelling notebooks spell cell types differently ("Secretory progenitor
cells" vs "Secretory PCs", "Inflammed" vs "Inflamed"). `src/palette.py` maps
every alias to one canonical name; colours (`get_color`) and pseudobulk
selection (`create_pseudobulk`) both go through it, so add new spellings there.

### Clustering parameters

From `src/preprocessing.run_clustering_pipeline` and the `1_Clustering` notebooks:

| Parameter                 | Value                    | Note                                       |
|---------------------------|--------------------------|--------------------------------------------|
| HVGs                      | 2,000 for D2_DZ, D2_Lapa, D4_DZ; 3,000 for D4_Lapa, D4_AS, G6, D10_Lapa | `flavor="seurat"` |
| Scaling                   | `max_value=10`           | applied to both the HVG subset and the full object |
| Cell cycle                | S and G2M scores regressed out | Tirosh genes in `src.markers.cell_cycle_markers` |
| kNN graph                 | `n_neighbors=100`, `use_rep="X"` | built on the scaled HVG matrix, not on PCs, so `n_pcs=30` has no effect. D4_AS also tries `n_neighbors=30`. |
| Leiden                    | resolution 3 tried first, 1.0 kept, `n_iterations=2` | D4_AS was tuned separately (1.3–2.0 tried); see its notebook |

### Notebook header

Every notebook starts with this, so it works from any working directory:

```python
import sys
from pathlib import Path
_p = Path(".").resolve()
while not (_p / "src" / "config.py").exists() and _p != _p.parent:
    _p = _p.parent
sys.path.insert(0, str(_p))
from src.config import DATA_DIR, FIGURES_DIR, DATASETS
```

Use the `src.config` paths for every read and write; avoid `os.chdir` and bare
relative paths.

## `src/` reference

| Module            | What it holds                                                          |
|-------------------|------------------------------------------------------------------------|
| `config.py`       | all paths, the `DATASETS` registry; creates output directories on import |
| `markers.py`      | marker gene dictionaries, `load_cellassign_matrix`, `filter_marker_dict` |
| `palette.py`      | cell-type and EEC-subtype colours, name normalisation                  |
| `preprocessing.py`| `run_qc_pipeline`, `run_clustering_pipeline`                           |
| `cell_assign.py`  | size factors, CellAssign wrapper, prediction annotation               |
| `dge.py`          | pseudobulk aggregation and PyDESeq2                                    |
| `gsea.py`         | `GSEAEnrichmentAnalyzer` (gseapy prerank)                              |
| `decoupler.py`    | `run_decoupler_safe` (CollecTRI + ULM on HVGs)                         |
| `plotting.py`     | print-ready UMAP panels (`umap_panel`, `gene_umaps`)                   |
| `figure_style.py` | slide-sized figure presets                                             |
| `vis.py`          | Veres-style panels, doublet-score UMAPs. Note: importing it sets small (5 pt) matplotlib fonts globally. |
| `loupe.py`        | h5ad → Loupe Browser `.cloupe` export (needs `loupepy`)                |

## Not in the pipeline

`pipeline.yaml` lists notebooks kept for the record but never run: per-dataset
sense checks (`1.1_*`), the superseded D10 EEC notebooks (`1.2`–`1.4`), the
pre-consolidation DGE notebooks in `dge/archive/`, and `exploratory/`.
