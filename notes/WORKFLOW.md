# EGFR Single-Cell RNA-seq Project — Workflow Guide

## 1. Project Overview

This project analyses **scRNA-seq data from human duodenal organoids** treated with EGFR pathway inhibitors, investigating how EGFR signalling affects intestinal stem cell differentiation and interferon responses.

### Datasets

| Key        | Timepoint | Condition            | Description                          |
|------------|-----------|----------------------|--------------------------------------|
| **D2_DZ**  | Day 2     | DMSO (control)       | Vehicle control                      |
| **D2_Lapa**| Day 2     | Lapatinib            | EGFR/HER2 dual inhibitor             |
| **D4_DZ**  | Day 4     | DMSO (control)       | Vehicle control                      |
| **D4_Lapa**| Day 4     | Lapatinib            | EGFR/HER2 dual inhibitor             |
| **D4_AS**  | Day 4     | Afatinib + SHP099    | Pan-ERBB + SHP2 inhibitor combo      |
| **D10_Lapa**| Day 10   | Lapatinib            | EGFR/HER2 dual inhibitor             |
| **G6**     | Day 10    | DMSO (control)       | Vehicle control for D10               |

### Cell Types Identified

- **ISCs** — Intestinal stem cells (LGR5+, AXIN2+, ASCL2+)
- **PCs** — Progenitor cells (DLL1+, NEUROG3+, CD44+)
- **Proliferating PCs** — Actively dividing progenitors (MKI67+)
- **Secretory PCs** — Secretory-lineage progenitors (INSM1+, NEUROD1+, SOX4+, ATOH1+)
- **Inflamed Secretory PCs** — REG3A/REG1B/REG1A-expressing secretory progenitors
- **Enterocytes** — Absorptive epithelial cells (KRT20+, FABP2+, ALPI+)
- **Goblet cells** — Mucus-secreting cells (MUC2+, FCGBP+, GFI1+)
- **EECs** — Enteroendocrine cells (CHGA+, CHGB+) with subtypes (X, D, I, K, EC)
- **NEUROG3+ progenitor cells** — Endocrine-committed progenitors

---

## 2. Data Flow

```
Seurat (R)
  │  r-src/  contains conversion scripts
  ▼
Raw H5AD files  (data/*.h5ad)
  │
  ├──  Step 0: Quality Control
  │    MT/ribo annotation, QC metrics, Scrublet doublet detection
  │    → analysis/qc/doublet_*.h5ad
  │
  ├──  Step 1: Clustering
  │    HVG selection → scaling → cell cycle regression →
  │    PCA → neighbors → Leiden → UMAP
  │    → analysis/clustered_*.h5ad
  │
  ├──  Step 2: Doublet Plotting (visual QC)
  │
  ├──  Step 3: Marker Plotting
  │    Dotplots, stacked violins, UMAP overlays for known markers
  │
  ├──  Step 4: CellAssign
  │    Probabilistic cell type assignment using marker gene matrix
  │    → analysis/cellassign_objects/*.h5ad
  │
  ├──  Step 5: Manual Labelling
  │    Scientific judgement: override/refine CellAssign predictions
  │    → analysis/manual_labelled*/*.h5ad
  │
  ├──  Step 6: Plotting & Visualisation
  │    Veres-style panels, population ratio bars
  │
  ├──  Step 7: Decoupler TF Activity
  │    CollecTRI resource, ULM method → TF activity scores
  │
  └──  DGE: Pseudobulk DESeq2
       Per cell-type pseudobulk → PyDESeq2 → volcano plots
       → data/dge/pydeseq-output/*.csv
```

---

## 3. Key Biological Decisions

### Marker Gene Choices
All marker definitions live in `src/markers.py`. The main epithelial markers (`cell_type_markers`) drive both the CellAssign classification and the visual validation plots.

### Manual Label Overrides
Step 5 is inherently manual — it's where you inspect CellAssign output and correct misclassifications based on biological knowledge. Common overrides:
- Cluster-level reassignment (e.g. "cluster 19 = EECs based on CHGA expression")
- Gene-expression thresholds (e.g. "MKI67+ cells in cluster 3 = Proliferating PCs")
- Inflamed secretory PCs identified by REG3A/REG1B/REG1A positivity

### Resolution Tuning
Leiden resolution varies per dataset (typically 1.0–2.0). Higher resolution captures more granular populations but risks overclustering. This is tuned iteratively using sense-check notebooks (Step 1.1).

---

## 4. Technical Decisions

| Parameter          | Typical Value    | Rationale                              |
|--------------------|------------------|----------------------------------------|
| `n_top_genes`      | 2000–3000        | HVGs for dimensionality reduction      |
| `n_neighbors`      | 30–100           | Graph construction (higher = smoother) |
| `n_pcs`            | 30               | PCA components for neighbors           |
| `leiden_resolution` | 1.0–2.0         | Clustering granularity (dataset-specific) |
| `max_value` (scale)| 10               | Clip scaled expression                 |
| Cell cycle         | Regressed out    | S_score + G2M_score via Tirosh genes   |
| Pseudobulk         | Sum per sample   | Standard for DESeq2 input              |
| DESeq2 design      | `~condition`     | Single-factor comparison               |

---

## 5. AnnData Structure & Coding Patterns

### AnnData Slots
```python
adata.X              # Log-normalised expression (after Seurat export)
adata.layers["counts"]  # Raw integer counts (used for CellAssign, DGE)
adata.obs             # Cell metadata (participant, Condition, leiden, manual_label, ...)
adata.var             # Gene metadata (mt, ribo, highly_variable, ...)
adata.obsm["X_pca"]   # PCA embedding
adata.obsm["X_umap"]  # UMAP embedding
adata.obsp            # distances, connectivities (neighbor graph)
adata.uns             # Colours, dendrogram, etc.
```

### Standard Notebook Header
```python
import sys
from pathlib import Path
sys.path.insert(0, str(Path("../..").resolve()))

from src.config import *
from src.markers import cell_type_markers
```

### pathlib for All Paths
```python
# Good
adata = sc.read_h5ad(ANALYSIS_DIR / "clustered_dataset.h5ad")

# Avoid
os.chdir("../../data")
adata = sc.read_h5ad("clustered_dataset.h5ad")
```

---

## 6. `src/` Module Reference

| Module             | Contents                                                    |
|--------------------|-------------------------------------------------------------|
| `src/config.py`    | `PROJECT_ROOT`, `DATA_DIR`, `ANALYSIS_DIR`, `DATASETS` dict |
| `src/markers.py`   | `cell_type_markers`, `EEC_markers`, `cell_cycle_markers`, `filter_present_genes()`, `filter_marker_dict()` |
| `src/palette.py`   | `celltype_palette`, `get_color()`, `normalize_celltype_name()` |
| `src/preprocessing.py` | `run_qc_pipeline()`, `run_clustering_pipeline()`        |
| `src/cell_assign.py` | `compute_size_factors()`, `run_cellassign()`, `annotate_predictions()`, `plot_cluster_assignment_heatmap()` |
| `src/dge.py`       | `create_pseudobulk()`, `merge_pseudobulk()`, `build_clinical_df()`, `run_pydeseq2()`, `plot_volcano()` |
| `src/vis.py`       | Veres-style panels (`plot_veres_panel()`, `plot_veres_panel_highlight()`), doublet UMAPs, marker UMAP panels |
| `src/gsea.py`      | GSEA/enrichment workflows (gseapy integration)              |
| `src/decoupler.py` | TF activity scoring (`run_decoupler_safe()` — CollecTRI + ULM, 5K HVGs) |
| `src/loupe.py`     | AnnData to Loupe Browser `.cloupe` conversion               |
| `src/utils.py`     | **Shim** — re-exports from `markers` and `palette` for backwards compatibility |

---

## 7. DGE Workflow

All DGE comparisons are defined in `notebooks/dge/dge_comparisons.yaml` and run through a single parameterised notebook (`notebooks/dge/pseudobulk_deseq2.ipynb`).

To run a comparison:
1. Open `pseudobulk_deseq2.ipynb`
2. Set `COMPARISON = "d2_lapa_iscs"` (or any key from the YAML)
3. Run all cells

Or via papermill:
```bash
papermill pseudobulk_deseq2.ipynb output.ipynb -p COMPARISON d2_lapa_iscs
```

Archived per-comparison notebooks are in `notebooks/dge/archive/` for reference.

---

## 8. Environment Setup

### Three Conda Environments

| Environment         | Used For                                | Spec File                          |
|---------------------|-----------------------------------------|------------------------------------|
| `single-cell-env`   | Steps 0–5 (QC → Labelling → Plotting)  | `conda/single-cell-env-minimal.yaml` |
| `pydeseq-env`       | DGE notebooks (pseudobulk DESeq2)       | `conda/pydeseq-env.yaml`          |
| `decoupler-env`     | TF activity & enrichment analysis       | `conda/decoupler-env.yaml`        |

CellAssign (Step 4) runs under `single-cell-env` since it includes `scvi-tools`.

### Creating from Spec
```bash
conda env create -f conda/single-cell-env-minimal.yaml
conda env create -f conda/pydeseq-env.yaml
conda env create -f conda/decoupler-env.yaml
```

### Pre-packed Environments
Packed tarballs are stored in `conda/` and as release assets on the project remote.
To install from a tarball:
```bash
mkdir -p ~/envs/single-cell-env
tar -xzf conda/egfr-single-cell-0.1.0.tar.gz -C ~/envs/single-cell-env
conda activate ~/envs/single-cell-env
conda-unpack
```

### Packing an Environment
```bash
conda pack -n single-cell-env -o conda/egfr-single-cell-<version>.tar.gz
conda pack -n pydeseq-env -o conda/pydeseq-env-osx-arm64.tar.gz
conda pack -n decoupler-env -o conda/decoupler-env-osx-arm64.tar.gz
```

---

## 9. Directory Structure

```
single-cell/
├── src/                         # Reusable Python modules
│   ├── config.py                # Paths & dataset registry
│   ├── markers.py               # Marker gene definitions
│   ├── palette.py               # Colour palette
│   ├── preprocessing.py         # QC & clustering pipelines
│   ├── cell_assign.py           # CellAssign wrapper
│   ├── dge.py                   # Pseudobulk DGE pipeline
│   ├── vis.py                   # Visualisation (Veres panels, doublets, markers)
│   ├── gsea.py                  # GSEA/enrichment workflows
│   ├── decoupler.py             # TF activity scoring (CollecTRI + ULM)
│   ├── loupe.py                 # AnnData → Loupe Browser conversion
│   └── utils.py                 # Backwards-compat shim
├── notebooks/
│   ├── d2-dz/                   # Dataset-specific analysis notebooks
│   ├── d2-lapa/
│   ├── d4-as/
│   ├── d4-dz/
│   ├── d4-lapa/
│   ├── d10-lapa/                # Most active — includes runner scripts:
│   │   ├── _run_8.py            #   Donor panel figure generation
│   │   └── _run_decoupler_panels.py  # TF activity 2x2 panels
│   ├── g6/
│   ├── cellrank/                # CellRank trajectory inference
│   ├── dge/                     # DGE: config YAML + parameterised notebook
│   ├── re-labelling/            # Manual relabelling pipeline
│   └── exploratory/             # Exploratory / one-off analysis
├── figures/
│   ├── d10-for-publication/     # Publication figures
│   │   ├── decoupler/           #   TF activity panels (all datasets)
│   │   └── eec-gene-panels/     #   EEC marker gene UMAPs
│   ├── d10-lapa-donors/         # Donor-level composition plots
│   ├── decoupler/               # Per-dataset TF activity UMAPs
│   ├── gsea-publications/       # GSEA enrichment panels
│   └── ...                      # marker-genes, qc, clusters, etc.
├── analysis/                    # Processed data (gitignored)
│   ├── cellassign_objects/      #   CellAssign predictions
│   ├── manual_labelled/         #   Manual labels v1 (D4, G6)
│   ├── manual_labelled_2/       #   Manual labels v2 (D2) + EEC subsets
│   ├── loupe/                   #   .cloupe files
│   └── qc/                      #   QC outputs
├── data/                        # DGE outputs (gitignored; raw h5ad on external volume)
├── conda/                       # Environment specs + packed tarballs
├── utilities/                   # CellAssign marker CSVs
├── r-src/                       # R conversion scripts
├── notes/                       # README, this workflow guide
└── archive/                     # Old scripts, summaries, duplicate figures
```
