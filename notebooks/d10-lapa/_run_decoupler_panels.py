#!/usr/bin/env python
"""
Generate TF activity 2x2 panels for a single dataset.

Creates a 2x2 figure: cell type UMAP (top-left) + HNF4A, IRF1, STAT1
activity scores (RdBu_r, centred at 0). Uses 5K HVGs and CollecTRI ULM.

Usage (run each dataset as a separate process to manage memory):
    python _run_decoupler_panels.py D10_Lapa initial_cellassign_prediction
    python _run_decoupler_panels.py D2_DZ manual_label
    python _run_decoupler_panels.py D2_Lapa manual_label
    python _run_decoupler_panels.py D4_DZ manual_label
    python _run_decoupler_panels.py D4_Lapa manual_label

For D10_Lapa, use_hvg=False is recommended (28K genes, 33K cells fits in
memory) to maximise TF coverage. For larger datasets, 5K HVGs is the
sweet spot for IRF1 recovery (>=25 targets) without OOM.
"""
import sys, os, warnings, gc
warnings.filterwarnings("ignore")
os.environ["MPLBACKEND"] = "Agg"

from pathlib import Path
_p = Path(__file__).resolve().parent
while not (_p / "src" / "config.py").exists() and _p != _p.parent:
    _p = _p.parent
sys.path.insert(0, str(_p))

from src.config import DATASETS, FIGURES_DIR
from src.palette import celltype_palette
import scanpy as sc
import decoupler as dc
import matplotlib.pyplot as plt

ds_key = sys.argv[1]       # e.g. "D2_DZ"
label_key = sys.argv[2]    # e.g. "manual_label"

SAVE_DIR = FIGURES_DIR / "d10-for-publication" / "decoupler"
SAVE_DIR.mkdir(parents=True, exist_ok=True)
tfs = ["HNF4A", "IRF1", "STAT1"]

print(f"=== {ds_key} ===", flush=True)
adata = sc.read_h5ad(DATASETS[ds_key]["labelled"])
print(f"Loaded {adata.n_obs} x {adata.n_vars}", flush=True)

# Force fresh 5K HVGs
if "highly_variable" in adata.var.columns:
    adata.var.drop(columns=["highly_variable"], inplace=True)

if "counts" in adata.layers:
    adata_tmp = adata.copy()
    adata_tmp.X = adata_tmp.layers["counts"]
    sc.pp.highly_variable_genes(adata_tmp, n_top_genes=5000, flavor="seurat_v3")
    adata.var["highly_variable"] = adata_tmp.var["highly_variable"]
    del adata_tmp
    gc.collect()
else:
    sc.pp.highly_variable_genes(adata, n_top_genes=5000)

n_hvg = adata.var["highly_variable"].sum()
print(f"Subsetting to {n_hvg} HVGs", flush=True)
adata = adata[:, adata.var["highly_variable"]].copy()
gc.collect()

# Palette
adata.obs[label_key] = adata.obs[label_key].astype("category")
cats = adata.obs[label_key].cat.categories
adata.uns[f"{label_key}_colors"] = [celltype_palette.get(c, "#808080") for c in cats]

# ULM
collectri = dc.op.collectri(organism="human")
collectri_sub = collectri[collectri["target"].isin(adata.var_names)].copy()
print(f"{collectri_sub['source'].nunique()} TFs, {len(collectri_sub)} interactions", flush=True)

dc.mt.ulm(data=adata, net=collectri_sub, tmin=25, verbose=True)
score = dc.pp.get_obsm(adata, key="score_ulm")
print(f"{score.n_vars} TFs computed", flush=True)

for tf in tfs:
    print(f"  {tf}: {'YES' if tf in score.var_names else 'NO'}", flush=True)

# Plot
fig, axes = plt.subplots(2, 2, figsize=(10, 9))
sc.pl.umap(score, color=label_key, ax=axes[0, 0], show=False, frameon=False, title="")

for ax, tf in zip([axes[0, 1], axes[1, 0], axes[1, 1]], tfs):
    if tf in score.var_names:
        sc.pl.umap(score, color=tf, ax=ax, show=False, frameon=False,
                   cmap="RdBu_r", vcenter=0, title=f"{tf} score")
        ax.set_title(f"{tf} score", fontsize=10)
    else:
        ax.set_title(f"{tf} (n/a)", fontsize=10)
        ax.axis("off")

plt.tight_layout()
fname = f"{ds_key.lower()}_tf_activity_2x2.pdf"
fig.savefig(SAVE_DIR / fname, dpi=600, bbox_inches="tight")
plt.close(fig)
print(f"Saved {fname}", flush=True)
