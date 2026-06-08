"""
Step 1 — D4_Lapa exploration of inflamed secretory progenitors.

Design decisions (per user):
  - KEEP the existing PCA/UMAP embedding (reuse X_pca for the graph, X_umap for plots).
  - Keep all cells, but FLAG freemuxlet Doublet/Negative so we can sanity-check artifacts.
  - X is scaled/clipped (-7..10); recompute log1p-CP10k from layers['counts'] for scoring/plots.

Outputs:
  - outputs/d4_lapa_reclustered.h5ad  (adds leiden_r* + lognorm layer + lineage scores)
  - figures/01_*.png
"""
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

sc.settings.verbosity = 1
HERE = Path(__file__).resolve().parent.parent
FIG = HERE / "figures"; OUT = HERE / "outputs"
FIG.mkdir(exist_ok=True); OUT.mkdir(exist_ok=True)

F = "/Users/standale/Library/CloudStorage/Dropbox-BCH/Stanley Dale/David and Stanley/EGFRi/data/data-objects/3.final-h5ad/d4_lapa_manual_labels.h5ad"
ad = sc.read_h5ad(F)
print("loaded", ad.shape)

# ---- recompute log1p CP10k from counts into a layer (X stays scaled) ----
ad.layers["lognorm"] = ad.layers["counts"].copy()
tmp = sc.AnnData(ad.layers["counts"].copy(), obs=ad.obs[[]].copy(), var=ad.var[[]].copy())
sc.pp.normalize_total(tmp, target_sum=1e4)
sc.pp.log1p(tmp)
ad.layers["lognorm"] = tmp.X
del tmp

# ---- flag demultiplex artifacts ----
ad.obs["is_artifact"] = ad.obs["participant"].isin(["Doublet", "Negative"]).map({True: "artifact", False: "singlet"}).astype("category")

# ---- reuse existing embedding: rebuild neighbors on X_pca, recluster ----
n_pcs = min(30, ad.obsm["X_pca"].shape[1])
sc.pp.neighbors(ad, n_neighbors=30, n_pcs=n_pcs, use_rep="X_pca")
resolutions = [0.2, 0.5, 1.0, 1.5, 2.0]
for r in resolutions:
    key = f"leiden_r{r}"
    sc.tl.leiden(ad, resolution=r, key_added=key, flavor="igraph", n_iterations=2, directed=False)
    print(key, "->", ad.obs[key].nunique(), "clusters")

# ---- lineage signatures (defining the secretory-biased PC compartment) ----
lineage = {
    "ISC":           ["LGR5", "OLFM4", "ASCL2", "AXIN2", "SMOC2", "RNF43"],
    "Proliferation": ["MKI67", "TOP2A", "CDK1", "PCNA", "UBE2C", "BIRC5"],
    "Secretory_prog":["ATOH1", "NEUROG3", "DLL1", "SPDEF", "INSM1", "NEUROD1", "SOX4", "HES6"],
    "Absorptive":    ["FABP2", "FABP1", "APOA1", "APOA4", "ALPI", "KRT20", "RBP2"],
    "Goblet":        ["MUC2", "FCGBP", "TFF3", "SPINK4", "ZG16", "CLCA1"],
    "EEC":           ["CHGA", "CHGB", "PAX4", "PCSK1", "SCG5"],
    "Paneth_like":   ["DEFA5", "DEFA6", "LYZ", "REG3A"],
}
def present(genes): return [g for g in genes if g in ad.var_names]
for name, genes in lineage.items():
    g = present(genes)
    sc.tl.score_genes(ad, g, score_name=f"sig_{name}", use_raw=False, layer="lognorm")
    print(f"sig_{name}: {len(g)}/{len(genes)} genes -> {g}")

ad.write(OUT / "d4_lapa_reclustered.h5ad")
print("wrote", OUT / "d4_lapa_reclustered.h5ad")

# ---- figures ----
# A) UMAP: resolutions + existing manual_label + donor + artifact
panels = [f"leiden_r{r}" for r in resolutions] + ["manual_label", "participant", "is_artifact"]
fig, axes = plt.subplots(2, 4, figsize=(22, 10))
for ax, key in zip(axes.flat, panels):
    sc.pl.umap(ad, color=key, ax=ax, show=False, legend_loc="on data" if key.startswith("leiden") else "right margin",
               legend_fontsize=6, size=6, title=key)
for ax in axes.flat[len(panels):]: ax.axis("off")
fig.tight_layout(); fig.savefig(FIG / "01_clustering_overview.png", dpi=130); plt.close(fig)

# B) UMAP: lineage signature scores
sigs = [f"sig_{k}" for k in lineage]
fig, axes = plt.subplots(2, 4, figsize=(22, 10))
for ax, key in zip(axes.flat, sigs):
    sc.pl.umap(ad, color=key, ax=ax, show=False, cmap="viridis", size=6, title=key)
for ax in axes.flat[len(sigs):]: ax.axis("off")
fig.tight_layout(); fig.savefig(FIG / "01_lineage_scores.png", dpi=130); plt.close(fig)

# C) key individual markers (lognorm)
markers = ["LGR5","OLFM4","MKI67","ATOH1","NEUROG3","DLL1","MUC2","FABP2","CHGA",
           "REG3A","REG1A","LCN2","DUOX2","CXCL1","S100A9","CD74","HLA-DRA","TFF3"]
markers = [m for m in markers if m in ad.var_names]
ncol=6; nrow=int(np.ceil(len(markers)/ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(4*ncol, 3.6*nrow))
for ax, g in zip(axes.flat, markers):
    sc.pl.umap(ad, color=g, ax=ax, show=False, cmap="magma", size=6, layer="lognorm", title=g)
for ax in axes.flat[len(markers):]: ax.axis("off")
fig.tight_layout(); fig.savefig(FIG / "01_key_markers.png", dpi=130); plt.close(fig)

# D) lineage score means per cluster (r1.0) table
r="leiden_r1.0"
tab = ad.obs.groupby(r)[[f"sig_{k}" for k in lineage]].mean()
tab.to_csv(OUT / "01_lineage_means_r1.0.csv")
print(tab.round(3).to_string())
print("DONE")
