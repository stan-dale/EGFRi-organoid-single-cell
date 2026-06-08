"""
Same per-gene UMAP grid as step 5, but on the D4 CONTROL (D4_DZ) data, for a
side-by-side read of whether the inflammation programme is lapatinib-driven.

NOTE: this is the FULL control object on its own global UMAP (no secretory
subsetting yet), so it is not strictly apples-to-apples with the lapatinib
secretory-subset embedding in 05 — it's a first look at expression/localisation.
First panel = manual_label for orientation (incl. the 228 REG-gated inflamed cells).
"""
import scanpy as sc
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
FIG = HERE / "figures"

F = "/Users/standale/Library/CloudStorage/Dropbox-BCH/Stanley Dale/David and Stanley/EGFRi/data/data-objects/3.final-h5ad/d4_dz_manual_labels.h5ad"
ad = sc.read_h5ad(F)
print("loaded control", ad.shape)

# log1p CP10k from counts -> lognorm layer (X in these objects is scaled)
ad.layers["lognorm"] = ad.layers["counts"].copy()
tmp = sc.AnnData(ad.layers["counts"].copy(), obs=ad.obs[[]].copy(), var=ad.var[[]].copy())
sc.pp.normalize_total(tmp, target_sum=1e4); sc.pp.log1p(tmp)
ad.layers["lognorm"] = tmp.X; del tmp

GENES = [
    "CCL20", "CXCL8", "CXCL1", "CXCL2", "CXCL3", "NFKBIA", "IER3", "GADD45B", "JUN", "PPP1R15A",
    "LCN2", "NOS2", "DUOX2", "DUOXA2", "REG3A",
    "CD74", "HLA-DRA", "CIITA", "STAT1", "IFIT3",
]
present = [g for g in GENES if g in ad.var_names]
missing = [g for g in GENES if g not in ad.var_names]
if missing:
    print("not in control var_names (skipped):", missing)

panels = ["manual_label"] + present
ncol = 6
nrow = int(np.ceil(len(panels) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(4.0 * ncol, 3.6 * nrow))
for ax, key in zip(axes.flat, panels):
    if key == "manual_label":
        sc.pl.umap(ad, color=key, ax=ax, show=False, size=6, title="manual_label",
                   legend_loc="right margin", legend_fontsize=6)
    else:
        sc.pl.umap(ad, color=key, ax=ax, show=False, cmap="magma", size=6,
                   layer="lognorm", title=key, colorbar_loc="right")
for ax in axes.flat[len(panels):]:
    ax.axis("off")
fig.suptitle("D4 CONTROL (D4_DZ) — same inflammation markers (lognorm), full UMAP", y=1.004, fontsize=14)
fig.tight_layout()
fig.savefig(FIG / "06_control_inflammation_gene_umaps.png", dpi=140, bbox_inches="tight")
plt.close(fig)
print("wrote", FIG / "06_control_inflammation_gene_umaps.png", "with", len(present), "genes")
