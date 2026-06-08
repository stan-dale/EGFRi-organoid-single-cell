"""
Per-gene UMAP grid (same style as 01_key_markers / the Alarmins CCL20/CXCL8 panels):
one panel per gene, magma, lognorm expression, on the SECRETORY-SUBSET embedding
(X_umap_sub) so it lines up with the subclusters.

Curated inflammation panel centred on the chemokine / NF-kB programme that CCL20 &
CXCL8 belong to, plus the ISP-core and a few MHC/IFN genes for context.
Swap GENES for any list you want.
"""
import scanpy as sc
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
FIG = HERE / "figures"; OUT = HERE / "outputs"

ad = sc.read_h5ad(OUT / "d4_lapa_secretory.h5ad")

GENES = [
    # chemokine / NF-kB / immediate-early (the CCL20/CXCL8 programme, ~ subcluster c11)
    "CCL20", "CXCL8", "CXCL1", "CXCL2", "CXCL3", "NFKBIA", "IER3", "GADD45B", "JUN", "PPP1R15A",
    # ISP-core / LND inflammation (~ c2)
    "LCN2", "NOS2", "DUOX2", "DUOXA2", "REG3A",
    # MHC-II / interferon (~ c3)
    "CD74", "HLA-DRA", "CIITA", "STAT1", "IFIT3",
]
present = [g for g in GENES if g in ad.var_names]
missing = [g for g in GENES if g not in ad.var_names]
if missing:
    print("not in var_names (skipped):", missing)

ncol = 5
nrow = int(np.ceil(len(present) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(4.0 * ncol, 3.6 * nrow))
for ax, g in zip(axes.flat, present):
    sc.pl.embedding(ad, "X_umap_sub", color=g, ax=ax, show=False, cmap="magma",
                    size=10, layer="lognorm", title=g, colorbar_loc="right")
for ax in axes.flat[len(present):]:
    ax.axis("off")
fig.suptitle("D4_Lapa secretory subset — inflammation markers (lognorm)", y=1.005, fontsize=14)
fig.tight_layout()
fig.savefig(FIG / "05_inflammation_gene_umaps.png", dpi=140, bbox_inches="tight")
plt.close(fig)
print("wrote", FIG / "05_inflammation_gene_umaps.png", "with", len(present), "genes")
