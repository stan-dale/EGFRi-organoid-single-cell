"""
OLFM4 on D4 control (D4_DZ) vs D4 lapatinib (D4_Lapa), side by side on a single
SHARED log-normalised colour scale so expression level is directly comparable.

Each object is on its own UMAP embedding (layout is not comparable), but the
colour mapping is identical (same vmin/vmax, one shared colorbar), so the
intensity read across the two is apples-to-apples. lognorm = log1p CP10k from
layers['counts'] (X is scaled/clipped in these objects).

    python compare_olfm4.py
"""
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import scanpy as sc

import figure_style as fs

HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"

GENE = "OLFM4"
CONTROL = "/Users/standale/Library/CloudStorage/Dropbox-BCH/Stanley Dale/David and Stanley/EGFRi/data/data-objects/3.final-h5ad/d4_dz_manual_labels.h5ad"
LAPA = HERE / "outputs" / "d4_lapa_reclustered.h5ad"


def lognorm_vec(ad, gene):
    """log1p CP10k expression for one gene, from layers['counts']."""
    if "lognorm" in ad.layers:
        x = ad[:, gene].layers["lognorm"]
    else:
        tmp = sc.AnnData(ad.layers["counts"].copy(),
                         obs=ad.obs[[]].copy(), var=ad.var[[]].copy())
        sc.pp.normalize_total(tmp, target_sum=1e4)
        sc.pp.log1p(tmp)
        x = tmp[:, gene].X
    return np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()


ctrl = sc.read_h5ad(CONTROL)
lapa = sc.read_h5ad(LAPA)
print("control", ctrl.shape, "| lapatinib", lapa.shape)

c_expr = lognorm_vec(ctrl, GENE)
l_expr = lognorm_vec(lapa, GENE)
c_xy = ctrl.obsm["X_umap"]
l_xy = lapa.obsm["X_umap"]

# Shared colour scale: vmin=0, vmax = 99th pct of the pooled expression.
vmax = float(np.percentile(np.concatenate([c_expr, l_expr]), 99))
norm = Normalize(vmin=0, vmax=vmax)
cmap = "magma"
print(f"shared scale: vmin=0, vmax={vmax:.2f} (p99 pooled)")

fs.setup_style()
fig, axes = fs.make_figure(slot="wide", ncols=2, nrows=1)
panels = [("D4 control (D4_DZ)", c_xy, c_expr), ("D4 lapatinib (D4_Lapa)", l_xy, l_expr)]
for ax, (title, xy, expr) in zip(axes, panels):
    order = np.argsort(expr)  # draw high-expressing cells on top
    ax.scatter(xy[order, 0], xy[order, 1], c=expr[order], cmap=cmap, norm=norm,
               s=6, linewidths=0, rasterized=True)
    ax.set_title(title)
    ax.set_facecolor("none")
    ax.set_aspect("equal")
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])

sm = ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
cbar = fig.colorbar(sm, ax=list(axes), fraction=0.025, pad=0.02)
cbar.set_label(f"{GENE} (log1p CP10k)")
fig.suptitle(f"{GENE} — D4 control vs lapatinib (shared scale)", y=0.98)

out = FIG / "compare_OLFM4_control_vs_lapa_wide.png"
p = fs.save_for_slide(fig, str(out))
plt.close(fig)
from PIL import Image
with Image.open(p) as im:
    w, h = im.size
print(f"wrote {p}  ({w}x{h} px)")
