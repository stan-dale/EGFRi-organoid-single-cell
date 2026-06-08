"""
Higher-contrast version of the OLFM4-style control-vs-lapatinib comparison, for
any gene. Same SHARED log-normalised scale (so the two panels stay comparable),
but the colour scale is clipped to a tunable lower percentile (PCT) instead of
p99 — most cells then span the full colormap and the control->lapatinib shift is
much more visible. Each panel is annotated with % expressing and median(lognorm
among expressing) so the difference is quantitative, not just visual.

Does NOT overwrite compare_olfm4.py's p99 output (different filename suffix).

    python compare_genes_contrast.py
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

GENES = ["OLFM4", "LGR5"]
PCT = 90          # colour-scale clip percentile (lower => higher contrast). p99 was the original.
CMAP = "magma"

CONTROL = "/Users/standale/Library/CloudStorage/Dropbox-BCH/Stanley Dale/David and Stanley/EGFRi/data/data-objects/3.final-h5ad/d4_dz_manual_labels.h5ad"
LAPA = HERE / "outputs" / "d4_lapa_reclustered.h5ad"


def lognorm_vec(ad, gene):
    if "lognorm" in ad.layers:
        x = ad[:, gene].layers["lognorm"]
    else:
        tmp = sc.AnnData(ad.layers["counts"].copy(),
                         obs=ad.obs[[]].copy(), var=ad.var[[]].copy())
        sc.pp.normalize_total(tmp, target_sum=1e4)
        sc.pp.log1p(tmp)
        x = tmp[:, gene].X
    return np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()


def stats(expr):
    pos = expr > 0
    pct = 100 * pos.mean()
    med = float(np.median(expr[pos])) if pos.any() else 0.0
    return pct, med


print("loading objects...")
ctrl = sc.read_h5ad(CONTROL)
lapa = sc.read_h5ad(LAPA)
c_xy = ctrl.obsm["X_umap"]
l_xy = lapa.obsm["X_umap"]
from PIL import Image

for gene in GENES:
    c_expr = lognorm_vec(ctrl, gene)
    l_expr = lognorm_vec(lapa, gene)

    # vmax from EXPRESSING cells only — robust for sparse genes (LGR5) where the
    # PCT percentile over all cells would be 0 and collapse the scale.
    pooled = np.concatenate([c_expr, l_expr])
    nonzero = pooled[pooled > 0]
    vmax = float(np.percentile(nonzero, PCT)) if nonzero.size else 1.0
    if vmax <= 0:
        vmax = float(pooled.max()) or 1.0
    norm = Normalize(vmin=0, vmax=vmax)
    c_pct, c_med = stats(c_expr)
    l_pct, l_med = stats(l_expr)
    print(f"{gene}: shared vmax(nonzero p{PCT})={vmax:.2f} | "
          f"control {c_pct:.1f}% expr med {c_med:.2f} | lapa {l_pct:.1f}% expr med {l_med:.2f}")

    fs.setup_style()
    fig, axes = fs.make_figure(slot="wide", ncols=2, nrows=1)
    panels = [
        (f"D4 control (D4_DZ)\n{c_pct:.1f}% expressing · median {c_med:.2f}", c_xy, c_expr),
        (f"D4 lapatinib (D4_Lapa)\n{l_pct:.1f}% expressing · median {l_med:.2f}", l_xy, l_expr),
    ]
    for ax, (title, xy, expr) in zip(axes, panels):
        order = np.argsort(expr)
        ax.scatter(xy[order, 0], xy[order, 1], c=expr[order], cmap=CMAP, norm=norm,
                   s=6, linewidths=0, rasterized=True)
        ax.set_title(title, fontsize=13)
        ax.set_facecolor("none"); ax.set_aspect("equal")
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xticks([]); ax.set_yticks([])

    sm = ScalarMappable(norm=norm, cmap=CMAP); sm.set_array([])
    cbar = fig.colorbar(sm, ax=list(axes), fraction=0.025, pad=0.02, extend="max")
    cbar.set_label(f"{gene} (log1p CP10k)")
    fig.suptitle(f"{gene} — D4 control vs lapatinib (shared scale, clipped at p{PCT})", y=0.98)

    out = FIG / f"compare_{gene}_control_vs_lapa_contrast_wide.png"
    p = fs.save_for_slide(fig, str(out))
    plt.close(fig)
    with Image.open(p) as im:
        w, h = im.size
    print(f"  wrote {p}  ({w}x{h} px)")
