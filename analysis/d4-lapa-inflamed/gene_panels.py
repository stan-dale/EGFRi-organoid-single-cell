"""
Generate per-gene UMAP panels on the D4_Lapa embedding, pre-shaped for the Figma
slide slots via figure_style. Each gene gets its own panel (magma, lognorm), laid
out in a single row sized to the chosen slot.

    python gene_panels.py            # builds the configured GROUPS below

Edit GROUPS to add more figures: (filename_stem, [genes...], slot).
"""
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc

import figure_style as fs

HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
DATA = HERE / "outputs" / "d4_lapa_reclustered.h5ad"

# (output stem, genes, slot, ncols) — one PNG per entry, one panel per gene.
# ncols=None -> single row (one column per gene).
GROUPS = [
    ("panel_LGR5_OLFM4", ["LGR5", "OLFM4"], "wide", None),
    ("panel_chemokine", ["CXCL1", "CXCL8", "CCL20", "LCN2"], "single", 2),
    ("panel_MHCII", ["CD74", "HLA-DRA"], "wide", None),
    ("panel_AMP", ["REG1A", "REG3A"], "wide", None),
]


def make_gene_panel(ad, genes, slot, out_path, ncols=None):
    present = [g for g in genes if g in ad.var_names]
    missing = [g for g in genes if g not in ad.var_names]
    if missing:
        print("  skipped (not in var_names):", missing)

    ncols = ncols or len(present)
    nrows = int(np.ceil(len(present) / ncols))

    fs.setup_style()
    fig, axes = fs.make_figure(slot=slot, ncols=ncols, nrows=nrows)
    axes = np.atleast_1d(axes).ravel()

    for ax, g in zip(axes, present):
        sc.pl.embedding(
            ad, "X_umap", color=g, ax=ax, show=False,
            cmap="magma", size=10, layer="lognorm", title=g,
            frameon=False, colorbar_loc="right",
        )
        ax.set_facecolor("none")
        ax.set_aspect("equal")

    for ax in axes[len(present):]:  # blank any unused grid cells
        ax.axis("off")

    p = fs.save_for_slide(fig, str(out_path))
    plt.close(fig)
    return p


if __name__ == "__main__":
    ad = sc.read_h5ad(DATA)
    print("loaded", ad.shape, "from", DATA.name)
    from PIL import Image
    for stem, genes, slot, ncols in GROUPS:
        out = FIG / f"{stem}_{slot}.png"
        p = make_gene_panel(ad, genes, slot, out, ncols=ncols)
        with Image.open(p) as im:
            w, h = im.size
        print(f"wrote {p}  ({w}x{h} px, slot={slot})")
