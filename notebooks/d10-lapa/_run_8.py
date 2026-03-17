#!/usr/bin/env python
"""Run 8_Donor_panels notebook logic — memory-optimised version."""
import matplotlib
matplotlib.use("Agg")

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gc
import sys
from pathlib import Path

_p = Path(__file__).resolve().parent
while not (_p / "src" / "config.py").exists() and _p != _p.parent:
    _p = _p.parent
sys.path.insert(0, str(_p))

from src.config import CELLASSIGN_DIR, LABELLED_DIR, FIGURES_DIR
from src.palette import celltype_palette
from src.vis import setup_matplotlib_params, plot_veres_panel, plot_veres_panel_highlight

setup_matplotlib_params()

# ---------- Config ----------
LABEL_KEY = "initial_cellassign_prediction"
DONOR_KEY = "participant"
REAL_DONORS = ["H329", "H896", "H897"]
EEC_LABEL_KEY = "granular_EEC_label_v3"

SAVE_DIR = FIGURES_DIR / "d10-lapa-donors"
SAVE_DIR.mkdir(parents=True, exist_ok=True)

LABEL_ORDER = [
    "EECs", "Secretory PCs", "Goblet cells",
    "Proliferating PCs", "ISCs", "PCs", "Enterocytes",
]
EEC_LABEL_ORDER = [
    "X cells", "D cells", "I cells",
    "K cells", "Enterochromaffin cells", "Early EECs",
]
# EEC subtype palette — matches Veres cycle from 1.4_EEC_subanalysis
# (Set1_9 in label_order, with Early EECs overridden to grey)
eec_palette = {
    "X cells":                "#e41a1c",  # red
    "D cells":                "#377eb8",  # blue
    "I cells":                "#4daf4a",  # green
    "K cells":                "#984ea3",  # purple
    "Enterochromaffin cells": "#ff7f00",  # orange
    "Early EECs":             "#999999",  # grey
}

from palettable.colorbrewer.qualitative import Dark2_3
from matplotlib.colors import to_hex
donor_palette = {
    d: to_hex([c / 255.0 for c in rgb])
    for d, rgb in zip(REAL_DONORS, Dark2_3.colors)
}

DPI = 600


def load_slim(path, donor_key, real_donors):
    """Load h5ad keeping only obs + obsm (no X matrix) to save memory."""
    import anndata as ad
    import h5py

    adata = sc.read_h5ad(path, backed="r")
    mask = adata.obs[donor_key].isin(real_donors).values
    obs = adata.obs[mask].copy()
    umap = adata.obsm["X_umap"][mask].copy()
    adata.file.close()

    # Drop unused categorical levels (e.g. Doublet/Negative)
    for col in obs.columns:
        if hasattr(obs[col], "cat"):
            obs[col] = obs[col].cat.remove_unused_categories()

    # Build a minimal AnnData with no X
    n = mask.sum()
    slim = ad.AnnData(
        X=np.zeros((n, 0), dtype=np.float32),
        obs=obs,
    )
    slim.obsm["X_umap"] = umap
    return slim


def plot_veres_triptych(
    adata, split_key, donors, label_key, label_order,
    palette, title=None, save=None, dpi=600,
):
    """Three Veres-style UMAP scatters side-by-side with one shared legend."""
    from src.vis import make_label_params, prepare_for_scatter
    import matplotlib.gridspec as gridspec

    cats = label_order
    lp = make_label_params(cats, palette=palette)

    mm = 25.4
    panel_w = 89 / 2 / mm  # single Veres panel width
    fig_w = panel_w * 3 + 0.6  # 3 panels + legend room
    fig_h = panel_w * 0.75
    fig = plt.figure(figsize=(fig_w, fig_h), dpi=dpi)

    # gridspec: 3 scatter columns + 1 narrow legend column
    gs = gridspec.GridSpec(1, 4, figure=fig, width_ratios=[1, 1, 1, 0.4],
                           wspace=0.05)

    s_black, s_white, s_type = 4, 2, 1.5

    for i, donor in enumerate(donors):
        ad = adata[adata.obs[split_key] == donor].copy()
        X = ad.obsm["X_umap"]
        labels = ad.obs[label_key].astype(str).values
        proj, rgb = prepare_for_scatter(X, labels, lp)

        ax = fig.add_subplot(gs[0, i])
        ax.scatter(proj[:, 0], proj[:, 1], c="k", s=s_black, edgecolor="none", rasterized=True)
        ax.scatter(proj[:, 0], proj[:, 1], c="w", s=s_white, edgecolor="none", rasterized=True)
        ax.scatter(proj[:, 0], proj[:, 1], c=rgb, s=s_type, alpha=0.7, edgecolor="none", rasterized=True)
        ax.set_xticks([]); ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xlabel(f"{donor} ({ad.n_obs:,})", fontsize=6, fontweight="bold", labelpad=4)
        ax.set_ylabel("")
        ax.set_title("")

    # Shared legend
    ax_leg = fig.add_subplot(gs[0, 3])
    ax_leg.axis("off")
    present = [c for c in cats if c in set(adata.obs[label_key].astype(str))]
    for j, lb in enumerate(present):
        y = 1.0 - j * 0.09
        ax_leg.scatter(0.05, y, s=18, c=lp[lb]["color"].reshape(1, -1),
                       transform=ax_leg.transAxes, clip_on=False, zorder=5)
        ax_leg.text(0.18, y, lb, fontsize=5, va="center",
                    transform=ax_leg.transAxes, clip_on=False)

    if title:
        fig.suptitle(title, fontsize=7, fontweight="extra bold", y=1.02)

    if save:
        fig.savefig(save, dpi=dpi, transparent=True, bbox_inches="tight")
    return fig


# ---------- Load data (slim — no expression matrix) ----------
print("Loading D10 Lapa (slim)...")
adata = load_slim(
    CELLASSIGN_DIR / "d10_lapa_predictions.h5ad",
    DONOR_KEY, REAL_DONORS,
)
print(f"D10 Lapa: {adata.n_obs:,} cells, donors: {REAL_DONORS}")
print(adata.obs.groupby(DONOR_KEY, observed=True)[LABEL_KEY].value_counts().unstack(fill_value=0))
print()

print("Loading EEC subset (slim)...")
adata_eec = load_slim(
    LABELLED_DIR / "knn_EECs_egfDuod_D10_Lapa_DZ.h5ad",
    DONOR_KEY, REAL_DONORS,
)
print(f"EEC subset: {adata_eec.n_obs:,} cells")
print(adata_eec.obs[EEC_LABEL_KEY].value_counts())
print()
gc.collect()

# ---------- 1. Secretory PC highlight panels ----------
print("=== Secretory PC highlights ===")
for donor in REAL_DONORS:
    ad = adata[adata.obs[DONOR_KEY] == donor].copy()
    fig = plot_veres_panel_highlight(
        ad, label_key=LABEL_KEY, highlight_label="Secretory PCs",
        stage_text=f"D10 Lapa — {donor}", palette=celltype_palette,
        save=SAVE_DIR / f"secretory_pcs_{donor}.pdf", dpi=DPI,
    )
    plt.close(fig)
    print(f"  {donor}: {ad.n_obs:,} cells")

# ---------- 2. Enterocyte highlight panels ----------
print("=== Enterocyte highlights ===")
for donor in REAL_DONORS:
    ad = adata[adata.obs[DONOR_KEY] == donor].copy()
    fig = plot_veres_panel_highlight(
        ad, label_key=LABEL_KEY, highlight_label="Enterocytes",
        stage_text=f"D10 Lapa — {donor}", palette=celltype_palette,
        save=SAVE_DIR / f"enterocytes_{donor}.pdf", dpi=DPI,
    )
    plt.close(fig)
    print(f"  {donor}: {ad.n_obs:,} cells")

# ---------- 3. EEC by donor (Veres) ----------
print("=== EEC by donor ===")
fig = plot_veres_panel(
    adata_eec, label_key=DONOR_KEY, stage_text="D10 EECs — by donor",
    label_order=REAL_DONORS, ratio_order=REAL_DONORS,
    palette=donor_palette, save=SAVE_DIR / "eec_by_donor.pdf", dpi=DPI,
)
plt.close(fig)
print("  done")

# ---------- 4. Overall D10 — all donors UMAP (scanpy) ----------
print("=== Overall D10 — all donors (scanpy) ===")
sc.pl.umap(
    adata, color=DONOR_KEY, palette=donor_palette,
    title="D10 Lapa — All donors", frameon=False, show=False,
    save="_d10_all_donors.pdf",
)
plt.close("all")
print("  done")

# ---------- 5. Overall D10 — all donors UMAP (Veres) ----------
print("=== Overall D10 — all donors (Veres) ===")
fig = plot_veres_panel(
    adata, label_key=DONOR_KEY, stage_text="D10 Lapa — All donors",
    label_order=REAL_DONORS, ratio_order=REAL_DONORS,
    palette=donor_palette, save=SAVE_DIR / "d10_all_donors_veres.pdf", dpi=DPI,
)
plt.close(fig)
print("  done")

# ---------- 6. Overall D10 — per-donor cell types (scanpy) ----------
print("=== Overall D10 — per-donor cell types (scanpy) ===")
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
for ax, donor in zip(axes, REAL_DONORS):
    ad = adata[adata.obs[DONOR_KEY] == donor].copy()
    sc.pl.umap(
        ad, color=LABEL_KEY, palette=celltype_palette,
        title=f"D10 Lapa — {donor} ({ad.n_obs:,} cells)",
        frameon=False, show=False, ax=ax,
    )
plt.tight_layout()
plt.savefig(SAVE_DIR / "d10_per_donor_celltypes_scanpy.pdf", dpi=DPI, bbox_inches="tight")
plt.close(fig)
print("  done")

# ---------- 7. Overall D10 — per-donor cell types (Veres) ----------
print("=== Overall D10 — per-donor cell types (Veres) ===")
for donor in REAL_DONORS:
    ad = adata[adata.obs[DONOR_KEY] == donor].copy()
    fig = plot_veres_panel(
        ad, label_key=LABEL_KEY, stage_text=f"D10 Lapa — {donor}",
        label_order=LABEL_ORDER, ratio_order=LABEL_ORDER,
        palette=celltype_palette, save=SAVE_DIR / f"d10_celltypes_{donor}_veres.pdf", dpi=DPI,
    )
    plt.close(fig)
    print(f"  {donor}")

# ---------- 7b. Overall D10 — per-donor cell types (Veres triptych) ----------
print("=== Overall D10 — per-donor cell types (Veres triptych) ===")
fig = plot_veres_triptych(
    adata, split_key=DONOR_KEY, donors=REAL_DONORS,
    label_key=LABEL_KEY, label_order=LABEL_ORDER,
    palette=celltype_palette, title="D10 Lapa — Cell types by donor",
    save=SAVE_DIR / "d10_per_donor_celltypes_triptych.pdf", dpi=DPI,
)
plt.close(fig)
print("  done")

del adata; gc.collect()

# ---------- 8. EEC per-donor subtypes (scanpy) ----------
print("=== EEC per-donor subtypes (scanpy) ===")
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
for ax, donor in zip(axes, REAL_DONORS):
    ad = adata_eec[adata_eec.obs[DONOR_KEY] == donor].copy()
    sc.pl.umap(
        ad, color=EEC_LABEL_KEY, palette=eec_palette,
        title=f"EEC subtypes — {donor} ({ad.n_obs:,} cells)",
        frameon=False, show=False, ax=ax,
    )
plt.tight_layout()
plt.savefig(SAVE_DIR / "eec_per_donor_subtypes_scanpy.pdf", dpi=DPI, bbox_inches="tight")
plt.close(fig)
print("  done")

# ---------- 9. EEC per-donor subtypes (Veres) ----------
print("=== EEC per-donor subtypes (Veres) ===")
for donor in REAL_DONORS:
    ad = adata_eec[adata_eec.obs[DONOR_KEY] == donor].copy()
    fig = plot_veres_panel(
        ad, label_key=EEC_LABEL_KEY, stage_text=f"D10 EECs — {donor}",
        label_order=EEC_LABEL_ORDER, ratio_order=EEC_LABEL_ORDER,
        palette=eec_palette, save=SAVE_DIR / f"eec_subtypes_{donor}_veres.pdf", dpi=DPI,
    )
    plt.close(fig)
    print(f"  {donor}: {ad.n_obs:,} EECs")

# ---------- 9b. EEC per-donor subtypes (Veres triptych) ----------
print("=== EEC per-donor subtypes (Veres triptych) ===")
fig = plot_veres_triptych(
    adata_eec, split_key=DONOR_KEY, donors=REAL_DONORS,
    label_key=EEC_LABEL_KEY, label_order=EEC_LABEL_ORDER,
    palette=eec_palette, title="D10 EEC subtypes by donor",
    save=SAVE_DIR / "eec_per_donor_subtypes_triptych.pdf", dpi=DPI,
)
plt.close(fig)
print("  done")

# ---------- Summary ----------
print()
saved = sorted(SAVE_DIR.glob("*.pdf"))
print(f"Saved {len(saved)} figures to {SAVE_DIR}:")
for f in saved:
    print(f"  {f.name}")
