#!/usr/bin/env python
"""Step 1b — integration QC: UMAP panels + batch-mixing diagnostic.

Reads a unified per-timepoint object (from 20_integrate_timepoint.py) and writes
diagnostic UMAPs coloured by treatment (dataset), donor (participant), condition,
and Leiden cluster — the visuals needed to judge whether scVI integrated the two
treatment arms or whether DZ/Lapa look like unmerged batches.

Also reprints the leiden × participant cross-tab + max single-donor fraction, so
runs launched before that diagnostic was added still get it.

Reading guide:
  - participant panel  -> if donors intermix, batch correction worked.
  - dataset panel      -> DZ/Lapa sharing a core but diverging elsewhere = biology;
                          totally disjoint islands = suspicious.
  - shared donors (H329/H896 at D2) co-localising in common states is the
    cleanest evidence that the DZ/Lapa split is biological, not a run artifact.

Usage:
    $HOME/venvs/single-cell-gpu/bin/python 21_unified_umaps.py D2
"""
import os
import sys
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

HOME = Path.home()
UNIFIED_DIR = HOME / "single-cell" / "analysis" / "integration" / "unified"
FIG_DIR = HOME / "single-cell" / "figures" / "d2-d4-unified"


def banner(msg: str):
    print(f"\n{'=' * 8}  {msg}  {'=' * 8}", flush=True)


def main():
    timepoint = (sys.argv[1] if len(sys.argv) > 1
                 else os.environ.get("TIMEPOINT", "")).upper()
    if timepoint not in {"D2", "D4"}:
        sys.exit(f"timepoint must be D2 or D4; got {timepoint!r}")

    in_path = UNIFIED_DIR / f"{timepoint.lower()}_unified.h5ad"
    if not in_path.exists():
        sys.exit(f"not found: {in_path} (run 20_integrate_timepoint.py first)")
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    banner(f"load {in_path}")
    adata = ad.read_h5ad(in_path)
    print(f"{adata.shape}  obs cols: {list(adata.obs.columns)}")
    if "X_umap" not in adata.obsm:
        sys.exit("X_umap missing in obsm — integration step did not write a UMAP.")

    # ── batch-mixing diagnostic (recovered for pre-diagnostic runs) ──────────
    if {"leiden", "participant"}.issubset(adata.obs.columns):
        banner("leiden × participant (BATCH-CORRECTION diagnostic)")
        ct = pd.crosstab(adata.obs["leiden"], adata.obs["participant"])
        print(ct.to_string())
        frac = ct.div(ct.sum(axis=1), axis=0)
        print("\nmax single-donor fraction per cluster (near 1.0 = donor-pure):")
        print(frac.max(axis=1).round(2).to_string())
    if {"leiden", "dataset"}.issubset(adata.obs.columns):
        banner("leiden × dataset (treatment composition)")
        print(pd.crosstab(adata.obs["leiden"], adata.obs["dataset"]).to_string())

    # ── UMAP panels ──────────────────────────────────────────────────────────
    color_keys = [k for k in ("dataset", "participant", "condition", "leiden")
                  if k in adata.obs.columns]
    banner(f"UMAP panels: {color_keys}")

    # one combined multi-panel figure
    n = len(color_keys)
    fig, axes = plt.subplots(1, n, figsize=(6 * n, 5.5))
    if n == 1:
        axes = [axes]
    for ax, key in zip(axes, color_keys):
        legend_loc = "on data" if key == "leiden" else "right margin"
        sc.pl.umap(adata, color=key, ax=ax, show=False, size=4,
                   legend_loc=legend_loc, title=f"{timepoint} — {key}")
    fig.tight_layout()
    combined = FIG_DIR / f"{timepoint.lower()}_umap_panels.png"
    fig.savefig(combined, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {combined}")

    # treatment-split: DZ vs Lapa side by side on the same coordinates
    if "dataset" in adata.obs.columns:
        datasets = list(adata.obs["dataset"].cat.categories) \
            if hasattr(adata.obs["dataset"], "cat") \
            else sorted(adata.obs["dataset"].unique())
        fig, axes = plt.subplots(1, len(datasets), figsize=(6 * len(datasets), 5.5),
                                 sharex=True, sharey=True)
        if len(datasets) == 1:
            axes = [axes]
        for ax, ds in zip(axes, datasets):
            mask = (adata.obs["dataset"] == ds).values
            ax.scatter(adata.obsm["X_umap"][~mask, 0], adata.obsm["X_umap"][~mask, 1],
                       s=2, c="lightgrey", linewidths=0, rasterized=True)
            ax.scatter(adata.obsm["X_umap"][mask, 0], adata.obsm["X_umap"][mask, 1],
                       s=2, c="firebrick", linewidths=0, rasterized=True)
            ax.set_title(f"{timepoint} — {ds} ({int(mask.sum())} cells)")
            ax.set_xticks([]); ax.set_yticks([])
        fig.tight_layout()
        split = FIG_DIR / f"{timepoint.lower()}_umap_by_treatment.png"
        fig.savefig(split, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"wrote {split}")

    print("OK")


if __name__ == "__main__":
    main()
