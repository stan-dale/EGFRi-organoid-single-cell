#!/usr/bin/env python
"""Step 1 (unified labelling) — build ONE scVI latent per timepoint.

Subsets the canonical integrated.h5ad to a single timepoint (control + Lapa
joined), batch-corrects on participant with scVI, then clusters the integrated
latent (Leiden) and computes a UMAP. The saved object keeps the FULL gene set
(scVI trains on HVGs only) so downstream marker scoring in step 2 can reach any
gene.

This is the new "manual labelling primary" approach — there is NO reference /
query split and NO label transfer here. The existing 4-class scANVI predictions
are merged in step 2 as an advisory prior only.

Timepoint -> datasets:
    D2 : D2_DZ + D2_Lapa
    D4 : D4_DZ + D4_Lapa   (D4_AS is excluded — control + Lapa only)

Usage:
    $HOME/venvs/single-cell-gpu/bin/python 20_integrate_timepoint.py D2
    TIMEPOINT=D4 sbatch 20_integrate_timepoint.sbatch

Outputs:
    ~/single-cell/analysis/integration/unified/d{2,4}_unified.h5ad
        obsm["X_scVI"], obs["leiden"], obsm["X_umap"], layers["counts"], log-norm X
"""
import os
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc
import scvi
import torch

# ── config ───────────────────────────────────────────────────────────────
HOME = Path.home()
DATA_PATH = HOME / "single-cell" / "data" / "integrated.h5ad"
OUT_DIR = HOME / "single-cell" / "analysis" / "integration" / "unified"

# control + Lapa only, per timepoint (D4_AS deliberately excluded)
TIMEPOINT_DATASETS = {
    "D2": ["D2_DZ", "D2_Lapa"],
    "D4": ["D4_DZ", "D4_Lapa"],
}

N_HVG = int(os.environ.get("N_HVG", "3000"))
N_LATENT = int(os.environ.get("N_LATENT", "30"))
N_LAYERS = int(os.environ.get("N_LAYERS", "2"))
SCVI_EPOCHS = int(os.environ.get("SCVI_EPOCHS", "200"))
BATCH_KEY = os.environ.get("BATCH_KEY", "participant")
LEIDEN_RES = float(os.environ.get("LEIDEN_RES", "1.0"))
N_NEIGHBORS = int(os.environ.get("N_NEIGHBORS", "30"))
SEED = int(os.environ.get("SEED", "0"))


def banner(msg: str):
    print(f"\n{'=' * 8}  {msg}  {'=' * 8}", flush=True)


def run_leiden(adata, resolution, key="leiden"):
    """Leiden on the existing neighbour graph, preferring the fast igraph
    flavour but falling back to the legacy default if igraph is unavailable."""
    try:
        sc.tl.leiden(
            adata,
            resolution=resolution,
            key_added=key,
            flavor="igraph",
            n_iterations=2,
            directed=False,
        )
    except (TypeError, ImportError) as e:
        print(f"  igraph leiden unavailable ({e}); using default flavour", flush=True)
        sc.tl.leiden(adata, resolution=resolution, key_added=key)


def main():
    timepoint = (sys.argv[1] if len(sys.argv) > 1
                 else os.environ.get("TIMEPOINT", "")).upper()
    if timepoint not in TIMEPOINT_DATASETS:
        sys.exit(f"timepoint must be one of {list(TIMEPOINT_DATASETS)}; got {timepoint!r}")
    keep_datasets = TIMEPOINT_DATASETS[timepoint]

    scvi.settings.seed = SEED
    accelerator = "gpu" if torch.cuda.is_available() else "cpu"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / f"{timepoint.lower()}_unified.h5ad"

    print(f"scvi={scvi.__version__}  torch={torch.__version__}  accelerator={accelerator}")
    print(f"timepoint={timepoint}  datasets={keep_datasets}")
    print(f"data={DATA_PATH}")
    print(f"out={out_path}")
    print(f"params: N_HVG={N_HVG} N_LATENT={N_LATENT} N_LAYERS={N_LAYERS} "
          f"SCVI_EPOCHS={SCVI_EPOCHS} BATCH_KEY={BATCH_KEY} "
          f"LEIDEN_RES={LEIDEN_RES} N_NEIGHBORS={N_NEIGHBORS} SEED={SEED}")

    banner("load integrated.h5ad")
    t0 = time.time()
    a = ad.read_h5ad(DATA_PATH)
    print(f"loaded {a.shape} in {time.time() - t0:.1f}s  (X dtype={a.X.dtype})")

    banner(f"subset to {timepoint} ({' + '.join(keep_datasets)})")
    adata = a[a.obs["dataset"].isin(keep_datasets)].copy()
    del a
    print(f"subset: {adata.shape}")
    print("dataset counts:")
    print(adata.obs["dataset"].value_counts().to_string())
    print(f"\n{BATCH_KEY} counts:")
    print(adata.obs[BATCH_KEY].value_counts().to_string())
    if "counts" not in adata.layers:
        sys.exit("ERROR: layers['counts'] missing — required for scVI.")

    banner(f"select top {N_HVG} HVGs (seurat_v3, counts, batch-aware on {BATCH_KEY})")
    t0 = time.time()
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=N_HVG,
        flavor="seurat_v3",
        layer="counts",
        batch_key=BATCH_KEY,
        subset=False,          # keep full gene set in the saved object
    )
    n_hvg = int(adata.var["highly_variable"].sum())
    print(f"HVGs flagged: {n_hvg}  ({time.time() - t0:.1f}s)")

    banner("train scVI on the HVG subset")
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    scvi.model.SCVI.setup_anndata(adata_hvg, layer="counts", batch_key=BATCH_KEY)
    vae = scvi.model.SCVI(adata_hvg, n_latent=N_LATENT, n_layers=N_LAYERS)
    t0 = time.time()
    vae.train(
        max_epochs=SCVI_EPOCHS,
        accelerator=accelerator,
        devices=1,
        train_size=0.9,
        check_val_every_n_epoch=1,
    )
    print(f"scVI train took {time.time() - t0:.1f}s")

    banner("attach latent to full object")
    adata.obsm["X_scVI"] = vae.get_latent_representation()
    print(f"X_scVI: {adata.obsm['X_scVI'].shape}")

    banner(f"neighbours (use_rep=X_scVI, n_neighbors={N_NEIGHBORS}) + Leiden (res={LEIDEN_RES}) + UMAP")
    t0 = time.time()
    sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=N_NEIGHBORS, random_state=SEED)
    run_leiden(adata, resolution=LEIDEN_RES, key="leiden")
    sc.tl.umap(adata, random_state=SEED)
    print(f"graph/cluster/umap took {time.time() - t0:.1f}s")
    print("leiden cluster sizes:")
    print(adata.obs["leiden"].value_counts().sort_index().to_string())

    import pandas as pd
    banner("leiden × dataset cross-tab (treatment composition per cluster)")
    print(pd.crosstab(adata.obs["leiden"], adata.obs["dataset"]).to_string())

    banner("leiden × participant cross-tab (BATCH-CORRECTION diagnostic)")
    # If donors mix within clusters, scVI integrated them. A cluster that is
    # ~100% one donor is under-corrected. NB: H439 is Lapa-only and H897 is
    # DZ-only at D2, so a donor-pure cluster can still be real biology.
    ct_part = pd.crosstab(adata.obs["leiden"], adata.obs[BATCH_KEY])
    print(ct_part.to_string())
    # Per-cluster dominant-donor fraction — quick "is anything donor-pure?" scan
    frac = ct_part.div(ct_part.sum(axis=1), axis=0)
    print("\nmax single-donor fraction per cluster (closer to 1.0 = less mixed):")
    print(frac.max(axis=1).round(2).to_string())

    banner(f"write {out_path}")
    adata.uns["unified_integration"] = {
        "timepoint": timepoint,
        "datasets": keep_datasets,
        "n_hvg": n_hvg,
        "n_latent": N_LATENT,
        "n_layers": N_LAYERS,
        "scvi_epochs": SCVI_EPOCHS,
        "batch_key": BATCH_KEY,
        "leiden_resolution": LEIDEN_RES,
        "n_neighbors": N_NEIGHBORS,
        "seed": SEED,
    }
    adata.write_h5ad(out_path, compression="gzip")
    print(f"wrote {adata.shape} -> {out_path}")
    print("OK")


if __name__ == "__main__":
    main()
