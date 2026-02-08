"""
QC and clustering pipelines extracted from the project's
0_Quality_Control.ipynb and 1_Clustering.ipynb notebooks.

These functions wrap the standard Scanpy workflow used identically
across all 7 datasets.
"""

import scanpy as sc
import numpy as np
from pathlib import Path


# ── Quality Control ──

def annotate_qc_vars(adata):
    """Add mitochondrial and ribosomal gene annotations + QC metrics."""
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo"], inplace=True, log1p=True
    )
    return adata


def run_scrublet(adata):
    """Run scrublet doublet detection."""
    sc.pp.scrublet(adata)
    return adata


def run_qc_pipeline(adata, dataset_name, save_dir=None):
    """
    Full QC pipeline: annotate metrics, plot, run scrublet.

    Parameters
    ----------
    adata : AnnData
    dataset_name : str
        Used for plot titles and filenames.
    save_dir : Path, optional
        If provided, saves the QC-annotated object here.

    Returns
    -------
    adata : AnnData  (modified in-place and returned)
    """
    print(f"Running QC for {dataset_name}...")
    annotate_qc_vars(adata)

    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        show=False,
    )
    sc.pl.scatter(
        adata, "total_counts", "n_genes_by_counts",
        color="pct_counts_mt", show=False,
    )

    run_scrublet(adata)
    sc.pl.scrublet_score_distribution(adata, show=False)

    if save_dir is not None:
        save_dir = Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)
        out = save_dir / f"doublet_{dataset_name}.h5ad"
        adata.write_h5ad(out)
        print(f"  Saved: {out}")

    return adata


# ── Clustering ──

def run_clustering_pipeline(
    adata,
    dataset_name,
    cell_cycle_markers=None,
    n_top_genes=2000,
    n_neighbors=100,
    n_pcs=30,
    leiden_resolution=1.0,
    regress_cell_cycle=True,
    save_dir=None,
):
    """
    Full clustering pipeline: HVG → scale → cell-cycle regression →
    PCA → neighbors → Leiden → UMAP.

    Operates on a HVG subset for dimensionality reduction, then
    transfers embeddings back to the full object.

    Parameters
    ----------
    adata : AnnData
        Should already have .layers["counts"] and log-normalised .X.
    dataset_name : str
    cell_cycle_markers : dict, optional
        Must contain "S_genes" and "G2M_genes" lists.
    n_top_genes : int
    n_neighbors : int
    n_pcs : int
    leiden_resolution : float
    regress_cell_cycle : bool
    save_dir : Path, optional

    Returns
    -------
    adata : AnnData  (with X_pca, X_umap, leiden added)
    """
    from .markers import cell_cycle_markers as default_cc

    if cell_cycle_markers is None:
        cell_cycle_markers = default_cc

    print(f"Clustering {dataset_name}...")

    # HVG selection
    sc.pp.highly_variable_genes(
        adata, flavor="seurat", n_top_genes=n_top_genes
    )
    var = adata[:, adata.var["highly_variable"]].copy()

    # Scale
    sc.pp.scale(var, max_value=10)
    sc.pp.scale(adata, max_value=10)

    # Cell cycle regression
    if regress_cell_cycle:
        s_genes = [g for g in cell_cycle_markers["S_genes"] if g in var.var_names]
        g2m_genes = [g for g in cell_cycle_markers["G2M_genes"] if g in var.var_names]
        if s_genes and g2m_genes:
            sc.tl.score_genes_cell_cycle(var, s_genes=s_genes, g2m_genes=g2m_genes)
            sc.pp.regress_out(var, ["S_score", "G2M_score"])
            print(f"  Regressed cell cycle ({len(s_genes)} S, {len(g2m_genes)} G2M genes)")

    # PCA
    sc.tl.pca(var, n_comps=50, svd_solver="arpack", use_highly_variable=True)

    # Neighbors + Leiden + UMAP
    sc.pp.neighbors(var, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep="X", metric="euclidean")
    sc.tl.leiden(var, resolution=leiden_resolution, n_iterations=2)
    sc.tl.umap(var)

    # Transfer back to full object
    adata.obsm["X_pca"] = var.obsm["X_pca"]
    adata.obsm["X_umap"] = var.obsm["X_umap"]
    adata.obs["leiden"] = var.obs["leiden"].astype("category")
    adata.obsp["distances"] = var.obsp["distances"]
    adata.obsp["connectivities"] = var.obsp["connectivities"]

    if "leiden_colors" in adata.uns:
        del adata.uns["leiden_colors"]

    print(f"  {adata.obs['leiden'].nunique()} clusters at resolution {leiden_resolution}")

    if save_dir is not None:
        save_dir = Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)
        out = save_dir / f"clustered_{dataset_name}.h5ad"
        adata.write_h5ad(out)
        print(f"  Saved: {out}")

    return adata
