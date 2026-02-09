"""
Decoupler TF activity inference with HVG filtering and memory management.

Usage from a notebook:
    from src.decoupler import run_decoupler_safe
    adata, score = run_decoupler_safe(adata, celltype_key="manual_label")
"""

import gc

import numpy as np
import scanpy as sc
import decoupler as dc

from src.palette import celltype_palette


def run_decoupler_safe(
    adata,
    celltype_key="manual_label",
    tmin=25,
    use_hvg=True,
    max_genes=3000,
):
    """
    Run decoupler ULM with HVG filtering and memory management.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix (modified in-place with scores in obsm).
    celltype_key : str
        Column in ``adata.obs`` used for cell type colours.
    tmin : int
        Minimum number of targets per TF to keep.
    use_hvg : bool
        If True, subset to highly variable genes before scoring.
    max_genes : int
        Number of top HVGs to retain when ``use_hvg=True``.

    Returns
    -------
    adata : AnnData
        Input object (potentially HVG-subsetted) with ``score_ulm``
        and ``padj_ulm`` added to ``.obsm``.
    score_adata : AnnData
        Cells x TFs AnnData extracted from ``adata.obsm['score_ulm']``.
    """
    print(f"  Original: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # ── 1. Optional HVG subsetting ──
    if use_hvg:
        if "highly_variable" not in adata.var.columns:
            print("  Computing HVGs...")
            if "counts" in adata.layers:
                adata_tmp = adata.copy()
                adata_tmp.X = adata_tmp.layers["counts"]
                sc.pp.highly_variable_genes(
                    adata_tmp, n_top_genes=max_genes, flavor="seurat_v3"
                )
                adata.var["highly_variable"] = adata_tmp.var["highly_variable"]
                del adata_tmp
                gc.collect()
            else:
                sc.pp.highly_variable_genes(adata, n_top_genes=max_genes)

        n_hvg = adata.var["highly_variable"].sum()
        print(f"  Subsetting to {n_hvg:,} HVGs")
        adata = adata[:, adata.var["highly_variable"]].copy()
        gc.collect()

    # ── 2. Cell type colours ──
    if celltype_key and celltype_key in adata.obs.columns:
        adata.obs[celltype_key] = adata.obs[celltype_key].astype("category")
        categories = adata.obs[celltype_key].cat.categories
        adata.uns[f"{celltype_key}_colors"] = [
            celltype_palette.get(ct, "#808080") for ct in categories
        ]

    # ── 3. Load CollecTRI & run ULM ──
    print("  Loading CollecTRI...")
    collectri = dc.op.collectri(organism="human")
    collectri_sub = collectri[collectri["target"].isin(adata.var_names)].copy()
    print(f"  {collectri_sub['source'].nunique()} TFs, "
          f"{len(collectri_sub):,} interactions after intersecting")

    print(f"  Running ULM (tmin={tmin})...")
    dc.mt.ulm(data=adata, net=collectri_sub, tmin=tmin, verbose=True)

    # ── 4. Extract score AnnData ──
    score_adata = dc.pp.get_obsm(adata, key="score_ulm")
    n_tfs = score_adata.n_vars
    print(f"  Computed activity for {n_tfs} TFs")

    del collectri, collectri_sub
    gc.collect()

    return adata, score_adata
