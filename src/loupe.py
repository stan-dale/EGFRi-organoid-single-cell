"""
Utilities for converting AnnData objects to Loupe Browser (.cloupe) files.

Requires the ``loupe-env`` conda environment (Python >=3.11, loupepy).

Workarounds for two loupepy 1.1.1 bugs:
  1. ``get_count_matrix()`` ignores the ``layer`` parameter and always uses
     ``adata.X``.  We swap raw counts into ``.X`` before converting.
  2. The converter binary path is passed to ``os.system()`` without quoting,
     so paths containing spaces (e.g. macOS "Application Support") break.
     We accept an explicit ``converter_path`` that should point to a
     space-free symlink.
"""

import gc
import os
from pathlib import Path

import anndata as ad
import scanpy as sc
from scipy.sparse import csr_matrix

# Default converter location (space-free symlink)
_DEFAULT_CONVERTER = Path.home() / ".local" / "bin" / "loupe_converter"


def ensure_converter_symlink() -> Path:
    """Create a space-free symlink to the loupepy converter binary.

    Returns the symlink path.  Does nothing if it already exists.
    """
    target = (
        Path.home() / "Library" / "Application Support" / "loupepy" / "loupe_converter"
    )
    link = _DEFAULT_CONVERTER
    link.parent.mkdir(parents=True, exist_ok=True)
    if not link.exists() and target.exists():
        os.symlink(target, link)
    return link


def convert_h5ad_to_cloupe(
    h5ad_path,
    out_path,
    converter_path=None,
):
    """Memory-optimised h5ad → .cloupe conversion.

    Loads the full AnnData, extracts only what loupepy needs (raw counts,
    categorical obs columns, 2-D projections), frees the original object,
    then calls ``create_loupe_from_anndata``.

    Parameters
    ----------
    h5ad_path : str or Path
        Input ``.h5ad`` file.
    out_path : str or Path
        Output ``.cloupe`` file.
    converter_path : str or Path, optional
        Path to the ``loupe_converter`` binary.  Defaults to
        ``~/.local/bin/loupe_converter`` (the space-free symlink).
    """
    from loupepy import create_loupe_from_anndata

    if converter_path is None:
        converter_path = str(ensure_converter_symlink())

    h5ad_path, out_path = Path(h5ad_path), Path(out_path)

    adata = sc.read_h5ad(h5ad_path)
    print(f"  Loaded: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

    # Prefer raw counts; fall back to X
    counts = adata.layers.get("counts", adata.X)
    if not isinstance(counts, csr_matrix):
        counts = csr_matrix(counts)

    # Keep only categorical obs (loupepy drops the rest anyway)
    cat_cols = [c for c in adata.obs.columns if adata.obs[c].dtype.name == "category"]
    obs_df = adata.obs[cat_cols].copy()

    # Keep only 2-D projections (loupepy drops higher-dim)
    obsm = {
        k: v.copy()
        for k, v in adata.obsm.items()
        if v.ndim == 2 and v.shape[1] == 2
    }

    var_df = adata.var[[]].copy()
    barcodes = adata.obs.index.copy()

    del adata
    gc.collect()

    # Put counts into .X (works around loupepy layer bug)
    adata_min = ad.AnnData(X=counts, obs=obs_df, var=var_df, obsm=obsm)
    adata_min.obs.index = barcodes
    del counts, obs_df, var_df, obsm
    gc.collect()

    print(
        f"  Minimal: {len(cat_cols)} categories, "
        f"{len(adata_min.obsm)} projections"
    )

    create_loupe_from_anndata(
        adata_min,
        output_cloupe=str(out_path),
        loupe_converter_path=str(converter_path),
        force=True,
    )

    del adata_min
    gc.collect()

    size_mb = out_path.stat().st_size / 1e6
    print(f"  -> {out_path.name}  ({size_mb:.0f} MB)")
