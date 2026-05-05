"""
scArches-style transfer learning utilities.

Designed for memory-constrained execution (~16 GB RAM): the helpers stream
one dataset at a time during preprocessing, train SCVI/SCANVI on the D10
reference alone, then map each D2/D4 query independently against the frozen
reference model.

Mirrors the shape of `src/cell_assign.py` (compute -> run -> annotate).
"""

from __future__ import annotations

import gc
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse

from src.config import DATASETS


# ── Per-dataset preprocessing ──────────────────────────────────────────────


def _resolve_dataset_path(dataset_key: str, prefer: str = "auto") -> Path:
    """Resolve a dataset key to its h5ad path.

    `prefer="auto"` picks `labelled` when present (gives access to existing
    label columns like `initial_cellassign_prediction`), else `clustered`.
    """
    entry = DATASETS[dataset_key]
    if prefer in entry:
        return Path(entry[prefer])
    if prefer == "auto":
        for stage in ("labelled", "clustered", "raw"):
            if stage in entry and Path(entry[stage]).exists():
                return Path(entry[stage])
    raise FileNotFoundError(f"No h5ad available for dataset_key={dataset_key}")


def prepare_slim_h5ad(
    dataset_key: str,
    out_path: str | Path,
    gene_list: Sequence[str] | None = None,
    layer: str = "counts",
    prefer_stage: str = "auto",
    keep_obs_cols: Sequence[str] | None = None,
) -> Path:
    """Load a single dataset, slim it, write to disk.

    Drops `.X` (we re-derive normalisation downstream), drops obsm/obsp
    other than `X_umap` (preserved as `X_umap_orig`), and ensures
    `.layers[layer]` is sparse CSR. If `gene_list` is given, subsets to
    those genes (zero-fill is NOT done here; that's the scArches step's
    job via prepare_query_anndata).

    Parameters
    ----------
    dataset_key : str
        Key into `src.config.DATASETS`.
    out_path : Path
        Where to write the slim h5ad.
    gene_list : list[str], optional
        If set, subset `var_names` to this set (intersection).
    layer : str
        Counts layer name. Required to be present.
    prefer_stage : str
        Which DATASETS entry to read ("clustered", "labelled", "auto").
    keep_obs_cols : list[str], optional
        Which obs columns to keep. Defaults to a small standard set plus
        any cellassign-related columns from the source.
    """
    src_path = _resolve_dataset_path(dataset_key, prefer=prefer_stage)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(src_path)

    if layer not in adata.layers:
        raise KeyError(
            f"{dataset_key}: .layers['{layer}'] missing in {src_path}; cannot proceed."
        )

    # Default obs to keep
    default_keep = {
        "orig.ident", "nCount_RNA", "nFeature_RNA", "freemuxlet.identity",
        "participant", "Condition", "Time_point", "Treatment", "ident",
        "leiden", "doublet_score", "predicted_doublet",
        "initial_cellassign_prediction", "cellassign_confidence",
        "manual_label",
    }
    if keep_obs_cols is not None:
        keep = set(keep_obs_cols) & set(adata.obs.columns)
    else:
        keep = default_keep & set(adata.obs.columns)
    drop_obs = [c for c in adata.obs.columns if c not in keep]
    adata.obs = adata.obs.drop(columns=drop_obs)

    # Preserve original UMAP, drop other obsm/obsp
    if "X_umap" in adata.obsm:
        adata.obsm["X_umap_orig"] = adata.obsm["X_umap"]
    for key in list(adata.obsm.keys()):
        if key != "X_umap_orig":
            del adata.obsm[key]
    for key in list(adata.obsp.keys()):
        del adata.obsp[key]
    for key in list(adata.uns.keys()):
        del adata.uns[key]

    # Subset to gene list if provided
    if gene_list is not None:
        keep_genes = [g for g in gene_list if g in adata.var_names]
        missing = len(gene_list) - len(keep_genes)
        if missing:
            print(f"[{dataset_key}] {missing} genes from gene_list not present; kept {len(keep_genes)}")
        adata = adata[:, keep_genes].copy()

    # Force sparse CSR counts; drop other layers and .X
    counts = adata.layers[layer]
    if not sparse.isspmatrix_csr(counts):
        counts = sparse.csr_matrix(counts)
    new = ad.AnnData(
        X=counts,                    # placeholder; consumers should use .layers[layer]
        obs=adata.obs.copy(),
        var=adata.var[[]].copy(),    # drop var columns; HVG is recomputed downstream
        obsm={k: v for k, v in adata.obsm.items()},
    )
    new.layers[layer] = counts
    new.obs["dataset"] = dataset_key

    new.write_h5ad(out_path, compression="gzip")
    print(f"[{dataset_key}] wrote {out_path} -- {new.n_obs} cells x {new.n_vars} genes")

    del adata, new, counts
    gc.collect()
    return out_path


def intersect_var_names(dataset_keys: Iterable[str], prefer_stage: str = "auto") -> list[str]:
    """Read just `var_names` from each dataset (backed mode) and intersect.

    No counts loaded into RAM.
    """
    common: set[str] | None = None
    for key in dataset_keys:
        path = _resolve_dataset_path(key, prefer=prefer_stage)
        backed = sc.read_h5ad(path, backed="r")
        names = set(backed.var_names.tolist())
        backed.file.close()
        common = names if common is None else (common & names)
        print(f"[{key}] {len(names)} genes")
    out = sorted(common) if common else []
    print(f"intersection: {len(out)} genes")
    return out


# ── Reference HVG selection & training ─────────────────────────────────────


def select_reference_hvgs(
    reference_slim_path: str | Path,
    n_top_genes: int = 3000,
    flavor: str = "seurat_v3",
    layer: str = "counts",
) -> list[str]:
    """Run HVG on the reference slim file. Returns gene list.

    Uses raw counts via `flavor="seurat_v3"` (the scvi-tools recommended
    default). Loads only the reference into RAM.
    """
    adata = sc.read_h5ad(reference_slim_path)
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        flavor=flavor,
        layer=layer,
        subset=False,
    )
    hvgs = adata.var_names[adata.var["highly_variable"]].tolist()
    print(f"selected {len(hvgs)} HVGs from {reference_slim_path}")
    del adata
    gc.collect()
    return hvgs


def _prepare_reference_for_scanvi(
    adata: ad.AnnData,
    label_col: str,
    unknown_label: str,
    confidence_col: str | None,
    confidence_threshold: float | None,
) -> ad.AnnData:
    """Build the supervised label column used by SCANVI."""
    labels = adata.obs[label_col].astype(str).copy()
    if confidence_threshold is not None and confidence_col is not None:
        if confidence_col in adata.obs:
            below = adata.obs[confidence_col] < confidence_threshold
            labels[below] = unknown_label
            print(f"masked {int(below.sum())} reference cells below confidence {confidence_threshold}")
        else:
            print(f"warning: confidence column '{confidence_col}' not found; skipping threshold")
    adata.obs["scanvi_label_input"] = labels.astype("category")
    return adata


def train_reference(
    reference_slim_path: str | Path,
    model_dir: str | Path,
    label_col: str = "initial_cellassign_prediction",
    batch_key: str = "participant",
    categorical_covariates: Sequence[str] = ("Treatment",),
    n_latent: int = 30,
    n_layers: int = 2,
    layer: str = "counts",
    unknown_label: str = "Unknown",
    confidence_col: str | None = "cellassign_confidence",
    reference_confidence_threshold: float | None = None,
    max_epochs_scvi: int | None = None,
    max_epochs_scanvi: int | None = 20,
    seed: int = 0,
    accelerator: str = "auto",
) -> dict:
    """Train SCVI then SCANVI on the reference and persist both.

    Parameters
    ----------
    accelerator : str
        Forwarded to scvi-tools' `train()`. Use "gpu" on HPC.

    Returns
    -------
    dict with paths {"scvi": ..., "scanvi": ..., "reference_h5ad": ...}.
    """
    import scvi
    from scvi.model import SCVI, SCANVI

    scvi.settings.seed = seed
    model_dir = Path(model_dir)
    model_dir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(reference_slim_path)

    # Build supervised label column on the reference
    adata = _prepare_reference_for_scanvi(
        adata,
        label_col=label_col,
        unknown_label=unknown_label,
        confidence_col=confidence_col,
        confidence_threshold=reference_confidence_threshold,
    )

    # SCVI setup
    setup_kwargs = dict(layer=layer, batch_key=batch_key)
    cat_covs = [c for c in categorical_covariates if c in adata.obs.columns]
    if cat_covs:
        setup_kwargs["categorical_covariate_keys"] = cat_covs
    SCVI.setup_anndata(adata, **setup_kwargs)

    scvi_model = SCVI(adata, n_latent=n_latent, n_layers=n_layers)
    scvi_model.train(max_epochs=max_epochs_scvi, accelerator=accelerator)
    scvi_model.save(model_dir / "scvi", overwrite=True)

    # Convert to SCANVI
    scanvi_model = SCANVI.from_scvi_model(
        scvi_model,
        adata=adata,
        labels_key="scanvi_label_input",
        unlabeled_category=unknown_label,
    )
    scanvi_model.train(max_epochs=max_epochs_scanvi, accelerator=accelerator)
    scanvi_model.save(model_dir / "scanvi", overwrite=True)

    # Persist the reference adata used during training (latent + label input)
    adata.obsm["X_scVI"] = scvi_model.get_latent_representation()
    adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation()
    adata.obs["scanvi_prediction_ref"] = scanvi_model.predict()

    ref_h5ad = model_dir / "reference_with_latent.h5ad"
    adata.write_h5ad(ref_h5ad, compression="gzip")

    paths = {
        "scvi": str(model_dir / "scvi"),
        "scanvi": str(model_dir / "scanvi"),
        "reference_h5ad": str(ref_h5ad),
    }
    print("reference trained:", paths)

    del scvi_model, scanvi_model, adata
    gc.collect()
    return paths


# ── Per-query mapping ──────────────────────────────────────────────────────


def _annotate_query_obs(
    adata: ad.AnnData,
    soft: pd.DataFrame,
    ambiguous_threshold: float,
    proba_prefix: str,
) -> ad.AnnData:
    """Write argmax, raw, ambiguous-flagged label, confidence, and per-class
    probability columns into `adata.obs`."""
    proba = soft.reindex(adata.obs_names)
    raw = proba.idxmax(axis=1).astype(str)
    conf = proba.max(axis=1).astype(float)

    flagged = raw.copy()
    if ambiguous_threshold is not None:
        below = conf < ambiguous_threshold
        flagged[below] = "Ambiguous_" + raw[below]

    adata.obs["scanvi_prediction_raw"] = raw.values
    adata.obs["scanvi_prediction"] = flagged.values
    adata.obs["scanvi_confidence"] = conf.values
    for col in proba.columns:
        adata.obs[f"{proba_prefix}{col}"] = proba[col].to_numpy()
    return adata


def map_query(
    query_slim_path: str | Path,
    ref_model_dir: str | Path,
    predictions_h5ad_path: str | Path,
    predictions_csv_path: str | Path,
    layer: str = "counts",
    max_epochs_query: int = 100,
    ambiguous_threshold: float = 0.5,
    proba_prefix: str = "scanvi_prob_",
    seed: int = 0,
    accelerator: str = "auto",
) -> dict:
    """Map one query against the frozen SCANVI reference (scArches surgery).

    Returns dict with {"h5ad": ..., "csv": ..., "n_missing_genes": ...}.
    """
    import scvi
    from scvi.model import SCANVI

    scvi.settings.seed = seed
    ref_scanvi_dir = Path(ref_model_dir) / "scanvi"

    predictions_h5ad_path = Path(predictions_h5ad_path)
    predictions_csv_path = Path(predictions_csv_path)
    predictions_h5ad_path.parent.mkdir(parents=True, exist_ok=True)
    predictions_csv_path.parent.mkdir(parents=True, exist_ok=True)

    query = sc.read_h5ad(query_slim_path)

    # SCANVI requires the same setup keys; queries get "Unknown" labels
    # and may have new levels of categorical covariates.
    if "scanvi_label_input" not in query.obs:
        query.obs["scanvi_label_input"] = "Unknown"
    query.obs["scanvi_label_input"] = query.obs["scanvi_label_input"].astype("category")

    # Harmonise gene set against the reference (zero-fill missing genes)
    SCANVI.prepare_query_anndata(query, str(ref_scanvi_dir))
    n_missing = int((query.X.sum(axis=0) == 0).sum()) if hasattr(query.X, "sum") else 0

    # Surgery: load query data using the reference architecture, then a
    # short fine-tune so the latent space accommodates the query batch
    qmodel = SCANVI.load_query_data(query, str(ref_scanvi_dir))
    qmodel.train(
        max_epochs=max_epochs_query,
        plan_kwargs={"weight_decay": 0.0},
        accelerator=accelerator,
    )

    # Predict
    soft = qmodel.predict(soft=True)
    if not isinstance(soft, pd.DataFrame):
        soft = pd.DataFrame(soft, index=query.obs_names)
    query = _annotate_query_obs(
        query,
        soft=soft,
        ambiguous_threshold=ambiguous_threshold,
        proba_prefix=proba_prefix,
    )

    # Latent + UMAP on the integrated latent
    query.obsm["X_scANVI"] = qmodel.get_latent_representation()
    sc.pp.neighbors(query, use_rep="X_scANVI", n_neighbors=30)
    sc.tl.umap(query)
    query.obsm["X_umap_scvi"] = query.obsm["X_umap"]
    if "X_umap_orig" in query.obsm:
        query.obsm["X_umap"] = query.obsm["X_umap_orig"]

    # Save outputs
    query.write_h5ad(predictions_h5ad_path, compression="gzip")
    soft.to_csv(predictions_csv_path, index=True)
    print(f"wrote {predictions_h5ad_path} and {predictions_csv_path}")

    out = {
        "h5ad": str(predictions_h5ad_path),
        "csv": str(predictions_csv_path),
        "n_missing_genes": n_missing,
    }
    del qmodel, query, soft
    gc.collect()
    return out


# ── Lightweight aggregation for plotting ───────────────────────────────────


def combine_predictions_obs(
    query_h5ad_paths: Sequence[str | Path],
    reference_h5ad_path: str | Path | None = None,
    out_csv: str | Path | None = None,
    obs_cols: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Read `obs` only from each prediction h5ad and concat into a thin
    DataFrame for plotting/population-ratio analyses.

    Counts/.X are NOT loaded.
    """
    default_cols = [
        "dataset", "participant", "Time_point", "Treatment",
        "scanvi_prediction", "scanvi_prediction_raw", "scanvi_confidence",
        "initial_cellassign_prediction", "cellassign_confidence",
    ]
    cols = list(obs_cols) if obs_cols else default_cols

    frames: list[pd.DataFrame] = []
    paths = list(query_h5ad_paths)
    if reference_h5ad_path is not None:
        paths = [reference_h5ad_path] + paths

    for p in paths:
        backed = sc.read_h5ad(p, backed="r")
        present = [c for c in cols if c in backed.obs.columns]
        # Read obs without loading counts
        df = backed.obs[present].copy()
        # Pull UMAP coords from obsm if present (cheap)
        for key in ("X_umap_scvi", "X_umap_orig", "X_umap"):
            if key in backed.obsm:
                arr = backed.obsm[key][:]
                df[f"{key}_0"] = arr[:, 0]
                df[f"{key}_1"] = arr[:, 1]
                break
        backed.file.close()
        frames.append(df)

    out = pd.concat(frames, axis=0)
    if out_csv is not None:
        Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(out_csv, index=True)
        print(f"wrote {out_csv} ({len(out)} cells)")
    return out
