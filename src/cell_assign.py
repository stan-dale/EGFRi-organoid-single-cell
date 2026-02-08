"""
CellAssign wrapper functions extracted from the project's
4_Cell_Assign.ipynb notebooks.

Identical workflow is used across 6 datasets.
"""

import numpy as np
import pandas as pd
from scipy import sparse
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns


def compute_size_factors(adata, layer="counts"):
    """
    Compute library-size normalisation factors for CellAssign.

    Adds 'size_factor' to adata.obs.
    """
    lib_size = adata.layers[layer].sum(1)
    if sparse.issparse(lib_size):
        lib_size = np.asarray(lib_size).flatten()
    elif hasattr(lib_size, "A1"):
        lib_size = lib_size.A1
    else:
        lib_size = np.asarray(lib_size).flatten()
    adata.obs["size_factor"] = lib_size / np.mean(lib_size)
    return adata


def run_cellassign(adata, marker_csv_path, layer="counts", seed=0, drop_genes=None):
    """
    Run CellAssign on *adata* using a binary marker matrix CSV.

    Parameters
    ----------
    adata : AnnData
        Must have size_factor in .obs and *layer* in .layers.
    marker_csv_path : str or Path
        CSV with genes as rows, cell types as columns, 0/1 values.
    layer : str
    seed : int
    drop_genes : list[str], optional
        Gene names to drop from the marker matrix (e.g. genes not
        present in the dataset).

    Returns
    -------
    predictions : DataFrame  (cells x cell types, probability values)
    """
    import scvi
    from scvi.external import CellAssign

    scvi.settings.seed = seed

    markers = pd.read_csv(marker_csv_path, index_col=0)
    if drop_genes:
        markers = markers.drop(index=[g for g in drop_genes if g in markers.index])

    # Subset adata to marker genes
    filtered = adata[:, markers.index].copy()

    CellAssign.setup_anndata(filtered, layer=layer, size_factor_key="size_factor")
    model = CellAssign(filtered, markers)
    model.train()

    predictions = model.predict()
    return predictions


def annotate_predictions(adata, predictions, key="initial_cellassign_prediction"):
    """
    Store CellAssign predictions in adata.obs.

    Adds per-cell-type probability columns, the argmax prediction,
    and the confidence score.
    """
    for col in predictions.columns:
        adata.obs[f"cell_assign_prob_{col}"] = predictions[col].to_numpy()
    adata.obs[key] = predictions.idxmax(axis=1).to_numpy()
    adata.obs["cellassign_confidence"] = predictions.max(axis=1).to_numpy()
    return adata


def plot_cluster_assignment_heatmap(adata, predictions, cluster_key="leiden"):
    """
    Heatmap of mean CellAssign probability per Leiden cluster.
    """
    prob_cols = [c for c in adata.obs.columns if c.startswith("cell_assign_prob_")]
    cluster_means = adata.obs.groupby(cluster_key)[prob_cols].mean()
    cluster_means.columns = [c.replace("cell_assign_prob_", "") for c in cluster_means.columns]

    fig, ax = plt.subplots(figsize=(10, 6))
    sns.heatmap(cluster_means, cmap="viridis", annot=True, fmt=".2f", ax=ax)
    ax.set_title("Mean CellAssign Probability per Cluster")
    ax.set_xlabel("CellAssign cell type")
    ax.set_ylabel("Cluster")
    plt.tight_layout()
    plt.show()
    return fig
