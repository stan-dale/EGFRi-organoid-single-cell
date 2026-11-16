"""QC, clustering and CellAssign helpers on synthetic data."""

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

from src.cell_assign import annotate_predictions, compute_size_factors
from src.preprocessing import annotate_qc_vars, run_clustering_pipeline


def test_qc_vars_count_mitochondrial_fraction(adata):
    annotate_qc_vars(adata)
    assert adata.var["mt"].sum() == 3
    assert adata.var["ribo"].sum() == 2
    # 3 MT genes at mean 3 vs the rest at ~1 (markers higher): a few percent, never 0 or 100.
    pct = adata.obs["pct_counts_mt"]
    assert 0 < pct.median() < 20


def test_clustering_recovers_planted_cell_types(adata, tmp_path):
    run_clustering_pipeline(adata, "synthetic", n_top_genes=60, n_neighbors=15,
                            n_pcs=10, leiden_resolution=0.3, save_dir=tmp_path)
    assert adata.obsm["X_umap"].shape == (adata.n_obs, 2)
    assert "connectivities" in adata.obsp
    ari = adjusted_rand_score(adata.obs["true_type"], adata.obs["leiden"])
    assert ari > 0.9, f"clusters do not match planted types (ARI={ari:.2f})"
    assert (tmp_path / "clustered_synthetic.h5ad").exists()


def test_size_factors_average_to_one(adata):
    compute_size_factors(adata)
    assert np.isclose(adata.obs["size_factor"].mean(), 1.0)
    assert (adata.obs["size_factor"] > 0).all()


def test_annotate_predictions_takes_argmax(adata):
    probs = pd.DataFrame({"ISCs": 0.1, "EECs": 0.9}, index=adata.obs_names)
    probs.iloc[0] = [0.8, 0.2]
    annotate_predictions(adata, probs)
    assert adata.obs["initial_cellassign_prediction"].iloc[0] == "ISCs"
    assert (adata.obs["initial_cellassign_prediction"].iloc[1:] == "EECs").all()
    assert np.isclose(adata.obs["cellassign_confidence"].iloc[0], 0.8)
    assert "cell_assign_prob_EECs" in adata.obs
