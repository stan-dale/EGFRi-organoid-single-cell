"""Shared fixtures: a small synthetic organoid dataset with known structure.

The dataset mimics what the pipeline expects from data/data-objects/raw/*.h5ad
(see r-src/build_h5ad_from_10x.py): raw integer counts in layers["counts"],
log-normalised .X, and obs columns participant / Condition.

Planted structure lets tests check answers, not just that code runs:
  * three cell types, each over-expressing its markers from src.markers
  * DE_GENE is ~4x higher in "Lapa" than "Dz" in every cell type
  * MT- genes make up a known fraction of counts
"""

import os
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

matplotlib.use("Agg")
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

CELL_TYPES = {
    "ISCs": ["LGR5", "AXIN2", "ASCL2", "OLFM4"],
    "Enterocytes": ["KRT20", "FABP2", "FABP1", "ALPI"],
    "EECs": ["CHGA", "CHGB", "NEUROD1", "PCSK1"],
}
DE_GENE = "IFI27"
DONORS = ["donor1", "donor2", "donor3"]
CONDITIONS = ["Dz", "Lapa"]


def make_adata(n_cells_per_group=40, n_background=500, seed=0):
    # n_background >= 500: sc.pp.calculate_qc_metrics reports top-500-gene proportions.
    import anndata as ad
    import scanpy as sc
    from src.markers import cell_cycle_markers

    rng = np.random.default_rng(seed)
    markers = [g for genes in CELL_TYPES.values() for g in genes]
    cc = cell_cycle_markers["S_genes"][:15] + cell_cycle_markers["G2M_genes"][:15]
    mt = ["MT-CO1", "MT-ND1", "MT-ATP6"]
    ribo = ["RPS3", "RPL7"]
    # Every gene the project's marker lists name, so marker-plotting code has data.
    import src.markers as M
    named = [g for d in (M.cell_type_markers, M.EEC_markers, M.EEC_marker_panels)
             for v in d.values() for g in ([v] if isinstance(v, str) else v)]
    named += M.Inflammed_Secretory_PC_markers + ["MKI67"]
    genes = list(dict.fromkeys(markers + cc + mt + ribo + [DE_GENE] + named
                               + [f"GENE{i}" for i in range(n_background)]))
    gidx = {g: i for i, g in enumerate(genes)}

    # Background genes get a log-normal spread of means, as real genes do. With a
    # flat background the markers would sit alone in the top expression bins and
    # Seurat-flavour HVG selection (dispersion z-scored within bins) would drop them.
    base = rng.lognormal(mean=0.0, sigma=1.5, size=len(genes))

    rows, obs = [], []
    for ct, ct_markers in CELL_TYPES.items():
        for donor in DONORS:
            for cond in CONDITIONS:
                mu = base.copy()
                mu[[gidx[g] for g in mt]] = 3.0
                mu[[gidx[g] for g in ct_markers]] = 25.0
                mu[gidx[DE_GENE]] = 20.0 if cond == "Lapa" else 5.0
                counts = rng.poisson(mu, size=(n_cells_per_group, len(genes)))
                rows.append(counts)
                obs += [{"participant": donor, "Condition": cond, "true_type": ct}] * n_cells_per_group

    X = np.vstack(rows).astype(np.float32)
    obs = pd.DataFrame(obs, index=[f"cell{i}" for i in range(len(X))]).astype("category")
    adata = ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=pd.DataFrame(index=genes))
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


@pytest.fixture
def adata():
    return make_adata()


@pytest.fixture
def labelled(adata):
    """Synthetic data with a manual_label column, as the 5_Labelling notebooks write."""
    adata.obs["manual_label"] = adata.obs["true_type"].astype("category")
    return adata
