"""Figure helpers produce files without needing real data."""

import matplotlib.pyplot as plt
import numpy as np
import pytest

from src.palette import celltype_palette


@pytest.fixture
def embedded(labelled):
    rng = np.random.default_rng(0)
    centres = {"ISCs": (0, 0), "Enterocytes": (5, 0), "EECs": (0, 5)}
    labelled.obsm["X_umap"] = np.array(
        [centres[t] for t in labelled.obs["true_type"]]) + rng.normal(size=(labelled.n_obs, 2))
    labelled.obs["leiden"] = labelled.obs["true_type"].cat.codes.astype(str).astype("category")
    labelled.obs["doublet_score"] = rng.uniform(size=labelled.n_obs)
    return labelled


def test_umap_panel_saves_pdf(embedded, tmp_path):
    from src.plotting import umap_panel

    out = tmp_path / "panel.pdf"
    fig = umap_panel(embedded, "manual_label", palette=celltype_palette,
                     order=["ISCs", "Enterocytes", "EECs"], title="t", save=out, dpi=72)
    assert out.exists() and out.stat().st_size > 0
    plt.close(fig)


def test_umap_panel_rejects_uncoloured_labels(embedded):
    from src.plotting import umap_panel

    with pytest.raises(KeyError, match="no colour"):
        umap_panel(embedded, "manual_label", palette={"ISCs": "#000000"})


def test_gene_umaps(embedded):
    from src.plotting import gene_umaps

    gene_umaps(embedded, ["LGR5", "CHGA", "NOT_A_GENE"], show=False)
    plt.close("all")


def test_save_doublet_umaps_writes_one_png_per_cluster(embedded, tmp_path):
    from src.vis import save_doublet_umaps

    save_doublet_umaps(embedded, "synthetic", base_dir=tmp_path, dpi=50)
    pngs = sorted(p.name for p in (tmp_path / "synthetic").glob("*.png"))
    assert pngs == [f"doublet_score_cluster_{c}.png" for c in ("0", "1", "2")]


def test_veres_panel_saves(embedded, tmp_path):
    from src.vis import plot_veres_panel

    out = tmp_path / "veres.png"
    plot_veres_panel(embedded, "manual_label", palette=celltype_palette, save=out, dpi=72)
    assert out.exists()
    plt.close("all")
