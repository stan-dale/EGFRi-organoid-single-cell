"""Pseudobulk DESeq2 and GSEA on synthetic data with a planted DE gene."""

import numpy as np
import pytest

from conftest import DE_GENE, DONORS
from src.dge import (build_clinical_df, create_pseudobulk, merge_pseudobulk,
                     prepare_for_volcanoser, run_pydeseq2)


def split_by_condition(adata):
    return [adata[adata.obs["Condition"] == c].copy() for c in ("Dz", "Lapa")]


def test_pseudobulk_sums_counts_per_donor(labelled):
    pb = create_pseudobulk(labelled, "ISCs")
    assert sorted(pb.columns) == sorted(f"{d}_{c}" for d in DONORS for c in ("Dz", "Lapa"))
    mask = ((labelled.obs["manual_label"] == "ISCs") & (labelled.obs["participant"] == "donor1")
            & (labelled.obs["Condition"] == "Lapa")).to_numpy()
    expected = np.asarray(labelled.layers["counts"][mask].sum(0)).ravel()
    np.testing.assert_array_equal(pb["donor1_Lapa"].to_numpy(), expected)


def test_pseudobulk_matches_label_aliases(labelled):
    """D2_Lapa calls them 'Secretory progenitor cells'; DGE config asks for 'Secretory PCs'."""
    labelled.obs["manual_label"] = labelled.obs["manual_label"].astype(str).replace(
        "EECs", "Secretory progenitor cells")
    pb = create_pseudobulk(labelled, "Secretory PCs")
    assert pb.shape[1] == 6


def test_pseudobulk_reports_available_labels(labelled):
    with pytest.raises(ValueError, match="Labels present"):
        create_pseudobulk(labelled, "Goblet cells")


def test_pseudobulk_excludes_unassigned_donors(labelled):
    labelled.obs["participant"] = labelled.obs["participant"].astype(str)
    labelled.obs.iloc[:5, labelled.obs.columns.get_loc("participant")] = "Doublet"
    pb = create_pseudobulk(labelled, "ISCs")
    assert not any(c.startswith("Doublet") for c in pb.columns)


def test_clinical_df_parses_sample_ids(labelled):
    counts = merge_pseudobulk(*[create_pseudobulk(a, "ISCs") for a in split_by_condition(labelled)]).T
    clin = build_clinical_df(counts, {"Lapa": "Lapa", "Dz": "Dz"})
    assert set(clin["condition"]) == {"Lapa", "Dz"}
    assert set(clin["participant"]) == set(DONORS)


def test_deseq2_finds_planted_gene(labelled):
    """End to end as the DGE notebook runs it: two datasets -> merge -> DESeq2."""
    pbs = [create_pseudobulk(a, "ISCs") for a in split_by_condition(labelled)]
    counts = merge_pseudobulk(*pbs).T
    clin = build_clinical_df(counts, {"Lapa": "Lapa", "Dz": "Dz"})
    res = run_pydeseq2(counts, clin, contrast=["condition", "Lapa", "Dz"])

    top = res.dropna(subset=["padj"]).sort_values("padj").index[0]
    assert top == DE_GENE
    assert res.loc[DE_GENE, "log2FoldChange"] == pytest.approx(2.0, abs=0.4)  # 20 vs 5

    v = prepare_for_volcanoser(res)
    assert v.loc[DE_GENE, "neg_log10_padj"] > 5


def test_gsea_ranks_planted_gene_set(labelled, tmp_path):
    from src.gsea import GSEAEnrichmentAnalyzer

    pbs = [create_pseudobulk(a, "ISCs") for a in split_by_condition(labelled)]
    counts = merge_pseudobulk(*pbs).T
    res = run_pydeseq2(counts, build_clinical_df(counts, {"Lapa": "Lapa", "Dz": "Dz"}),
                       contrast=["condition", "Lapa", "Dz"])
    dge_csv = tmp_path / "dge.csv"
    res.to_csv(dge_csv)

    genes = list(res.index)
    gene_sets = {"PLANTED": [DE_GENE] + genes[-6:], "BACKGROUND": genes[50:60]}
    a = GSEAEnrichmentAnalyzer(dge_csv, tmp_path / "gsea", "synthetic", gene_sets)
    a.run_full_analysis()
    out = a.gsea_results.res2d.set_index("Term")
    assert (tmp_path / "gsea" / "synthetic_gsea_results.csv").exists()
    assert (tmp_path / "gsea" / "figures" / "synthetic_gsea_barplot.pdf").exists()
    assert set(out.index) == {"PLANTED", "BACKGROUND"}
