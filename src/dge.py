"""
Pseudobulk differential gene expression (DGE) pipeline using PyDESeq2.

Extracts the common workflow from all DGE notebooks into reusable functions.
Designed to be driven by a single parameterised notebook + YAML config.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
from pathlib import Path


# ── Pseudobulk creation ──

def create_pseudobulk(
    adata,
    cell_type,
    label_key="manual_label",
    condition_key="Condition",
    participant_key="participant",
    layer="counts",
    exclude_participants=None,
):
    """
    Subset *adata* to a cell type and aggregate counts per sample
    into a pseudobulk matrix.

    Parameters
    ----------
    adata : AnnData
    cell_type : str
        Value in adata.obs[label_key] to subset to.
    label_key : str
    condition_key : str
    participant_key : str
    layer : str
    exclude_participants : list[str], optional
        Participant labels to remove (e.g. ["Doublet", "Negative"]).

    Returns
    -------
    pb_counts : DataFrame
        Genes (rows) x samples (columns).
    """
    if exclude_participants is None:
        exclude_participants = ["Doublet", "Negative"]

    # Remove excluded participants
    mask = ~adata.obs[participant_key].isin(exclude_participants)
    ad = adata[mask].copy()

    # Subset to cell type
    sub = ad[ad.obs[label_key] == cell_type].copy()
    if sub.n_obs == 0:
        raise ValueError(f"No cells found for {label_key}=={cell_type!r}")

    sub.obs["pb_id"] = (
        sub.obs[participant_key].astype(str)
        + "_"
        + sub.obs[condition_key].astype(str)
    )

    # Extract counts
    counts = sub.layers[layer]
    if issparse(counts):
        counts = counts.tocsr()

    counts_df = pd.DataFrame.sparse.from_spmatrix(
        counts, index=sub.obs_names, columns=sub.var_names
    )

    # Aggregate per sample
    pb = counts_df.groupby(sub.obs["pb_id"]).sum()
    pb = pb.T  # genes x samples

    print(f"  Pseudobulk for {cell_type}: {pb.shape[1]} samples, {pb.shape[0]} genes, "
          f"{sub.n_obs} cells")
    return pb


def merge_pseudobulk(*dfs):
    """
    Merge pseudobulk DataFrames from multiple conditions/datasets.

    Fills missing genes with 0 and converts to int.
    """
    combined = pd.concat(dfs, axis=1)
    combined = combined.fillna(0).astype(int)
    return combined


def build_clinical_df(counts_df, conditions_map):
    """
    Build a clinical metadata DataFrame from sample IDs.

    Parameters
    ----------
    counts_df : DataFrame
        Rows = samples (after transposing pseudobulk).
    conditions_map : dict
        Maps condition substring -> condition label.
        E.g. {"Lapa": "Lapa", "Dz": "Dz", "AS": "AS"}

    Returns
    -------
    clinical_df : DataFrame
    """
    samples = counts_df.index.tolist()
    conditions = []
    for s in samples:
        matched = False
        for substr, label in conditions_map.items():
            if substr in s:
                conditions.append(label)
                matched = True
                break
        if not matched:
            conditions.append("Unknown")

    participants = [s.split("_")[0] for s in samples]

    return pd.DataFrame(
        {"condition": conditions, "participant": participants},
        index=samples,
    )


# ── DESeq2 ──

def run_pydeseq2(counts_df, clinical_df, contrast, design="condition"):
    """
    Run PyDESeq2 on pseudobulk counts.

    Parameters
    ----------
    counts_df : DataFrame
        Rows = samples, columns = genes (transposed pseudobulk).
    clinical_df : DataFrame
        Must have columns matching *design*.
    contrast : list
        E.g. ["condition", "Lapa", "Dz"] (factor, tested, reference).
    design : str
        Column name in clinical_df for the design formula.

    Returns
    -------
    results_df : DataFrame
        DESeq2 results with log2FoldChange, padj, etc.
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    dds = DeseqDataSet(
        counts=counts_df,
        metadata=clinical_df,
        design_factors=design,
    )
    dds.deseq2()

    stats = DeseqStats(dds, contrast=contrast)
    stats.summary()

    return stats.results_df


# ── Volcanoser prep ──

def prepare_for_volcanoser(results_df):
    """
    Add neg-log10 columns to PyDESeq2 results for use with Volcanoser.

    Computes ``neg_log10_pvalue`` and ``neg_log10_padj`` from the raw
    *pvalue* and *padj* columns.  NaN inputs produce NaN outputs.

    Parameters
    ----------
    results_df : DataFrame
        PyDESeq2 results with 'pvalue' and 'padj' columns.

    Returns
    -------
    DataFrame
        Copy of *results_df* with two extra columns appended.
    """
    df = results_df.copy()
    df["neg_log10_pvalue"] = -np.log10(df["pvalue"])
    df["neg_log10_padj"] = -np.log10(df["padj"])
    return df
