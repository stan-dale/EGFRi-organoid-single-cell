"""Build a pipeline-compatible h5ad from a cellranger 10x output.

Reproduces the format that the R conversion (`r-src/2_convert.r`) produces, so
the same h5ad can drop into the pipeline starting at notebook 0_Quality_Control.
This is the entry point if you're starting from GEO data (cellranger output)
rather than the Dropbox Seurat objects.

Output structure (matches R-converted h5ad):
    adata.X                    log-normalized (log1p of normalize_total to 10_000)
    adata.layers["counts"]     raw integer counts
    adata.uns["X_name"]        "logcounts"
    adata.obs columns:
        orig.ident             sample id (e.g. "D10_Lapa")
        nCount_RNA             total counts per cell
        nFeature_RNA           number of expressed genes per cell
        freemuxlet.identity    donor demultiplexing result (NA if --freemuxlet not provided)
        participant            donor assignment (from freemuxlet)
        Condition              sample condition (e.g. "Lapa")
        Time_point             integer day (2, 4, 10)
        Treatment              treatment label (e.g. "Lapa", "AS")
        ident                  copy of orig.ident

Usage:
    python r-src/build_h5ad_from_10x.py \\
        --tenx-h5 /path/to/D10_Lapa/outs/filtered_feature_bc_matrix.h5 \\
        --sample D10_Lapa \\
        --condition Lapa \\
        --timepoint 10 \\
        --treatment Lapa \\
        --out data/data-objects/raw/egfDuod_D10_Lapa_DZ.h5ad

With donor demultiplexing (recommended — needed for per-donor analyses):
    python r-src/build_h5ad_from_10x.py \\
        --tenx-h5 D10_Lapa/outs/filtered_feature_bc_matrix.h5 \\
        --sample D10_Lapa \\
        --condition Lapa \\
        --timepoint 10 \\
        --treatment Lapa \\
        --freemuxlet D10_Lapa/freemuxlet_output.tsv \\
        --out data/data-objects/raw/egfDuod_D10_Lapa_DZ.h5ad

The --freemuxlet TSV is expected to have columns `barcode`, `freemuxlet.identity`,
and `participant` (matching the obs columns the pipeline expects). Adjust the
CSV/TSV reader below if your demuxlet output has a different schema.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


def build(args: argparse.Namespace) -> None:
    print(f"reading 10x: {args.tenx_h5}")
    adata = sc.read_10x_h5(args.tenx_h5)
    adata.var_names_make_unique()
    print(f"  loaded {adata.shape[0]:,} barcodes × {adata.shape[1]:,} genes")

    # Optional cell-filtering — the upstream Seurat pipeline filtered to fewer
    # cells than the raw 10x has (e.g. raw=2.6M barcodes → final=33,507 cells
    # for D10 Lapa). Use --min-counts / --min-features to approximate this
    # filtering. Exact reproduction of the published cell counts may require
    # additional filtering (mito %, doublets) applied by the original pipeline.
    if args.min_counts > 0:
        sc.pp.filter_cells(adata, min_counts=args.min_counts)
        print(f"  after min_counts={args.min_counts}: {adata.shape[0]:,} cells")
    if args.min_features > 0:
        sc.pp.filter_cells(adata, min_genes=args.min_features)
        print(f"  after min_features={args.min_features}: {adata.shape[0]:,} cells")
    if args.min_cells_per_gene > 0:
        sc.pp.filter_genes(adata, min_cells=args.min_cells_per_gene)
        print(f"  after min_cells_per_gene={args.min_cells_per_gene}: {adata.shape[1]:,} genes")
    # Drop scanpy-injected QC fields that the R-converted h5ad doesn't have.
    adata.obs.drop(columns=[c for c in ("n_genes", "n_counts") if c in adata.obs.columns], inplace=True)
    if "n_cells" in adata.var.columns:
        adata.var.drop(columns=["n_cells"], inplace=True)

    # Preserve raw counts in a layer; X will be log-normalized below.
    adata.layers["counts"] = adata.X.copy()

    # nCount_RNA / nFeature_RNA — match Seurat's per-cell QC fields.
    counts = adata.layers["counts"]
    n_count = np.asarray(counts.sum(axis=1)).ravel().astype(float)
    n_feature = np.asarray((counts > 0).sum(axis=1)).ravel().astype(int)

    adata.obs["orig.ident"] = args.sample
    adata.obs["nCount_RNA"] = n_count
    adata.obs["nFeature_RNA"] = n_feature

    # Donor demultiplexing — only present if a freemuxlet TSV was supplied.
    if args.freemuxlet:
        print(f"reading freemuxlet: {args.freemuxlet}")
        fm = pd.read_csv(args.freemuxlet, sep="\t")
        if "barcode" not in fm.columns:
            raise SystemExit(
                f"--freemuxlet TSV must have a 'barcode' column; got {list(fm.columns)}"
            )
        fm = fm.set_index("barcode")
        adata.obs = adata.obs.join(fm, how="left")
    else:
        # Use empty string rather than pd.NA — h5py can't write a column whose
        # dtype is purely-NA (it needs concrete strings or floats).
        adata.obs["freemuxlet.identity"] = ""
        adata.obs["participant"] = ""

    adata.obs["Condition"] = args.condition
    adata.obs["Time_point"] = args.timepoint
    adata.obs["Treatment"] = args.treatment
    adata.obs["ident"] = args.sample

    # Log-normalize (matches Seurat's NormalizeData(LogNormalize, scale.factor=10000)).
    sc.pp.normalize_total(adata, target_sum=10_000)
    sc.pp.log1p(adata)

    adata.uns["X_name"] = "logcounts"

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out)
    print(f"wrote {out}  ({adata.shape[0]:,} cells × {adata.shape[1]:,} genes)")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tenx-h5", required=True, type=Path,
                   help="Path to cellranger filtered_feature_bc_matrix.h5")
    p.add_argument("--sample", required=True,
                   help="Sample identifier (e.g. D10_Lapa) — written to orig.ident and ident")
    p.add_argument("--condition", required=True,
                   help="Condition label (e.g. Dz, Lapa, AS)")
    p.add_argument("--timepoint", required=True, type=int,
                   help="Day as integer (2, 4, or 10)")
    p.add_argument("--treatment", required=True,
                   help="Treatment label (e.g. Lapa, AS, Dz)")
    p.add_argument("--freemuxlet", type=Path, default=None,
                   help="Optional: TSV with columns barcode, freemuxlet.identity, participant")
    p.add_argument("--min-counts", type=int, default=0,
                   help="Drop barcodes with fewer than this many total UMIs (default 0 = keep all)")
    p.add_argument("--min-features", type=int, default=0,
                   help="Drop barcodes with fewer than this many expressed genes (default 0 = keep all)")
    p.add_argument("--min-cells-per-gene", type=int, default=3,
                   help="Drop genes expressed in fewer than this many cells "
                        "(default 3 — matches Seurat CreateSeuratObject's default)")
    p.add_argument("--out", required=True, type=Path,
                   help="Destination h5ad path (e.g. data/data-objects/raw/egfDuod_D10_Lapa_DZ.h5ad)")
    build(p.parse_args())


if __name__ == "__main__":
    main()
