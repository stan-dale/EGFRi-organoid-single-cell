"""
GSEA pre-rank enrichment analysis using gseapy.

Usage from a notebook:
    from src.gsea import GSEAEnrichmentAnalyzer
    from src.markers import pathway_gene_sets

    analyzer = GSEAEnrichmentAnalyzer(dge_file, output_dir, name,
                                      gene_sets=pathway_gene_sets)
    analyzer.run_full_analysis()
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gseapy as gp


class GSEAEnrichmentAnalyzer:
    """Run proper GSEA pre-rank enrichment using gseapy."""

    def __init__(self, dge_file, output_dir, name, gene_sets):
        """
        Parameters
        ----------
        dge_file : Path
            CSV of DESeq2 results (gene index, log2FoldChange, padj cols).
        output_dir : Path
            Where to write results and figures.
        name : str
            Human-readable label used in titles and filenames.
        gene_sets : dict[str, list[str]]
            Mapping of pathway name -> gene list.
        """
        self.dge_file = Path(dge_file)
        self.output_dir = Path(output_dir)
        self.name = name
        self.gene_sets = gene_sets
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.dge_df = None
        self.ranked_list = None
        self.gsea_results = {}

    def load_dge(self):
        """Load and preprocess DESeq2 results."""
        print(f"Loading DGE results from {self.dge_file.name}...")
        self.dge_df = pd.read_csv(self.dge_file, index_col=0)
        print(f"  Loaded {len(self.dge_df)} genes")

        required_cols = ["log2FoldChange", "padj"]
        before = len(self.dge_df)
        self.dge_df = self.dge_df.dropna(subset=required_cols)
        after = len(self.dge_df)
        print(f"  After filtering NA: {after} genes ({before - after} removed)")

    def create_ranked_list(self):
        """Create ranked gene list for GSEA prerank."""
        log10p = -np.log10(self.dge_df["padj"] + 1e-300)
        sign_fc = np.sign(self.dge_df["log2FoldChange"])
        self.dge_df["rank_metric"] = (
            log10p * sign_fc + 0.01 * self.dge_df["log2FoldChange"]
        )
        self.dge_df = self.dge_df.sort_values("rank_metric", ascending=False)
        self.ranked_list = self.dge_df["rank_metric"]

        print(f"\nRanked gene list: {len(self.ranked_list)} genes")
        print(f"  Top upregulated:   {self.ranked_list.index[0]}")
        print(f"  Top downregulated: {self.ranked_list.index[-1]}")

    def run_gsea(self):
        """Run GSEA prerank using gseapy."""
        print("\nRunning GSEA prerank...")
        pre_res = gp.prerank(
            rnk=self.ranked_list,
            gene_sets=self.gene_sets,
            min_size=5,
            max_size=500,
            permutation_num=1000,
            outdir=None,
            seed=42,
            verbose=False,
        )
        self.gsea_results = pre_res

        print("\nResults:")
        for _, row in pre_res.res2d.iterrows():
            print(f"  {row['Term']}:")
            print(f"    NES: {row['NES']:.3f}  "
                  f"FDR: {row['FDR q-val']:.4f}")

    def save_results(self):
        """Save GSEA results tables."""
        print(f"\nSaving results to {self.output_dir}...")
        results_table = self.gsea_results.res2d

        results_file = self.output_dir / f"{self.name}_gsea_results.csv"
        results_table.to_csv(results_file)
        print(f"  {results_file}")

        for _, row in results_table.iterrows():
            term = row["Term"]
            gene_set = self.gene_sets[term]
            genes_in_data = [g for g in gene_set if g in self.dge_df.index]
            if genes_in_data:
                gene_details = self.dge_df.loc[genes_in_data].copy()
                gene_details = gene_details.sort_values(
                    "log2FoldChange", ascending=False
                )
                detail_file = self.output_dir / f"{self.name}_{term}_genes.csv"
                gene_details.to_csv(detail_file)
                print(f"  {detail_file}")

    def plot_results(self):
        """Generate GSEA summary bar plot and per-pathway enrichment plots."""
        print("\nGenerating plots...")
        figures_dir = self.output_dir / "figures"
        figures_dir.mkdir(exist_ok=True)

        results_table = self.gsea_results.res2d

        # ── Summary bar plot ──
        fig, ax = plt.subplots(figsize=(8, 5))
        terms = results_table["Term"].tolist()
        nes_values = results_table["NES"].values
        colors = ["red" if x > 0 else "blue" for x in nes_values]

        ax.barh(terms, nes_values, color=colors, alpha=0.7)
        ax.axvline(0, color="black", linestyle="--", linewidth=1)
        ax.set_xlabel("Normalized Enrichment Score (NES)")
        ax.set_title(f"GSEA Results: {self.name}")
        plt.tight_layout()

        plot_file = figures_dir / f"{self.name}_gsea_barplot.pdf"
        plt.savefig(plot_file, bbox_inches="tight", dpi=150)
        plt.close()
        print(f"  {plot_file}")

        # ── Per-pathway enrichment plots ──
        for _, row in results_table.iterrows():
            term = row["Term"]
            try:
                enrichment_data = self.gsea_results.results[term]

                fig = plt.figure(figsize=(6, 5))
                gs = fig.add_gridspec(
                    3, 1, height_ratios=[1.5, 0.5, 0.5], hspace=0.3
                )

                # Running enrichment score
                ax1 = fig.add_subplot(gs[0])
                RES = enrichment_data["RES"]
                x = np.arange(len(RES))
                color = "green" if row["NES"] > 0 else "red"
                ax1.plot(x, RES, linewidth=2, color=color)
                ax1.axhline(y=0, color="gray", linestyle="--", linewidth=1)

                max_idx = np.argmax(np.abs(RES))
                ax1.plot(max_idx, RES[max_idx], "o", color="blue", markersize=5)

                sig_text = f"NES = {row['NES']:.3f}\nFDR = {row['FDR q-val']:.3f}"
                ax1.text(
                    0.02, 0.98, sig_text,
                    transform=ax1.transAxes, fontsize=9,
                    verticalalignment="top",
                    bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
                )
                ax1.set_ylabel("Enrichment Score", fontweight="bold")
                ax1.set_title(term, fontsize=11, fontweight="bold")
                ax1.set_xlim(0, len(RES))
                ax1.spines["top"].set_visible(False)
                ax1.spines["right"].set_visible(False)
                ax1.tick_params(labelbottom=False)

                # Hit barcode
                ax2 = fig.add_subplot(gs[1], sharex=ax1)
                hit_indices = enrichment_data["hits"]
                ax2.vlines(hit_indices, 0, 1, linewidth=0.5, color="black")
                ax2.set_ylim(0, 1)
                ax2.set_yticks([])
                ax2.set_ylabel("Genes", fontsize=8)
                ax2.spines["top"].set_visible(False)
                ax2.spines["right"].set_visible(False)
                ax2.spines["bottom"].set_visible(False)
                ax2.tick_params(labelbottom=False)
                ax2.set_facecolor("#f5f5f5")

                # Ranking metric
                ax3 = fig.add_subplot(gs[2], sharex=ax1)
                ranking = self.gsea_results.ranking
                ax3.fill_between(x, 0, ranking, color="gray", alpha=0.3)
                ax3.set_xlabel("Rank in Gene List", fontweight="bold")
                ax3.set_ylabel("Rank\nMetric", fontsize=8, rotation=0, labelpad=20)
                ax3.spines["top"].set_visible(False)
                ax3.spines["right"].set_visible(False)
                ax3.set_xlim(0, len(RES))

                plt.tight_layout()
                plot_file = figures_dir / f"{self.name}_{term}_gsea_plot.pdf"
                plt.savefig(plot_file, bbox_inches="tight", dpi=150)
                plt.close()
                print(f"  {plot_file}")

            except Exception as e:
                print(f"  WARNING: Could not create plot for {term}: {e}")

    def run_full_analysis(self):
        """Run the complete GSEA workflow."""
        print(f"\n{'=' * 60}")
        print(f"GSEA Analysis: {self.name}")
        print(f"{'=' * 60}")

        self.load_dge()
        self.create_ranked_list()
        self.run_gsea()
        self.save_results()
        self.plot_results()

        print(f"\nCompleted GSEA analysis for {self.name}\n")
