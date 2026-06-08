"""
Is the inflamed state different between D4 control and D4 lapatinib?

Rigorous, matched comparison:
  - Pool D4_DZ (control) + D4_Lapa (lapatinib), drop Doublet/Negative.
  - Normalise + score INFLAMMATION ONCE on the pooled object (shared background ->
    scores are comparable across conditions; AUCell & mean-z are rank/standardised
    so robust to this anyway).
  - Define the secretory compartment IDENTICALLY in both (negative selection).
  - Compare control vs lapatinib WITHIN the secretory compartment:
      Mann-Whitney U, Cliff's delta effect size, median shift, %REG+ cells.
  - Stratify by donor (H896/H897/H439 present in both) so a single line can't drive it.

Outputs: outputs/07_condition_comparison_stats.csv
         outputs/07_per_donor_medians.csv
         figures/07_control_vs_lapa_inflammation.png
"""
import scanpy as sc
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
FIG = HERE / "figures"; OUT = HERE / "outputs"
BASE = "/Users/standale/Library/CloudStorage/Dropbox-BCH/Stanley Dale/David and Stanley/EGFRi/data/data-objects/3.final-h5ad"

# ---- load + pool ----
parts = []
for f, cond in [("d4_dz", "control"), ("d4_lapa", "lapatinib")]:
    a = sc.read_h5ad(f"{BASE}/{f}_manual_labels.h5ad")
    a = a[a.obs["participant"].isin(["H896", "H897", "H439"])].copy()  # singlets only
    a.obs["cond"] = cond
    a.obs["donor"] = a.obs["participant"].astype(str)
    a.X = a.layers["counts"].copy()
    a.obs = a.obs[["cond", "donor", "manual_label"]].copy()
    a.var = a.var[[]].copy()
    parts.append(a)

ad = sc.concat(parts, join="inner", index_unique="-")
print("pooled", ad.shape, "| cond:", dict(ad.obs["cond"].value_counts()))

ad.layers["counts"] = ad.X.copy()
sc.pp.normalize_total(ad, target_sum=1e4)
sc.pp.log1p(ad)
ad.layers["lognorm"] = ad.X.copy()

# ---- secretory compartment, identical negative-selection definition ----
programs = {
    "Stem":       ["LGR5", "OLFM4", "ASCL2", "AXIN2", "SMOC2", "RNF43", "SOX9", "EPHB2"],
    "Absorptive": ["FABP2", "FABP1", "APOA1", "APOA4", "ALPI", "KRT20", "RBP2", "APOB", "ANPEP", "SI"],
    "Prolif":     ["MKI67", "TOP2A", "CDK1", "PCNA", "UBE2C", "BIRC5", "CCNB1", "CENPF"],
    "Secretory":  ["NEUROG3", "NEUROD1", "INSM1", "SOX4", "SPDEF", "DLL1", "DLL3", "HES6",
                   "BCL2", "PROX1", "MUC2", "TFF3", "SPINK4", "CHGA", "CHGB", "FCGBP"],
}
for n, gs in programs.items():
    sc.tl.score_genes(ad, [g for g in gs if g in ad.var_names], score_name=f"prog_{n}", use_raw=False)
ZMIN = 0.25
Z = pd.DataFrame({p: ad.obs[f"prog_{p}"] for p in programs}, index=ad.obs_names)
Z = (Z - Z.mean()) / Z.std()
argmax, amax = Z.idxmax(1), Z.max(1)
cls = pd.Series("Secretory", index=ad.obs_names, dtype=object)
cls[(argmax == "Absorptive") & (amax > ZMIN)] = "Absorptive"
cls[(argmax == "Stem") & (amax > ZMIN)] = "Stem"
ad.obs["lineage_class"] = cls
print("\nsecretory compartment by condition:")
print(pd.crosstab(ad.obs["cond"], ad.obs["lineage_class"]))

# ---- inflammation signatures ----
raw_sets = {
    "Inflamed_canonical": ["REG3A", "REG1B", "REG1A"],
    "ISP_core_LND":       ["LCN2", "NOS2", "DUOX2", "DUOXA2"],
    "ISP_Crohn":          ["LCN2", "DUOX2", "DUOXA2", "NOS2", "HLA-DRA", "HLA-DPA1",
                           "CD74", "CIITA", "CASP1", "CXCL1", "DMBT1"],
    "Chemokine_NFkB":     ["CCL20", "CXCL8", "CXCL1", "CXCL2", "CXCL3", "NFKBIA",
                           "IER3", "GADD45B", "JUN", "PPP1R15A"],
    "IFN_response":       ["STAT1", "STAT2", "IRF1", "IRF7", "ISG15", "IFIT1", "IFIT3",
                           "IFITM1", "IFITM3", "MX1", "OAS1", "B2M", "TAP1", "PSMB8", "PSMB9"],
}
det = np.asarray((ad.layers["counts"] > 0).sum(0)).ravel() > 0
det_genes = set(ad.var_names[det])
sets = {n: [g for g in gs if g in det_genes] for n, gs in raw_sets.items()}

for n, gs in sets.items():
    sc.tl.score_genes(ad, gs, score_name=f"sg_{n}", use_raw=False)

def gmat(genes):
    idx = [ad.var_names.get_loc(g) for g in genes]
    M = ad.layers["lognorm"][:, idx]
    return M.toarray() if sp.issparse(M) else np.asarray(M)
for n, gs in sets.items():
    M = gmat(gs); Zg = (M - M.mean(0)) / (M.std(0) + 1e-9)
    ad.obs[f"mz_{n}"] = Zg.mean(1)

def aucell(adata, gene_sets, max_rank_frac=0.05, chunk=3000):
    X = adata.layers["lognorm"]
    g = np.where(np.asarray((adata.layers["counts"] > 0).sum(0)).ravel() > 0)[0]
    pos = {gg: i for i, gg in enumerate(adata.var_names[g])}
    mr = int(np.ceil(max_rank_frac * len(g)))
    cols = {n: np.array([pos[x] for x in gs if x in pos]) for n, gs in gene_sets.items()}
    out = {n: np.zeros(adata.n_obs, np.float32) for n in gene_sets}
    Xd = X[:, g]
    for s in range(0, adata.n_obs, chunk):
        e = min(s + chunk, adata.n_obs)
        blk = Xd[s:e]; blk = blk.toarray() if sp.issparse(blk) else np.asarray(blk)
        order = np.argsort(-blk, axis=1, kind="stable")
        rk = np.empty_like(order); rk[np.arange(blk.shape[0])[:, None], order] = np.arange(blk.shape[1])[None, :]
        for n, c in cols.items():
            if len(c):
                out[n][s:e] = (mr - np.minimum(rk[:, c], mr)).sum(1) / (mr * len(c))
    return out
for n, v in aucell(ad, sets).items():
    ad.obs[f"auc_{n}"] = v

# REG+ flag (their original gating)
reg = [g for g in ["REG3A", "REG1B", "REG1A"] if g in ad.var_names]
ad.obs["REG_pos"] = (gmat(reg) > 0).any(1)

# ---- compare WITHIN the secretory compartment ----
sec = ad.obs[ad.obs["lineage_class"] == "Secretory"].copy()
def cliffs_delta(x, y):
    U, p = mannwhitneyu(x, y, alternative="two-sided")
    return 2 * U / (len(x) * len(y)) - 1, p

rows = []
for n in sets:
    for pre, lab in [("auc", "AUCell"), ("mz", "mean_z"), ("sg", "score_genes")]:
        col = f"{pre}_{n}"
        x = sec.loc[sec["cond"] == "lapatinib", col].values
        y = sec.loc[sec["cond"] == "control", col].values
        d, p = cliffs_delta(x, y)
        rows.append({"signature": n, "method": lab,
                     "median_lapa": np.median(x), "median_ctrl": np.median(y),
                     "cliffs_delta_lapa_vs_ctrl": d, "mwu_p": p})
stats = pd.DataFrame(rows)
stats.to_csv(OUT / "07_condition_comparison_stats.csv", index=False)
print("\n=== secretory compartment: lapatinib vs control (Cliff's delta, +=higher in lapa) ===")
print(stats[stats.method == "AUCell"][["signature", "median_ctrl", "median_lapa",
      "cliffs_delta_lapa_vs_ctrl", "mwu_p"]].round(4).to_string(index=False))

# REG+ fraction within secretory, by condition + donor
regfrac = sec.groupby(["cond", "donor"])["REG_pos"].mean().unstack("cond")
print("\n=== %REG+ within secretory compartment, by donor ===")
print((regfrac * 100).round(1).to_string())

# per-donor medians (AUCell) to show consistency of direction
pdm = sec.groupby(["donor", "cond"])[[f"auc_{n}" for n in sets]].median()
pdm.to_csv(OUT / "07_per_donor_medians.csv")
print("\n=== per-donor median AUCell (ISP_Crohn / Chemokine_NFkB) ===")
print(pdm[["auc_ISP_Crohn", "auc_Chemokine_NFkB"]].round(4).to_string())

# ---- figure: box/strip of AUCell by condition per signature, dots = donors ----
fig, axes = plt.subplots(1, len(sets), figsize=(4.2 * len(sets), 5.2))
donor_col = {"H896": "#1b9e77", "H897": "#d95f02", "H439": "#7570b3"}
for ax, n in zip(axes, sets):
    data = [sec.loc[sec["cond"] == c, f"auc_{n}"].values for c in ["control", "lapatinib"]]
    bp = ax.boxplot(data, positions=[0, 1], widths=0.5, showfliers=False, patch_artist=True)
    for patch, col in zip(bp["boxes"], ["#bdbdbd", "#fc8d62"]):
        patch.set_facecolor(col)
    # donor median dots
    for ci, c in enumerate(["control", "lapatinib"]):
        for dn, col in donor_col.items():
            m = sec.loc[(sec["cond"] == c) & (sec["donor"] == dn), f"auc_{n}"].median()
            ax.scatter(ci, m, color=col, s=55, edgecolor="k", zorder=3, label=dn if (ci == 0 and n == list(sets)[0]) else None)
    d = stats[(stats.signature == n) & (stats.method == "AUCell")]["cliffs_delta_lapa_vs_ctrl"].values[0]
    ax.set_title(f"{n}\nCliff's d={d:+.2f}", fontsize=10)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["control", "lapa"])
    ax.set_ylabel("AUCell" if n == list(sets)[0] else "")
axes[0].legend(title="donor median", fontsize=8, loc="upper left")
fig.suptitle("Inflammation within the secretory compartment: D4 control vs lapatinib (matched pipeline)", y=1.02)
fig.tight_layout(); fig.savefig(FIG / "07_control_vs_lapa_inflammation.png", dpi=140, bbox_inches="tight")
plt.close(fig)
print("\nwrote", FIG / "07_control_vs_lapa_inflammation.png")
print("DONE")
