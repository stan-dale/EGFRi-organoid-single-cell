"""
UMAP figures for publication.

umap_panel() draws the labelled panel: layout follows Veres et al. (2019) Nature
569:368 — a UMAP scatter with a text legend to its left and a population-ratio
bar beneath, sized to half of an 89 mm single-column figure. Colours come from a
{label: hex} dict rather than a positional cycle, so a cell type keeps its colour
across every figure it appears in.

gene_umaps() draws per-gene expression panels in the project's turbo style, the
one used for the TEAD, YAP1 and tafazzin figures.
"""

from __future__ import annotations

import warnings
from collections import Counter
from contextlib import contextmanager
from pathlib import Path

import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgb
from sklearn.neighbors import NearestNeighbors

from .figure_style import SLOT as SLIDE_SLOT

MM_PER_INCH = 25.4
PANEL_W_MM = 89 / 2        # half a single-column (89 mm) figure
PANEL_H_MM = 25.4
RATIO_BAR_FRAC = 0.07      # share of panel height given to the ratio bar

# Output sizes. "print" is the native panel; the rest borrow their widths from
# figure_style's slide slots. The panel is 1.75:1 and no slot matches that, so a
# slot sets the width only and the panel keeps its own aspect — scaling is
# uniform, and fonts and dot sizes scale with it, so the design is unchanged.
PANEL_SIZE_IN = (PANEL_W_MM / MM_PER_INCH, PANEL_H_MM / MM_PER_INCH)
SLOT_WIDTH_IN = {"print": PANEL_SIZE_IN[0], **{k: v[0] for k, v in SLIDE_SLOT.items()}}


def slot_scale(slot: str) -> float:
    """Uniform scale that makes the panel as wide as `slot`."""
    if slot not in SLOT_WIDTH_IN:
        raise ValueError(f"unknown slot {slot!r}; choose from {sorted(SLOT_WIDTH_IN)}")
    return SLOT_WIDTH_IN[slot] / PANEL_SIZE_IN[0]

# Layered scatter sizes: black ring, white halo, coloured fill.
DOT_OUTLINE, DOT_HALO, DOT_FILL = 4, 2, 1.5

# Legend and ratio-bar geometry, in axes coordinates.
LINE_STEP = 1 / 10.5
X_DOT, X_TEXT = 0.22, 0.30
BAR_X0, BAR_W = 0.30, 0.68
MARGIN = 0.012          # figure-fraction inset so nothing overhangs the edge
TITLE_BAND = 0.11       # figure fraction reserved above the panel for the title

PAPER_RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Arial", "DejaVu Sans"],
    "axes.titlesize": 5,
    "axes.labelsize": 5,
    "xtick.labelsize": 5,
    "ytick.labelsize": 5,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "pdf.fonttype": 42,     # editable text in Illustrator
}


@contextmanager
def paper_style(scale: float = 1.0):
    """Apply the print rcParams for one block, then restore what was there.

    Point sizes scale linearly, so a scaled panel looks identical to the native
    one enlarged rather than the same text in a bigger frame.
    """
    rc = {k: (v * scale if isinstance(v, (int, float)) else v)
          for k, v in PAPER_RC.items()}
    rc["pdf.fonttype"] = PAPER_RC["pdf.fonttype"]
    with plt.rc_context(rc):
        yield


def homogeneity_order(X: np.ndarray, labels: np.ndarray,
                      n_neighbors: int = 10) -> np.ndarray:
    """Draw order that puts cells with agreeing neighbours on top.

    Scores each cell by the fraction of its nearest neighbours in UMAP space
    sharing its label, then sorts ascending, so mixed boundary cells are painted
    first and clean cluster cores land over them.
    """
    knn = NearestNeighbors(n_neighbors=n_neighbors).fit(X).kneighbors(
        X, return_distance=False)
    # Column 0 is the cell itself, since we query the points we fitted on.
    same = (labels[knn] == labels[knn][:, [0]]).mean(axis=1)
    return np.argsort(same, kind="stable")


def resolve_colors(labels: np.ndarray, colors: dict) -> np.ndarray:
    """Map an (n,) array of label strings to an (n, 3) RGB array."""
    out = np.empty((len(labels), 3))
    for name, rgb in colors.items():
        out[labels == name] = rgb
    return out


def draw_scatter(ax, X: np.ndarray, colors: np.ndarray, alpha: float = 0.7,
                 scale: float = 1.0) -> None:
    """Paint the UMAP as a black ring, a white halo, then the coloured fill."""
    ax.set(xticks=[], yticks=[])
    ax.set_frame_on(False)
    ax.patch.set_visible(False)
    # Marker size is an area in points squared, so it scales by the square.
    area = scale ** 2
    ax.scatter(X[:, 0], X[:, 1], c="k", s=DOT_OUTLINE * area,
               edgecolor="none", rasterized=True)
    ax.scatter(X[:, 0], X[:, 1], c="w", s=DOT_HALO * area,
               edgecolor="none", rasterized=True)
    ax.scatter(X[:, 0], X[:, 1], c=colors, s=DOT_FILL * area, alpha=alpha,
               edgecolor="none", rasterized=True)


def draw_legend(ax, entries: list, colors: dict, n_cells: int,
                scale: float = 1.0) -> None:
    """Cell count, then a colour dot and name for each entry."""
    ax.set(xticks=[], yticks=[], xlim=(0, 1), ylim=(0, 1))
    ax.set_frame_on(False)

    line = 2
    ax.text(X_TEXT, 1 - line * LINE_STEP, f"({n_cells:,} cells)",
            va="center", fontsize=5 * scale, clip_on=False)
    line += 1.5

    for name in entries:
        y = 1 - line * LINE_STEP
        ax.scatter(X_DOT, y + 0.008, s=15 * scale ** 2, c=[colors[name]],
                   clip_on=False)
        ax.text(X_TEXT, y, name, va="center", fontsize=5 * scale, clip_on=False)
        line += 1


def draw_ratio_bar(ax, labels: np.ndarray, entries: list, colors: dict,
                   scale: float = 1.0) -> None:
    """Population fractions stacked in `entries` order, inside a bordered box."""
    ax.set(xticks=[], yticks=[], xlim=(0, 1), ylim=(0, 1))
    ax.set_frame_on(False)

    fracs = pd.Series(Counter(labels), dtype=float)
    fracs /= fracs.sum()

    unlisted = set(fracs.index) - set(entries)
    if unlisted:
        warnings.warn(f"left out of the ratio bar: {', '.join(sorted(unlisted))}")

    left = 0.0
    for name in entries:
        if name not in fracs:
            continue
        width = float(fracs[name])
        ax.add_patch(patches.Rectangle(
            (BAR_X0 + BAR_W * left, 0.02), BAR_W * width, 1.0,
            facecolor=colors[name], edgecolor="none", clip_on=False))
        left += width

    ax.text(0.175, 0.5, "Pop. ratios:", va="center", ha="center",
            fontsize=5 * scale, clip_on=False)
    ax.add_patch(patches.Rectangle(
        (BAR_X0, 0.02), BAR_W, 1.0, facecolor="none",
        edgecolor="k", linewidth=0.5 * scale, clip_on=False))


def umap_panel(adata, label_key: str, palette: dict, order: list | None = None,
               highlight: str | None = None, title: str | None = None,
               save=None, slot: str = "print", dpi: int = 600,
               n_neighbors: int = 10, muted: str = "#cccccc"):
    """Draw one Veres-style panel and return the figure.

    `palette` maps label to hex colour and must cover every label present.
    `order` fixes the legend and ratio-bar sequence; unlisted labels warn rather
    than silently vanishing. `highlight` greys every label but one, keeping the
    full ratio bar for context. `slot` sets the output width — "print" is the
    native 44.5 x 25.4 mm panel, "single" matches figure_style's slide slot.
    """
    scale = slot_scale(slot)
    labels = adata.obs[label_key].astype(str).to_numpy()
    X = adata.obsm["X_umap"]
    present = list(pd.unique(labels))

    entries = [n for n in (order if order is not None else present) if n in present]

    if highlight is None:
        unknown = sorted(set(present) - set(palette))
        if unknown:
            raise KeyError(
                f"no colour for {', '.join(unknown)} — add them to the palette")
        colors = {name: to_rgb(palette[name]) for name in present}
        legend_entries = entries
    else:
        if highlight not in present:
            raise ValueError(f"{highlight!r} is not a value in {label_key!r}")
        colors = {name: to_rgb(muted) for name in present}
        colors[highlight] = to_rgb(palette[highlight])
        legend_entries = [highlight]

    with paper_style(scale):
        w, h = PANEL_SIZE_IN[0] * scale, PANEL_SIZE_IN[1] * scale
        heights = ((1 - RATIO_BAR_FRAC) * h, RATIO_BAR_FRAC * h)  # row ratios
        widths = (w - heights[0], heights[0])

        fig = plt.figure(figsize=(w, h), dpi=dpi)
        # A small inset keeps the title and the ratio bar's border inside the
        # figure. Without it they overhang and get cut, since nothing trims.
        top = 1 - (TITLE_BAND if title else MARGIN)
        gs = gridspec.GridSpec(2, 2, fig, 0, MARGIN, 1, top,
                               hspace=0, wspace=0,
                               width_ratios=widths, height_ratios=heights)

        if title:
            # Its own band across the full width, so length never collides
            # with the scatter the way an in-legend title does.
            fig.text(X_TEXT * widths[0] / w, 1 - TITLE_BAND / 2, title,
                     va="center", fontsize=6 * scale)

        draw_legend(fig.add_subplot(gs[0, 0]), legend_entries, colors,
                    adata.n_obs, scale)

        row = homogeneity_order(X, labels, n_neighbors)
        draw_scatter(fig.add_subplot(gs[0, 1], zorder=-1), X[row],
                     resolve_colors(labels, colors)[row], scale=scale)

        draw_ratio_bar(fig.add_subplot(gs[1, :]), labels, entries, colors, scale)

        # Deliberately no bbox_inches="tight": trimming to content would give
        # every panel a different size, which is the opposite of a fixed slot.
        if save is not None:
            save = Path(save)
            save.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save, dpi=dpi, transparent=True)
    return fig


# Every gene panel is drawn at this width, so panels match across figures and the
# figure simply grows with the number of genes. SCANPY_PANEL_IN is the panel width
# scanpy's own 120000/n_cells dot size assumes.
GENE_PANEL_IN = 3.0
SCANPY_PANEL_IN = 4.0


def gene_umaps(adata, genes, ncols: int = 4, panel_in: float = GENE_PANEL_IN,
               size: float | None = None, cmap: str = "turbo", vmax="p100",
               use_raw: bool | None = None, save=None, show: bool = True):
    """Per-gene expression UMAPs in the project's turbo style.

    Follows the whole-atlas TEAD/YAP1 scripts: turbo, no frame. Genes absent
    from `adata` are dropped and reported.

    Two defaults differ from those scripts, because this runs on a small subset.
    `vmax` is "p100", not "p95" — a hormone expressed in 4% of cells has a 95th
    percentile of zero, so a p95 ceiling lands on background and the panel reads
    flat. `size` defaults to scanpy's 120000/n_cells, rescaled for the
    panel width, rather than the pinned 15 that suits ~33k cells in one panel. Pass `use_raw=True` to colour from
    log-normalised values when .X has been scaled.
    """
    import scanpy as sc

    genes = [genes] if isinstance(genes, str) else list(genes)
    lookup = adata.raw.var_names if (use_raw and adata.raw is not None) else adata.var_names
    present = [g for g in genes if g in lookup]
    missing = [g for g in genes if g not in lookup]
    if missing:
        print(f"not in var_names: {', '.join(missing)}")
    if not present:
        print("no genes to plot")
        return None

    if size is None:
        # Dot area is in points squared, so it scales with the panel area.
        size = 120000 / adata.n_obs * (panel_in / SCANPY_PANEL_IN) ** 2

    with plt.rc_context({"figure.figsize": (panel_in, panel_in)}):
        sc.pl.umap(adata, color=present, ncols=min(ncols, len(present)), cmap=cmap,
                   vmax=vmax, size=size, use_raw=use_raw, frameon=False, show=False)
        fig = plt.gcf()
        if save is not None:
            save = Path(save)
            save.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save, bbox_inches="tight")
    if show:
        # Render once. Returning the figure too would make the inline backend
        # draw it a second time.
        plt.show()
        plt.close(fig)
        return None
    return fig
