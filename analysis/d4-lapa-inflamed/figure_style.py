"""
figure_style — matplotlib/scanpy figures pre-shaped for the lab's Figma slide template.

The slide deck has two standardized image slots in Figma, with two fixed aspect
ratios. This module produces matplotlib figures sized exactly to those slots, with
transparent backgrounds and projection-readable fonts, so a saved PNG drops into
the slot without rescaling (which would otherwise distort UMAPs / blur text).

Two slots (inches), see ``SLOT``:
  - "single" (10.32 x 7.78, ~4:3)   -> 2064 x 1556 px at dpi=200  — one figure
  - "wide"   (15.00 x 7.00, ~2.14:1) -> 3000 x 1400 px at dpi=200  — wide composites

These are double the on-slide pixel density, so they stay crisp at both 1080p and 4K.

Typical use:
    import figure_style as fs
    fig, ax = fs.make_figure(slot="single")
    ...                                  # draw onto ax
    fs.save_for_slide(fig, "my_panel.png")

UMAP / dotplot / time-course starters are in the ``make_*_figure`` helpers, and
``configure_scanpy(slot=...)`` makes scanpy's own plotters (sc.pl.umap, etc.)
inherit the slot sizing.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# --------------------------------------------------------------------------- #
# 1. Slot dimensions (inches). At dpi=200 -> 2064x1556 and 3000x1400 px.
# --------------------------------------------------------------------------- #
SLOT = {
    "single": (10.32, 7.78),  # 4:3-ish — for one figure (UMAP, ROC, heatmap)
    "wide":   (15.00, 7.00),  # 2.14:1 — for wide composites (two UMAPs side by side)
}

_STYLE_APPLIED = False


# --------------------------------------------------------------------------- #
# Font discovery — prefer Inter, fall back silently to Arial.
# --------------------------------------------------------------------------- #
def _inter_available() -> bool:
    """True if an 'Inter' font is usable, registering it with matplotlib if needed."""
    # Already known to matplotlib's font cache?
    if any(f.name == "Inter" for f in fm.fontManager.ttflist):
        return True
    # Otherwise ask the OS (fc-list on linux/mac) and register the file if found.
    if shutil.which("fc-list"):
        try:
            out = subprocess.run(
                ["fc-list", ":", "file", "family"],
                capture_output=True, text=True, timeout=5,
            ).stdout
        except (subprocess.SubprocessError, OSError):
            return False
        for line in out.splitlines():
            path, _, families = line.partition(":")
            if "Inter" in families and path.strip().lower().endswith((".ttf", ".otf")):
                try:
                    fm.fontManager.addfont(path.strip())
                except (RuntimeError, OSError):
                    continue
                if any(f.name == "Inter" for f in fm.fontManager.ttflist):
                    return True
    return False


# --------------------------------------------------------------------------- #
# 2. setup_style — projection-readable academic-slide rcParams (idempotent).
# --------------------------------------------------------------------------- #
def setup_style() -> None:
    """Apply slide-friendly matplotlib rcParams. Safe to call repeatedly."""
    global _STYLE_APPLIED
    if _STYLE_APPLIED:
        return

    sans = ["Inter", "Arial", "DejaVu Sans"]
    if not _inter_available():
        # Inter not present — drop it so matplotlib doesn't emit findfont warnings.
        sans = [f for f in sans if f != "Inter"]

    matplotlib.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": sans,
        "font.size": 14,
        "axes.titlesize": 16,
        "axes.titleweight": "medium",   # ~500
        "axes.labelsize": 14,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 1.0,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "legend.fontsize": 12,
        "legend.frameon": False,
        "lines.linewidth": 1.5,
        "savefig.transparent": True,
        "savefig.dpi": 200,
    })
    _STYLE_APPLIED = True


# --------------------------------------------------------------------------- #
# 3. make_figure — figure sized to a slot, transparent patch + axes.
# --------------------------------------------------------------------------- #
def make_figure(slot: str = "single", dpi: int = 200, nrows: int = 1, ncols: int = 1):
    """Return (fig, ax) — or (fig, axes) when nrows*ncols > 1 — sized to ``slot``.

    The figure patch and every axis facecolor are transparent so the PNG carries
    no white background into the slide.
    """
    setup_style()
    if slot not in SLOT:
        raise ValueError(f"unknown slot {slot!r}; choose from {sorted(SLOT)}")

    fig, axes = plt.subplots(nrows, ncols, figsize=SLOT[slot], dpi=dpi)
    fig.patch.set_alpha(0)
    for ax in fig.get_axes():
        ax.set_facecolor("none")
    return fig, axes


# --------------------------------------------------------------------------- #
# 4. save_for_slide — PNG, transparent, NO bbox_inches="tight".
# --------------------------------------------------------------------------- #
def save_for_slide(fig, path) -> str:
    """Save ``fig`` as a transparent PNG at the exact slot size. Returns abs path.

    Deliberately does NOT use bbox_inches="tight": trimming whitespace would change
    the figure's aspect ratio away from the slot target.
    """
    path = Path(path)
    if path.suffix.lower() != ".png":
        path = path.with_suffix(".png")
    fig.savefig(path, transparent=True)  # no bbox_inches -> keeps slot aspect ratio
    return str(path.resolve())


# --------------------------------------------------------------------------- #
# 5. configure_scanpy — make sc.pl.* inherit slot sizing (scanpy optional).
# --------------------------------------------------------------------------- #
def configure_scanpy(slot: str = "single") -> None:
    """Point scanpy's figure params at the given slot. No-op if scanpy is absent."""
    setup_style()
    if slot not in SLOT:
        raise ValueError(f"unknown slot {slot!r}; choose from {sorted(SLOT)}")
    try:
        import scanpy as sc
    except ImportError:
        return
    sc.set_figure_params(
        figsize=SLOT[slot], dpi=200, dpi_save=200,
        transparent=True, frameon=False,
    )


# --------------------------------------------------------------------------- #
# 6. Convenience wrappers for common single-cell figure types.
# --------------------------------------------------------------------------- #
def make_umap_figure(slot: str = "single"):
    """UMAP starter: no spines, no ticks, no axis labels, equal aspect."""
    fig, ax = make_figure(slot=slot)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_aspect("equal")
    return fig, ax


def make_dotplot_figure(slot: str = "single"):
    """Dotplot starter: both axes visible (genes x groups)."""
    fig, ax = make_figure(slot=slot)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return fig, ax


def make_time_course_figure(slot: str = "wide"):
    """Time-course starter: wide slot, x-axis on, y cleanly labelled."""
    fig, ax = make_figure(slot=slot)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlabel("Time")
    ax.set_ylabel("Value")
    return fig, ax


# --------------------------------------------------------------------------- #
# 8. Usage example — generate two test PNGs and verify their pixel dimensions.
# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    import numpy as np
    from PIL import Image
    from sklearn.datasets import make_blobs

    here = Path(__file__).resolve().parent

    def _check(path: str, expected: tuple[int, int]) -> None:
        with Image.open(path) as im:
            got = im.size  # (width, height)
        ok = "OK" if got == expected else "MISMATCH"
        print(f"  {path}\n    {got[0]}x{got[1]} px (expected {expected[0]}x{expected[1]}) -> {ok}")
        assert got == expected, f"{path}: got {got}, expected {expected}"

    X, y = make_blobs(n_samples=1000, centers=4, n_features=2, random_state=0)

    # 1) single-slot synthetic UMAP scatter -> 2064x1556
    fig, ax = make_umap_figure(slot="single")
    ax.scatter(X[:, 0], X[:, 1], c=y, cmap="tab10", s=8, linewidths=0)
    ax.set_title("Synthetic UMAP (single slot)")
    p1 = save_for_slide(fig, str(here / "example_umap.png"))
    plt.close(fig)
    print("wrote", p1)
    _check(p1, (2064, 1556))

    # 2) wide-slot two-panel side-by-side scatter -> 3000x1400
    fig, axes = make_figure(slot="wide", nrows=1, ncols=2)
    for ax, (i, j, ttl) in zip(axes, [(0, 1, "view A: dims 1-2"), (1, 0, "view B: dims 2-1")]):
        ax.scatter(X[:, i], X[:, j], c=y, cmap="tab10", s=8, linewidths=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_aspect("equal")
        ax.set_title(ttl)
    p2 = save_for_slide(fig, str(here / "example_wide.png"))
    plt.close(fig)
    print("wrote", p2)
    _check(p2, (3000, 1400))

    print("\nboth example figures produced at the expected dimensions.")
