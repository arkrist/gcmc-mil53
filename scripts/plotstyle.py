"""Shared matplotlib style for the project figures.

Two fixed series colours (never cycled): slot 1 = our simulation, slot 2 =
reference/literature. Markers differ too, so the series are distinguishable
without colour (colour-blind readers, greyscale print).
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SIM = dict(color="#2a78d6", marker="o", mfc="#2a78d6", ms=7, lw=2, capsize=3)
REF = dict(color="#eb6834", marker="s", mfc="white", ms=8, lw=0, elinewidth=1.5, capsize=3)
EXTRA = ["#1baf7a", "#eda100", "#e87ba4"]  # further reference sets, in this order

plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 200, "savefig.bbox": "tight",
    "font.size": 11, "axes.titlesize": 12, "axes.labelsize": 11,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": "#e3e3e0", "grid.linewidth": 0.8,
    "axes.edgecolor": "#6b6b66", "xtick.color": "#3d3d3a", "ytick.color": "#3d3d3a",
    "legend.frameon": False,
})
