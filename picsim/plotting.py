"""Shared figure style."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams.update({
    "figure.dpi": 110, "font.size": 9, "axes.grid": True, "grid.alpha": 0.3,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
})


def fig(nrows=1, ncols=1, w=6.5, h=None):
    h = h or 2.6 * nrows
    return plt.subplots(nrows, ncols, figsize=(w, h), squeeze=False)
