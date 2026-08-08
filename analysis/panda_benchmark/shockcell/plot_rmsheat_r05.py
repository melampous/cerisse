#!/usr/bin/env python3
"""Figure 3.9-style RMS heatmaps at r/D_e = 0.5."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).parent
D = 0.0254
P_AMB = 99780.0
U_J = 414.2
XMAX = 8.0

# Retain the Figure 3.9 colour limits for direct visual comparison.
VMAX = {"p": 0.42, "u": 0.40}

plt.rcParams.update({
    "font.size": 9,
    "axes.linewidth": 0.6,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
})


def box_filter(values, x, width_over_D):
    width = max(3, int(round(width_over_D / (x[1] - x[0]))))
    if width % 2 == 0:
        width += 1
    result = np.convolve(values, np.ones(width) / width, mode="same")
    result[:width // 2] = values[:width // 2]
    result[-(width // 2):] = values[-(width // 2):]
    return result


def load_profile(mode, tag, quantity):
    data = np.load(HERE / f"r05_{mode}_L{tag}.npz")
    variance_key = (
        "pressure_variance" if quantity == "p" else "axial_velocity_variance"
    )
    scale = P_AMB if quantity == "p" else U_J
    rms = np.sqrt(np.maximum(data[variance_key], 0.0)) / scale
    return data["x"], rms, int(round(D / float(data["dx"])))


x_edges = np.linspace(0.0, XMAX, 801)
x_centres = 0.5 * (x_edges[:-1] + x_edges[1:])
heatmaps = {}

for mode in ("2d", "3d"):
    tags = "0123456" if mode == "2d" else "12345"
    for quantity in ("p", "u"):
        values = np.full((len(tags), len(x_centres)), np.nan)
        resolutions = []
        for row, tag in enumerate(tags):
            x, rms, resolution = load_profile(mode, tag, quantity)
            filtered = box_filter(rms, x, 0.05)
            values[row] = np.interp(x_centres, x, filtered)
            resolutions.append(resolution)
        heatmaps[(mode, quantity)] = values, resolutions


fig, axes = plt.subplots(2, 2, figsize=(7.2, 4.4), sharex=True)
for row, mode in enumerate(("2d", "3d")):
    for column, (quantity, colourmap) in enumerate((("p", "magma"), ("u", "viridis"))):
        axis = axes[row, column]
        values, resolutions = heatmaps[(mode, quantity)]
        image = axis.pcolormesh(
            x_edges,
            np.arange(len(resolutions) + 1),
            values,
            cmap=colourmap,
            vmin=0.0,
            vmax=VMAX[quantity],
            rasterized=True,
        )
        axis.set_yticks(np.arange(len(resolutions)) + 0.5)
        axis.set_yticklabels([str(value) for value in resolutions], fontsize=7.5)
        axis.tick_params(length=2.5)
        colourbar = fig.colorbar(image, ax=axis, pad=0.015, fraction=0.05)
        colourbar.set_label(
            "$p_{\\mathrm{rms}}/p_\\infty$"
            if quantity == "p"
            else "$u_{x,\\mathrm{rms}}/U_j$",
            fontsize=8,
            labelpad=2,
        )
        colourbar.ax.tick_params(length=2, labelsize=7.5)
        panel = f"({chr(97 + 2 * row + column)}) {mode.upper()}"
        axis.text(
            0.015,
            0.94,
            panel,
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=8.5,
            color="w",
            bbox=dict(facecolor="black", edgecolor="none", alpha=0.35, pad=1.0),
        )
        if column == 0:
            axis.set_ylabel("$D_e/\\Delta x_{\\min}$", labelpad=2)
        saturated = np.count_nonzero(values > VMAX[quantity]) / values.size
        print(f"{mode} {quantity}: fraction above Figure 3.9 limit = {saturated:.4%}")

for axis in axes[1]:
    axis.set_xlabel("$x/D_e$", labelpad=2)

fig.text(
    0.5,
    0.995,
    "$r/D_e=0.5$ (3D values are azimuthally averaged)",
    ha="center",
    va="top",
    fontsize=8.5,
)
fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.965), w_pad=1.0, h_pad=0.6)
fig.savefig(HERE / "paperfig_rmsheat_r05_2x2.png", dpi=300)
fig.savefig(HERE / "paperfig_rmsheat_r05_2x2.pdf")
print("wrote paperfig_rmsheat_r05_2x2")
