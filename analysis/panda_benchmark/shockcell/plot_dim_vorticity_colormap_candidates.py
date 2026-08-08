#!/usr/bin/env python3
"""Render Figure 3.20 colour-map candidates without replacing thesis assets."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


plt.rcParams.update({
    "font.size": 9,
    "axes.linewidth": 0.6,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
})

HERE = Path(__file__).parent
D, UJ = 0.0254, 414.2
S1 = D / UJ
XW, RW = 4.2, 1.25

d3 = np.load(HERE / "VBUD3D_L5.npz")
d2 = np.load(HERE / "VBUD2D_L4.npz")
x3, y3 = d3["x"], d3["y"]
r2, x2 = d2["r"], d2["xax"]


def mirror2d(field):
    """Mirror a scalar axisymmetric field from (r, x) to (x, +/-r)."""
    positive = field.T
    return np.concatenate([positive[:, ::-1], positive], axis=1)


y2 = np.concatenate([-r2[::-1], r2])
sch2 = mirror2d(d2["grho"])
sch3 = d3["grho"]
gref = np.nanpercentile(np.concatenate([sch2.ravel(), sch3.ravel()]), 99.9)
omth2 = mirror2d(d2["omth"]) * S1
omth3 = d3["omz"] * np.where(y3 >= 0.0, 1.0, -1.0)[None, :] * S1


def render(name, cmap, limit):
    fig, axs = plt.subplots(
        2, 2, figsize=(6.5, 4.4), sharex=True, sharey=True
    )
    rows = [
        ((np.exp(-8 * sch2 / gref), 0, 1, "gray"),
         (omth2, -limit, limit, cmap), x2, y2),
        ((np.exp(-8 * sch3 / gref), 0, 1, "gray"),
         (omth3, -limit, limit, cmap), x3, y3),
    ]
    heads = ["Numerical schlieren", r"$\omega_\theta D_e/U_j$"]
    images = {}
    for row, (*panels, xx, yy) in enumerate(rows):
        for col, (field, vmin, vmax, colour_map) in enumerate(panels):
            ax = axs[row, col]
            images[col] = ax.imshow(
                field.T,
                origin="lower",
                extent=[xx[0], xx[-1], yy[0], yy[-1]],
                vmin=vmin,
                vmax=vmax,
                cmap=colour_map,
                aspect="equal",
                interpolation="nearest",
                rasterized=True,
            )
            if row == 0:
                ax.set_title(heads[col], fontsize=8.5, pad=3)
            ax.text(
                0.02,
                0.95,
                f"({chr(97 + 2 * row + col)})",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=8,
                color="k",
                bbox={
                    "facecolor": "white",
                    "edgecolor": "none",
                    "alpha": 0.5,
                    "pad": 1.0,
                },
            )
            ax.set_xlim(0, XW)
            ax.set_ylim(-RW, RW)
            ax.tick_params(length=2.2, labelsize=7.5)
    for row in range(2):
        axs[row, 0].set_ylabel(r"$y/D_e$", labelpad=1)
    for ax in axs[1]:
        ax.set_xlabel(r"$x/D_e$", labelpad=2)
    fig.subplots_adjust(
        left=0.13,
        right=0.985,
        top=0.93,
        bottom=0.20,
        wspace=0.08,
        hspace=0.06,
    )
    for x0, key, label in (
        (0.13, 0, "Numerical schlieren"),
        (0.57, 1, r"$\omega_\theta D_e/U_j$"),
    ):
        cax = fig.add_axes([x0, 0.085, 0.34, 0.02])
        cb = fig.colorbar(images[key], cax=cax, orientation="horizontal")
        cb.set_label(label, fontsize=8, labelpad=1)
        cb.ax.tick_params(length=2, labelsize=7.5)
    stem = HERE / f"paperfig_dim_vortcomp_candidate_{name}"
    fig.savefig(stem.with_suffix(".png"), dpi=280)
    fig.savefig(stem.with_suffix(".pdf"), dpi=280)
    plt.close(fig)
    print(f"wrote {stem.name}")


# Existing RdBu colour language with tighter limits.
render("rdbu12", "RdBu_r", 12.0)
# Colour-vision-friendly diverging alternative with equally tight limits.
render("puor12", "PuOr_r", 12.0)
# Higher-contrast red/blue alternative with slightly wider limits.
render("seismic15", "seismic", 15.0)
