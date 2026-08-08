#!/usr/bin/env python3
"""Plot the scaled azimuthal-vorticity diagnostic for 2D and 3D fields."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
import numpy as np


HERE = Path(__file__).resolve().parent
DATA = np.load(HERE / "OMEGA_RHOR_NATIVE.npz")

plt.rcParams.update({
    "font.size": 9,
    "axes.linewidth": 0.6,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
})


def local_axis_mask(radius, levels):
    """Mask the two cells nearest the axis at each local AMR level."""
    local_spacing_over_diameter = 1.0 / (
        8.0 * 2.0 ** levels.astype(int)
    )
    return radius < 2.0 * local_spacing_over_diameter


r2 = DATA["r2"]
x2 = DATA["x2"]
y3 = DATA["y3"]
x3 = DATA["x3"]
q2 = DATA["q2"].astype(float)
q3 = DATA["q3"].astype(float)
level2 = DATA["level2"]
level3 = DATA["level3"]

q2 = np.ma.masked_where(
    local_axis_mask(r2[:, None], level2),
    q2,
)
q3 = np.ma.masked_where(
    local_axis_mask(np.abs(y3)[None, :], level3),
    q3,
)

# The 2D scalar is even under reflection about the axis.
r2_mirrored = np.concatenate([-r2[::-1], r2])
q2_mirrored = np.ma.concatenate(
    [q2[::-1, :], q2],
    axis=0,
).T

norm = SymLogNorm(
    linthresh=1.0,
    linscale=0.8,
    vmin=-200.0,
    vmax=200.0,
    base=10.0,
)
cmap = plt.get_cmap("RdBu_r").copy()
cmap.set_bad("white")

fig, axes = plt.subplots(
    2,
    1,
    figsize=(6.5, 4.35),
    sharex=True,
    sharey=True,
)

panels = [
    (
        axes[0],
        q2_mirrored,
        x2,
        r2_mirrored,
        r"(a) Axisymmetric L4",
    ),
    (
        axes[1],
        q3,
        x3,
        y3,
        r"(b) Cartesian L5 ($z=0$)",
    ),
]

image = None
for axis, field, x, radius, label in panels:
    image = axis.imshow(
        field.T,
        origin="lower",
        extent=[x[0], x[-1], radius[0], radius[-1]],
        cmap=cmap,
        norm=norm,
        aspect="equal",
        interpolation="nearest",
        rasterized=True,
    )
    axis.text(
        0.015,
        0.94,
        label,
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=8.0,
        bbox={
            "facecolor": "white",
            "edgecolor": "none",
            "alpha": 0.72,
            "pad": 1.4,
        },
    )
    axis.set_xlim(0.0, 4.2)
    axis.set_ylim(-1.25, 1.25)
    axis.set_ylabel(r"$y/D_e$", labelpad=2)
    axis.tick_params(length=2.3, labelsize=8)

axes[-1].set_xlabel(r"$x/D_e$", labelpad=2)

fig.subplots_adjust(
    left=0.105,
    right=0.985,
    top=0.985,
    bottom=0.185,
    hspace=0.08,
)
cbar_axis = fig.add_axes([0.25, 0.075, 0.50, 0.025])
cbar = fig.colorbar(
    image,
    cax=cbar_axis,
    orientation="horizontal",
    ticks=[-100, -10, -1, 0, 1, 10, 100],
    extend="both",
)
cbar.set_label(
    r"$\mathcal{Q}_{\theta}"
    r"=(\rho_jD_e^2/U_j)\,\omega_{\theta}/(\rho r)$",
    fontsize=8.4,
    labelpad=1,
)
cbar.ax.tick_params(length=2, labelsize=7.5)

for suffix in ("png", "pdf"):
    fig.savefig(
        HERE / f"paperfig_dim_omega_rhor.{suffix}",
        dpi=300,
    )
print("wrote paperfig_dim_omega_rhor")
