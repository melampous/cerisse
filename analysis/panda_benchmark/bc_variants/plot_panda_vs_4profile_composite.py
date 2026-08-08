#!/usr/bin/env python3
"""Panda M_j = 1.42 mean density and four exit-profile calculations.

The upper half of each panel contains the experimental cycle-mean density for
x/D_e >= 0.6. The unmeasured upstream part is filled by the corresponding
calculation sampled on the experimental grid for visual continuity. The lower
half contains the calculated time-mean density. Both are normalised by their
corresponding fully expanded jet density.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from scipy.interpolate import RegularGridInterpolator

HERE = Path(__file__).resolve().parent
SHOCK = HERE.parent / "shockcell"
SCR = HERE
RHO_J = 1.6413
XLIM = (0.0, 6.92)
VMIN, VMAX = 0.4, 1.45

# (label, simulation field)
VARIANTS = [
    ("Baseline shifted tanh", "SIM_L5_rhoMEAN_z0.npz"),
    ("Top hat", "SIM_tophat_rhoMEAN_z0.npz"),
    ("Wall-attached tanh", "SIM_walltanh_rhoMEAN_z0.npz"),
    ("Thin shifted tanh", "SIM_shift100_rhoMEAN_z0.npz"),
]

# ---- Panda native grid (r>=0 half) for the top; x<0.6 is unmeasured and is
# filled with the simulation data downsampled to Panda's coarse resolution.
grid = np.load(SCR / "panda_m142_grid.npz")
xu, ru, G = grid["xu"], grid["ru"], grid["G"]      # G[r, x]
jr = ru >= 0.0
r_top = ru[jr]                                       # 0 .. 0.70, dr~0.07
G_top = G[jr, :]                                     # (nr>=0, 80)
dx_p = float(np.median(np.diff(xu)))                 # ~0.08
x_fill = np.arange(xu[0] - dx_p, 0.0, -dx_p)[::-1]   # coarse cols for x<0.6
x_top = np.concatenate([x_fill, xu])

# shock-cell marks: x_fc (first closure) and x_1..x_5 (compression maxima)
MARKS = np.load(HERE / "shockcell_marks.npz")
MKEY = {"SIM_L5_rhoMEAN_z0.npz": "baseline", "SIM_tophat_rhoMEAN_z0.npz": "tophat",
        "SIM_walltanh_rhoMEAN_z0.npz": "walltanh",
        "SIM_shift100_rhoMEAN_z0.npz": "shift100"}
STROKE = [pe.withStroke(linewidth=1.4, foreground="k", alpha=0.85)]
CMAX = "white"


def draw_marks(ax, xfc, xmax, r_star):
    """x_fc as a white star, x_1..x_5 as white circles (r = r_star)."""
    ax.plot(xmax, np.full(len(xmax), r_star), marker="o", ls="", mfc=CMAX,
            mec="k", mew=0.7, ms=6.5, zorder=7, path_effects=STROKE)
    ax.plot([xfc], [r_star], marker="*", ls="", mfc=CMAX, mec="k", mew=0.6,
            ms=11, zorder=8, path_effects=STROKE)

plt.rcParams.update({
    "font.size": 11.2,
    "axes.labelsize": 11.4,
    "xtick.labelsize": 10.4,
    "ytick.labelsize": 10.4,
    "axes.linewidth": 0.7,
})
fig, axes = plt.subplots(4, 1, figsize=(7.2, 7.6), sharex=True,
                         gridspec_kw={"hspace": 0.22})

im = None
for ax, (lab, fn) in zip(axes, VARIANTS):
    d = np.load(SHOCK / fn)
    sx, sy, srho = d["x"], d["y"], d["rho"] / RHO_J
    # bottom half: r<=0 side of the simulation (full resolution)
    jn = sy <= 0.0
    sy_n, srho_n = sy[jn], srho[:, jn]
    im = ax.pcolormesh(sx, sy_n, srho_n.T, cmap="inferno", vmin=VMIN, vmax=VMAX,
                       shading="auto", rasterized=True)
    # top half: sim fill for x<0.6 (Panda coarse res) + Panda exp for x>=0.6
    interp = RegularGridInterpolator((sx, sy), srho, bounds_error=False,
                                     fill_value=np.nan)
    FX, FR = np.meshgrid(x_fill, r_top)
    fill = interp(np.column_stack([FX.ravel(), FR.ravel()])).reshape(FR.shape)
    top = np.concatenate([fill, G_top], axis=1)       # (nr, len(x_top))
    ax.pcolormesh(x_top, r_top, top, cmap="inferno", vmin=VMIN, vmax=VMAX,
                  shading="nearest", rasterized=True)
    ax.axhline(0.0, color="w", lw=0.8)
    ax.axvline(xu[0], color="w", lw=0.7, ls="--", alpha=0.9)
    # markers just off the axis: x_fc (star) + x_1..x_5 (circles)
    key = MKEY[fn]
    draw_marks(ax, float(MARKS[f"{key}_xfc"]), MARKS[f"{key}_xmax"], -0.06)
    draw_marks(ax, float(MARKS["experiment_xfc"]), MARKS["experiment_xmax"], 0.06)
    ax.set_xlim(*XLIM); ax.set_ylim(-0.85, 0.68)
    ax.set_ylabel(r"$r/D_e$")
    ax.set_title(f"Experiment (upper) / {lab} (lower)",
                 loc="left", fontsize=10.5, pad=2.0)
axes[-1].set_xlabel(r"$x/D_e$")

cax = fig.add_axes([0.25, 0.035, 0.5, 0.012])
cb = fig.colorbar(im, cax=cax, orientation="horizontal")
cb.set_label(r"$\langle\rho\rangle/\rho_j$", fontsize=11.0)
cb.ax.tick_params(labelsize=10.0)

fig.subplots_adjust(left=0.075, right=0.985, top=0.985, bottom=0.135)
for ext in ("png", "pdf"):
    fig.savefig(
        HERE / f"panda_vs_4profile_composite.{ext}",
        dpi=240,
        bbox_inches="tight",
        pad_inches=0.03,
    )
print("saved panda_vs_4profile_composite.png/pdf")
