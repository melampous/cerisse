#!/usr/bin/env python3
"""Panda M_j = 1.42 mean density vs three numerical-scheme calculations.

Same composition as plot_panda_vs_4profile_composite.py: the upper half of
each panel contains the experimental cycle-mean density for x/D_e >= 0.6,
with the unmeasured upstream part filled by the corresponding calculation
sampled on the experimental grid for visual continuity; the lower half is
the calculated time-mean density. Rows: baseline WENO-Z5 NS, TENO5 NS,
Euler (WENO-Z5). Both halves normalised by the fully expanded density."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator

HERE = Path(__file__).resolve().parent
SHOCK = HERE.parent / "shockcell"
RHO_J = 1.6413
XLIM = (0.0, 6.92)
VMIN, VMAX = 0.4, 1.45

VARIANTS = [
    ("WENO-Z5 NS (baseline)", "SIM_L5_rhoMEAN_z0.npz"),
    ("TENO5 NS", "SIM_teno_rhoMEAN_z0.npz"),
    ("Euler (WENO-Z5)", "SIM_euler_rhoMEAN_z0.npz"),
]

# baseline shock-cell metrics (3D L5 truestat, mean-pressure definitions)
XFC = 0.876
XN = [1.178, 2.479, 3.584, 4.588]
CFC = "#7fd8ff"

grid = np.load(HERE / "panda_m142_grid.npz")
xu, ru, G = grid["xu"], grid["ru"], grid["G"]      # G[r, x]
jr = ru >= 0.0
r_top = ru[jr]
G_top = G[jr, :]
dx_p = float(np.median(np.diff(xu)))
x_fill = np.arange(xu[0] - dx_p, 0.0, -dx_p)[::-1]
x_top = np.concatenate([x_fill, xu])

plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7})
fig, axes = plt.subplots(3, 1, figsize=(8.2, 6.6), sharex=True,
                         gridspec_kw={"hspace": 0.12})

im = None
for ax, (lab, fn) in zip(axes, VARIANTS):
    d = np.load(SHOCK / fn)
    sx, sy, srho = d["x"], d["y"], d["rho"] / RHO_J
    jn = sy <= 0.0
    im = ax.pcolormesh(sx, sy[jn], srho[:, jn].T, cmap="inferno",
                       vmin=VMIN, vmax=VMAX, shading="auto",
                       rasterized=True)
    interp = RegularGridInterpolator((sx, sy), srho, bounds_error=False,
                                     fill_value=np.nan)
    FX, FR = np.meshgrid(x_fill, r_top)
    fill = interp(np.column_stack([FX.ravel(), FR.ravel()])).reshape(FR.shape)
    top = np.concatenate([fill, G_top], axis=1)
    ax.pcolormesh(x_top, r_top, top, cmap="inferno", vmin=VMIN, vmax=VMAX,
                  shading="nearest", rasterized=True)
    ax.axhline(0.0, color="w", lw=0.8)
    ax.axvline(xu[0], color="w", lw=0.7, ls="--", alpha=0.9)
    ax.axvline(XFC, color=CFC, lw=0.9, ls="--", alpha=0.95)
    for xv in XN:
        ax.axvline(xv, color="w", lw=0.7, ls=":", alpha=0.95)
    ax.set_xlim(*XLIM); ax.set_ylim(-0.85, 0.68)
    ax.set_ylabel(r"$r/D_e$")
    ax.text(0.012, 0.06, f"Experiment (upper) / {lab} (lower)",
            transform=ax.transAxes, fontsize=9, va="bottom", color="w",
            bbox=dict(facecolor="0.15", edgecolor="none", alpha=0.55,
                      boxstyle="round,pad=0.25"))
axes[-1].set_xlabel(r"$x/D_e$")

import matplotlib.transforms as mtransforms
tr = mtransforms.blended_transform_factory(axes[0].transData,
                                           axes[0].transAxes)
axes[0].text(XFC, 1.04, r"$x_{fc}$", transform=tr, ha="center",
             va="bottom", fontsize=9, color="#1f4e9c")
for k, xv in enumerate(XN):
    axes[0].text(xv, 1.04, r"$x_%d$" % (k + 1), transform=tr, ha="center",
                 va="bottom", fontsize=9)

cax = fig.add_axes([0.25, 0.058, 0.5, 0.015])
cb = fig.colorbar(im, cax=cax, orientation="horizontal")
cb.set_label(r"$\langle\rho\rangle/\rho_j$", fontsize=9)
cb.ax.tick_params(labelsize=8)

fig.subplots_adjust(left=0.075, right=0.985, top=0.955, bottom=0.145)
for ext in ("png", "pdf"):
    fig.savefig(HERE / f"panda_vs_schemes_composite.{ext}", dpi=220)
print("saved panda_vs_schemes_composite.png/pdf")
