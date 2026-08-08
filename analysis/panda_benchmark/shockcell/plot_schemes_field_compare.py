#!/usr/bin/env python3
"""Whole-field scheme comparison at M_j = 1.42 (3D L5, z=0 streamwise
plane, identical grids). Three schemes stacked for direct comparison:
WENO-Z5 NS (baseline), TENO5 NS, Euler (WENO-Z5).

  paperfig_schemes_field_schlieren : instantaneous numerical schlieren
      exp(-8|grad rho|/|grad rho|_99.9), shows shock train, shear layer
      and turbulence together; teno/euler at t=9.300 ms, baseline at the
      nearest stored frame t=9.424 ms (instantaneous -> qualitative).
  paperfig_schemes_field_rhomean   : time-mean density <rho>/rho_j, shows
      the mean shock-cell train and its downstream decay.
Shared colour scale; x_fc and x_1..x_4 (baseline truestat) marked."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from pathlib import Path

import sys
plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7})
HERE = Path(__file__).parent
D, RHO_J = 0.0254, 1.6413
X24 = "--x24" in sys.argv
XLIM, RLIM = (24.0, 1.5) if X24 else (7.2, 1.2)
SFX = "_x24" if X24 else ""
XFC = 0.876
XN = [1.178, 2.479, 3.584, 4.588]
LAB = {"baseline": "WENO-Z5 NS (baseline)", "teno": "TENO5 NS",
       "euler": "Euler (WENO-Z5)"}
SCHEMES = ["baseline", "teno", "euler"]
ASP = "equal"
FIGSIZE = (11.0, 4.2) if X24 else (8.4, 5.6)


def markers(ax, top):
    ax.axvline(XFC, color="#7fd8ff", lw=0.9, ls="--", alpha=0.95)
    for xv in XN:
        ax.axvline(xv, color="w", lw=0.6, ls=":", alpha=0.7)
    if top:
        tr = mtransforms.blended_transform_factory(ax.transData,
                                                   ax.transAxes)
        ax.text(XFC, 1.04, r"$x_{fc}$", transform=tr, ha="center",
                va="bottom", fontsize=8.5, color="#1f4e9c")
        for k, xv in enumerate(XN):
            ax.text(xv, 1.04, r"$x_%d$" % (k + 1), transform=tr,
                    ha="center", va="bottom", fontsize=8.5)


# ============ Fig 1: instantaneous numerical schlieren ============
gs = {}
for s in SCHEMES:
    d = np.load(HERE / ("SCHEME_inst_%s%s.npz" % (s, SFX)))
    g, dx = d["grho"], float(d["dx"])
    nx, ny = g.shape
    x = (np.arange(nx) + 0.5) * dx / D
    y = (np.arange(ny) - ny // 2 + 0.5) * dx / D
    gs[s] = (x, y, g, float(d["time"]) * 1e3)
gref = np.nanpercentile(np.concatenate(
    [gs[s][2].ravel() for s in SCHEMES]), 99.9)

fig, axes = plt.subplots(3, 1, figsize=FIGSIZE, sharex=True,
                         gridspec_kw={"hspace": 0.10})
for ax, s in zip(axes, SCHEMES):
    x, y, g, tms = gs[s]
    ext = [x[0], x[-1], y[0], y[-1]]
    ax.imshow(np.exp(-8.0 * g.T / gref), origin="lower", extent=ext,
              vmin=0, vmax=1, cmap="gray", aspect=ASP,
              interpolation="nearest", rasterized=True)
    markers(ax, s == "baseline")
    ax.set_xlim(0, XLIM); ax.set_ylim(-RLIM, RLIM)
    ax.set_ylabel(r"$r/D_e$")
    ax.text(0.011, 0.06, "%s, $t=%.3f$ ms" % (LAB[s], tms),
            transform=ax.transAxes, fontsize=8.5, va="bottom", color="k",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.6,
                      boxstyle="round,pad=0.25"))
axes[-1].set_xlabel(r"$x/D_e$")
fig.subplots_adjust(left=0.075, right=0.985, top=0.94, bottom=0.09)
fig.savefig(HERE / ("paperfig_schemes_field_schlieren%s.png" % SFX), dpi=260)
fig.savefig(HERE / ("paperfig_schemes_field_schlieren%s.pdf" % SFX), dpi=260)
print("wrote paperfig_schemes_field_schlieren%s" % SFX)

# ============ Fig 2: mean density ============
FILES = {"baseline": "SIM_L5_rhoMEAN_z0%s.npz" % SFX,
         "teno": "SIM_teno_rhoMEAN_z0%s.npz" % SFX,
         "euler": "SIM_euler_rhoMEAN_z0%s.npz" % SFX}
fig, axes = plt.subplots(3, 1, figsize=FIGSIZE, sharex=True,
                         gridspec_kw={"hspace": 0.10})
im = None
for ax, s in zip(axes, SCHEMES):
    d = np.load(HERE / FILES[s])
    x, y, rho = d["x"], d["y"], d["rho"] / RHO_J
    im = ax.imshow(rho.T, origin="lower",
                   extent=[x[0], x[-1], y[0], y[-1]], cmap="inferno",
                   vmin=0.4, vmax=1.45, aspect="equal",
                   interpolation="nearest", rasterized=True)
    markers(ax, s == "baseline")
    ax.set_xlim(0, XLIM); ax.set_ylim(-RLIM, RLIM)
    ax.set_ylabel(r"$r/D_e$")
    ax.text(0.011, 0.06, LAB[s], transform=ax.transAxes, fontsize=8.5,
            va="bottom", color="w",
            bbox=dict(facecolor="0.15", edgecolor="none", alpha=0.55,
                      boxstyle="round,pad=0.25"))
axes[-1].set_xlabel(r"$x/D_e$")
fig.subplots_adjust(left=0.075, right=0.90, top=0.94, bottom=0.09)
cax = fig.add_axes([0.915, 0.20, 0.016, 0.6])
cb = fig.colorbar(im, cax=cax)
cb.set_label(r"$\langle\rho\rangle/\rho_j$", fontsize=9)
cb.ax.tick_params(labelsize=8)
fig.savefig(HERE / ("paperfig_schemes_field_rhomean%s.png" % SFX), dpi=260)
fig.savefig(HERE / ("paperfig_schemes_field_rhomean%s.pdf" % SFX), dpi=260)
print("wrote paperfig_schemes_field_rhomean%s" % SFX)
