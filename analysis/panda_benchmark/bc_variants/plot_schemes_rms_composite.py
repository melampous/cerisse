#!/usr/bin/env python3
"""Resolved density-fluctuation fields for the three numerical schemes:
rho_rms/rho_j on the z = 0 plane (full +-y), rows = baseline WENO-Z5 NS,
TENO5 NS, Euler (WENO-Z5). A more dissipative discretisation suppresses
the resolved fluctuation field; differences are directly visible in the
shear-layer band and its downstream filling. Common colour scale."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
SHOCK = HERE.parent / "shockcell"
RHO_J = 1.6413
XLIM = (0.0, 7.2)
XFC = 0.876
XN = [1.178, 2.479, 3.584, 4.588]

VARIANTS = [
    ("WENO-Z5 NS (baseline)", "SIM_L5_rhoRMS_z0.npz"),
    ("TENO5 NS", "SIM_teno_rhoRMS_z0.npz"),
    ("Euler (WENO-Z5)", "SIM_euler_rhoRMS_z0.npz"),
]

plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7})
fig, axes = plt.subplots(3, 1, figsize=(8.2, 6.2), sharex=True,
                         gridspec_kw={"hspace": 0.12})

fields = []
for _, fn in VARIANTS:
    d = np.load(SHOCK / fn)
    fields.append((d["x"], d["y"], d["rms"] / RHO_J))
VMAX = float(np.ceil(np.nanpercentile(
    np.concatenate([f[2].ravel() for f in fields]), 99.5) * 20) / 20)

im = None
for ax, (lab, _), (sx, sy, rms) in zip(axes, VARIANTS, fields):
    im = ax.pcolormesh(sx, sy, rms.T, cmap="magma", vmin=0.0, vmax=VMAX,
                       shading="auto", rasterized=True)
    ax.axvline(XFC, color="#7fd8ff", lw=0.9, ls="--", alpha=0.95)
    for xv in XN:
        ax.axvline(xv, color="w", lw=0.6, ls=":", alpha=0.8)
    ax.set_xlim(*XLIM); ax.set_ylim(-0.85, 0.85)
    ax.set_ylabel(r"$r/D_e$")
    ax.text(0.012, 0.06, lab, transform=ax.transAxes, fontsize=9,
            va="bottom", color="w",
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
cb.set_label(r"$\rho_{\mathrm{rms}}/\rho_j$", fontsize=9)
cb.ax.tick_params(labelsize=8)

fig.subplots_adjust(left=0.075, right=0.985, top=0.95, bottom=0.15)
for ext in ("png", "pdf"):
    fig.savefig(HERE / f"schemes_rhorms_composite.{ext}", dpi=220)
print("saved schemes_rhorms_composite.png/pdf  vmax=%.2f" % VMAX)
