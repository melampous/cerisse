#!/usr/bin/env python3
"""Signed instantaneous vorticity omega_z De/Uj of the three schemes over
0 < x/De < 8 (3D L5, z=0 streamwise plane). Sign shows rotation direction:
the two shear-layer sheets carry opposite-signed vorticity (red/blue).
This is instantaneous omega_z on the z=0 plane (shear-layer rollup and
turbulent eddies); teno/euler t=9.300 ms, baseline t=9.424 ms.
(instantaneous -> qualitative). Diverging colour scale, equal aspect."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from pathlib import Path

plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7})
HERE = Path(__file__).parent
D, UJ = 0.0254, 414.2
XLIM, RLIM = 8.0, 1.2
XFC = 0.876
XN = [1.178, 2.479, 3.584, 4.588]
LAB = {"baseline": "WENO-Z5 NS (baseline)", "teno": "TENO5 NS",
       "euler": "Euler (WENO-Z5)"}
SCHEMES = ["baseline", "teno", "euler"]
S1 = D / UJ
VLIM = 12.0

fig, axes = plt.subplots(3, 1, figsize=(8.6, 5.4), sharex=True,
                         gridspec_kw={"hspace": 0.10})
im = None
for ax, s in zip(axes, SCHEMES):
    d = np.load(HERE / ("SCHEME_inst_%s_x24.npz" % s))
    om = d["om"] * S1
    dx = float(d["dx"])
    nx, ny = om.shape
    x = (np.arange(nx) + 0.5) * dx / D
    y = (np.arange(ny) - ny // 2 + 0.5) * dx / D
    im = ax.imshow(om.T, origin="lower", extent=[x[0], x[-1], y[0], y[-1]],
                   vmin=-VLIM, vmax=VLIM, cmap="RdBu_r", aspect="equal",
                   interpolation="nearest", rasterized=True)
    ax.axvline(XFC, color="#2ca02c", lw=0.9, ls="--", alpha=0.9)
    for xv in XN:
        ax.axvline(xv, color="0.35", lw=0.6, ls=":", alpha=0.7)
    if s == "baseline":
        tr = mtransforms.blended_transform_factory(ax.transData,
                                                   ax.transAxes)
        ax.text(XFC, 1.05, r"$x_{fc}$", transform=tr, ha="center",
                va="bottom", fontsize=8.5, color="#1a7a1a")
        for k, xv in enumerate(XN):
            ax.text(xv, 1.05, r"$x_%d$" % (k + 1), transform=tr,
                    ha="center", va="bottom", fontsize=8.5)
    ax.set_xlim(0, XLIM); ax.set_ylim(-RLIM, RLIM)
    ax.set_ylabel(r"$r/D_e$")
    ax.text(0.011, 0.07, LAB[s], transform=ax.transAxes, fontsize=8.5,
            va="bottom", color="k",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.65,
                      boxstyle="round,pad=0.25"))
axes[-1].set_xlabel(r"$x/D_e$")
fig.subplots_adjust(left=0.075, right=0.895, top=0.945, bottom=0.09)
cax = fig.add_axes([0.91, 0.20, 0.016, 0.6])
cb = fig.colorbar(im, cax=cax)
cb.set_label(r"$\omega_z(\langle u\rangle)\,D_e/U_j$", fontsize=9)
cb.ax.tick_params(labelsize=8)
fig.savefig(HERE / "paperfig_schemes_field_vort_signed_inst.png", dpi=280)
fig.savefig(HERE / "paperfig_schemes_field_vort_signed_inst.pdf", dpi=280)
print("wrote paperfig_schemes_field_vort_signed_inst")
