#!/usr/bin/env python3
"""Mean-flow vorticity |omega_z(<u>)| De/Uj of the three schemes over
0 < x/De < 24 (3D L5, z=0 streamwise plane). This is the vorticity of the
TIME-MEAN velocity field (highlights the mean shear-layer sheet and its
spreading); it is NOT the time-mean of the instantaneous vorticity
magnitude <|omega|> (unavailable -- vorticity second moments are not
stored). Three schemes stacked, common colour scale, true equal aspect."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from pathlib import Path

plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7})
HERE = Path(__file__).parent
D, UJ = 0.0254, 414.2
XLIM, RLIM = 24.0, 1.5
XFC = 0.876
XN = [1.178, 2.479, 3.584, 4.588]
LAB = {"baseline": "WENO-Z5 NS (baseline)", "teno": "TENO5 NS",
       "euler": "Euler (WENO-Z5)"}
SCHEMES = ["baseline", "teno", "euler"]
S1 = D / UJ

fig, axes = plt.subplots(3, 1, figsize=(11.0, 4.2), sharex=True,
                         gridspec_kw={"hspace": 0.10})
im = None
for ax, s in zip(axes, SCHEMES):
    d = np.load(HERE / ("SCHEME_meanvort_%s_x24.npz" % s))
    om = np.abs(d["om_mean"]) * S1
    dx = float(d["dx"])
    nx, ny = om.shape
    x = (np.arange(nx) + 0.5) * dx / D
    y = (np.arange(ny) - ny // 2 + 0.5) * dx / D
    im = ax.imshow(om.T, origin="lower", extent=[x[0], x[-1], y[0], y[-1]],
                   vmin=0, vmax=8, cmap="viridis", aspect="equal",
                   interpolation="nearest", rasterized=True)
    ax.axvline(XFC, color="#7fd8ff", lw=0.9, ls="--", alpha=0.95)
    for xv in XN:
        ax.axvline(xv, color="w", lw=0.6, ls=":", alpha=0.6)
    if s == "baseline":
        tr = mtransforms.blended_transform_factory(ax.transData,
                                                   ax.transAxes)
        ax.text(XFC, 1.05, r"$x_{fc}$", transform=tr, ha="center",
                va="bottom", fontsize=8.5, color="#1f4e9c")
        for k, xv in enumerate(XN):
            ax.text(xv, 1.05, r"$x_%d$" % (k + 1), transform=tr,
                    ha="center", va="bottom", fontsize=8.5)
    ax.set_xlim(0, XLIM); ax.set_ylim(-RLIM, RLIM)
    ax.set_ylabel(r"$r/D_e$")
    ax.text(0.011, 0.08, LAB[s], transform=ax.transAxes, fontsize=8.5,
            va="bottom", color="w",
            bbox=dict(facecolor="0.15", edgecolor="none", alpha=0.55,
                      boxstyle="round,pad=0.25"))
axes[-1].set_xlabel(r"$x/D_e$")
fig.subplots_adjust(left=0.065, right=0.90, top=0.93, bottom=0.11)
cax = fig.add_axes([0.915, 0.20, 0.014, 0.55])
cb = fig.colorbar(im, cax=cax)
cb.set_label(r"$|\omega_z(\langle u\rangle)|\,D_e/U_j$", fontsize=9)
cb.ax.tick_params(labelsize=8)
fig.savefig(HERE / "paperfig_schemes_field_meanvort_x24.png", dpi=260)
fig.savefig(HERE / "paperfig_schemes_field_meanvort_x24.pdf", dpi=260)
print("wrote paperfig_schemes_field_meanvort_x24")
