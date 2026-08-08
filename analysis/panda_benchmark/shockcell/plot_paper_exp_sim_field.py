#!/usr/bin/env python3
"""Paper figure: split-composite time-mean density field, m142.
Upper half (r > 0): experiment — 37-phase mean of the Panda & Seasholtz
EPAPS M142DEN survey (rho/rho_j, their calibration rho_j = 1.653).
Lower half (r < 0): present LES, 3D l_max = 5, <rho>/rho_j on the z = 0
plane ([6.5, 10.16] ms window). Shared colormap and scale."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from pathlib import Path

plt.rcParams.update({
    "font.size": 9, "axes.linewidth": 0.6,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
})
HERE = Path(__file__).parent
RHO_J_SIM = 1.6413

exp = np.load(HERE / "PANDA_m142_meanfield.npz")
sim = np.load(HERE / "SIM_L5_rhoMEAN_z0.npz")

XR = (0.6, 6.92)
# experiment: interpolate scattered-regular points to a fine grid, keep r>=0
xg = np.linspace(XR[0], XR[1], 640)
rg = np.linspace(0.0, 0.66, 140)
XG, RG = np.meshgrid(xg, rg)
Ee = griddata((exp["x"], exp["r"]), exp["mean"], (XG, RG), method="cubic")
# simulation: lower half (y<0)
xs = sim["x"]; ys = sim["y"]
mS = (xs >= XR[0]) & (xs <= XR[1])
mY = ys <= 0.0
S = (sim["rho"][mS][:, mY] / RHO_J_SIM).T

V0, V1 = 0.35, 1.45
fig, ax = plt.subplots(figsize=(7.0, 2.7))
im1 = ax.imshow(Ee, origin="lower", extent=[XR[0], XR[1], 0, 0.66],
                vmin=V0, vmax=V1, cmap="RdYlBu_r", aspect="equal",
                interpolation="bilinear", rasterized=True)
ax.imshow(S, origin="lower", extent=[XR[0], XR[1], ys[mY][0], 0],
          vmin=V0, vmax=V1, cmap="RdYlBu_r", aspect="equal",
          interpolation="nearest", rasterized=True)
ax.axhline(0.0, color="k", lw=0.7)
ax.set_ylim(-0.85, 0.66)
ax.set_xlim(*XR)
ax.set_xlabel("$x/D_e$", labelpad=1.5)
ax.set_ylabel("$r/D_e$", labelpad=1.5)
ax.text(0.012, 0.955, "experiment (Panda & Seasholtz 1999)",
        transform=ax.transAxes, ha="left", va="top", fontsize=8,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=1.2))
ax.text(0.012, 0.06, "present LES (L5)",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=8,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=1.2))
cb = fig.colorbar(im1, ax=ax, pad=0.015, fraction=0.045)
cb.set_label("$\\langle\\rho\\rangle/\\rho_j$", labelpad=4)
cb.ax.tick_params(length=2.5, labelsize=8)
ax.tick_params(length=2.5)
fig.tight_layout()
fig.savefig(HERE / "paperfig_expsim_rho_m142.png", dpi=300)
fig.savefig(HERE / "paperfig_expsim_rho_m142.pdf", dpi=300)
print("wrote paperfig_expsim_rho_m142")
