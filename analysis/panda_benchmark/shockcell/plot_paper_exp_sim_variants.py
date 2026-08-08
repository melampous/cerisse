#!/usr/bin/env python3
"""Paper figure: split composites (upper: experiment, lower: LES) of the
time-mean density field for three inflow profiles: baseline (delta=254 um),
wall-tanh (a=200 um), tophat (delta=0). Shared colour scale."""
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
XR = (0.6, 6.92)
V0, V1 = 0.35, 1.45

exp = np.load(HERE / "PANDA_m142_meanfield.npz")
xg = np.linspace(XR[0], XR[1], 640)
rg = np.linspace(0.0, 0.66, 140)
XG, RG = np.meshgrid(xg, rg)
Ee = griddata((exp["x"], exp["r"]), exp["mean"], (XG, RG), method="cubic")

CASES = [
    ("SIM_L5_rhoMEAN_z0.npz",       "baseline ($\\delta=254$ µm)", "(a)"),
    ("SIM_walltanh_rhoMEAN_z0.npz", "wall-tanh ($a=200$ µm)",      "(b)"),
    ("SIM_tophat_rhoMEAN_z0.npz",   "top-hat ($\\delta=0$)",        "(c)"),
]
fig, axs = plt.subplots(3, 1, figsize=(7.0, 6.6), sharex=True)
for ax, (fn, lab, tag) in zip(axs, CASES):
    sim = np.load(HERE / fn)
    xs = sim["x"]; ys = sim["y"]
    mS = (xs >= XR[0]) & (xs <= XR[1])
    mY = ys <= 0.0
    S = (sim["rho"][mS][:, mY] / RHO_J_SIM).T
    im = ax.imshow(Ee, origin="lower", extent=[XR[0], XR[1], 0, 0.66],
                   vmin=V0, vmax=V1, cmap="RdYlBu_r", aspect="equal",
                   interpolation="bilinear", rasterized=True)
    ax.imshow(S, origin="lower", extent=[XR[0], XR[1], ys[mY][0], 0],
              vmin=V0, vmax=V1, cmap="RdYlBu_r", aspect="equal",
              interpolation="nearest", rasterized=True)
    ax.axhline(0.0, color="k", lw=0.7)
    ax.set_ylim(-0.85, 0.66)
    ax.set_xlim(*XR)
    ax.set_ylabel("$r/D_e$", labelpad=1.5)
    ax.text(0.012, 0.94, f"{tag} exp / {lab}", transform=ax.transAxes,
            ha="left", va="top", fontsize=8,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1.2))
    ax.tick_params(length=2.5)
axs[-1].set_xlabel("$x/D_e$", labelpad=1.5)
cb = fig.colorbar(im, ax=axs, pad=0.015, fraction=0.03)
cb.set_label("$\\langle\\rho\\rangle/\\rho_j$", labelpad=4)
cb.ax.tick_params(length=2.5, labelsize=8)
fig.savefig(HERE / "paperfig_expsim_rho_variants.png", dpi=300,
            bbox_inches="tight")
fig.savefig(HERE / "paperfig_expsim_rho_variants.pdf", dpi=300,
            bbox_inches="tight")
print("wrote paperfig_expsim_rho_variants")
