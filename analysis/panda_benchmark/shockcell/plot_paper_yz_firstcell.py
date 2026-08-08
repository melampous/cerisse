#!/usr/bin/env python3
"""Paper figure: yz cross-sections of <rho>/rho_j inside the first shock cell
(x/D = 0.2, 0.4, 0.6, 0.8; upstream of the axis closure at x_s = 0.88).
Baseline 3D l_max = 5, [6.5, 10.16] ms window. 1x4 panels, shared scale."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams.update({
    "font.size": 9, "axes.linewidth": 0.6,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
})
HERE = Path(__file__).parent
RHO_J = 1.6413
d = np.load(HERE / "SIM_L5_rhoMEAN_yz4.npz")
y, z = d["y"], d["z"]
STATIONS = [0.2, 0.4, 0.6, 0.8]
V0, V1 = 0.35, 1.45
R = 0.8

fig, axs = plt.subplots(1, 4, figsize=(10.4, 2.95), sharey=True)
for k, (xd, ax) in enumerate(zip(STATIONS, axs)):
    key = "x" + str(xd).replace(".", "p")
    F = d[key] / RHO_J                       # (y, z)
    im = ax.imshow(F, origin="lower", extent=[z[0], z[-1], y[0], y[-1]],
                   vmin=V0, vmax=V1, cmap="RdYlBu_r", aspect="equal",
                   interpolation="nearest", rasterized=True)
    th = np.linspace(0, 2*np.pi, 200)
    ax.plot(0.5*np.cos(th), 0.5*np.sin(th), color="k", lw=0.5, ls=":")
    ax.set_xlim(-R, R); ax.set_ylim(-R, R)
    ax.set_xlabel("$z/D_e$", labelpad=1.5)
    ax.text(0.03, 0.955, f"({chr(97+k)}) $x/D_e={xd}$", transform=ax.transAxes,
            ha="left", va="top", fontsize=8.5,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1.2))
    ax.tick_params(length=2.5)
axs[0].set_ylabel("$y/D_e$", labelpad=1.5)
cb = fig.colorbar(im, ax=axs, pad=0.012, fraction=0.025)
cb.set_label("$\\langle\\rho\\rangle/\\rho_j$", labelpad=4)
cb.ax.tick_params(length=2.5, labelsize=8)
fig.savefig(HERE / "paperfig_yz_firstcell.png", dpi=300, bbox_inches="tight")
fig.savefig(HERE / "paperfig_yz_firstcell.pdf", dpi=300, bbox_inches="tight")
print("wrote paperfig_yz_firstcell")
