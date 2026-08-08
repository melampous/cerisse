#!/usr/bin/env python3
"""Radial mean-density comparison for the WALL-TANH profile (3D L5) at
eight stations, same composition as paperfig_cmp_radial_rho (baseline):
M_j=1.42 survey filled black (primary), M_j=1.44 traverse open grey
(mirrored, secondary), two diametral cuts y-cut (red) and z-cut (blue)
as equal-period half-period-offset dashes, +-0.05 D_e envelope."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from pathlib import Path

plt.rcParams.update({
    "font.size": 9, "axes.linewidth": 0.6,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
})
HERE = Path(__file__).parent
D = 0.0254
RHO_J = 1.6413
C3D, C2D = "#b3251e", "#1f4e9c"


def boxfilt(a, x, wD):
    w = max(3, int(round(wD / (x[1] - x[0]))))
    if w % 2 == 0:
        w += 1
    out = np.convolve(a, np.ones(w) / w, mode="same")
    out[:w//2] = a[:w//2]; out[-(w//2):] = a[-(w//2):]
    return out


f5 = np.load(HERE.parent / "panda_fig5a_mj144_EPAPS.npz", allow_pickle=True)
RHOJ_EXP = float(f5["rhoj"])
fm = np.load(HERE.parent / "panda_m142den_field_EPAPS.npz")
xm, rm, Fm = fm["x"], fm["r"], fm["rho"]
sly = np.load(HERE / "SIM_walltanh_rhoMEAN_z0.npz")   # y-cut (z=0)
slz = np.load(HERE / "SIM_walltanh_rhoMEAN_y0.npz")   # z-cut (y=0)
xs3, rr3 = sly["x"], sly["y"]
cutY = sly["rho"] / RHO_J
cutZ = slz["rho"] / RHO_J

STATIONS = [0.6, 0.75, 0.9, 1.05, 1.25, 1.55, 2.0, 3.05]
OFF = 0.05

fig, axs = plt.subplots(2, 4, figsize=(7.2, 4.15), sharex=True, sharey=True)
for k, (stn, ax) in enumerate(zip(STATIONS, axs.flat)):
    ax.set_axisbelow(True)
    ax.grid(True, ls="--", lw=0.5, color="0.87")
    e = f5["x%s" % str(stn)]
    rex = e[:, 0]; rhex = e[:, 1] / RHOJ_EXP
    icol = np.argmin(np.abs(xm - stn))
    xcol = xm[icol]; rho42 = Fm[icol]
    jsel = np.abs(xs3 - stn) <= OFF + 1e-9
    stack = np.vstack([cutY[jsel], cutZ[jsel]])
    ax.fill_between(rr3, stack.min(0), stack.max(0), color=C3D, alpha=0.15,
                    lw=0, zorder=2)
    j = np.argmin(np.abs(xs3 - stn))
    ax.plot(rr3, cutY[j], color=C3D, lw=1.1, ls=(0, (4, 4)), zorder=4)
    ax.plot(rr3, cutZ[j], color=C2D, lw=1.1, ls=(4, (4, 4)), zorder=4)
    ax.plot(rex, rhex, "o", color="0.72", ms=2.5, mfc="none", mew=0.7,
            zorder=3)
    ax.plot(-rex, rhex, "o", color="0.72", ms=2.5, mfc="none", mew=0.7,
            zorder=3)
    ax.plot(rm, rho42, "o", color="k", ms=3.2, mew=0, zorder=5)
    ax.text(0.05, 0.06, "(%s) $x/D_e=%g$" % (chr(97 + k), stn),
            transform=ax.transAxes, ha="left", va="bottom", fontsize=7.5,
            zorder=6,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.75,
                      pad=1.2))
    ax.set_xlim(-0.85, 0.85)
    ax.set_xticks([-0.5, 0.0, 0.5])
    ax.tick_params(length=2.5)
for ax in axs[1]:
    ax.set_xlabel("$r/D_e$", labelpad=2)
for ax in axs[:, 0]:
    ax.set_ylabel("$\\langle\\rho\\rangle/\\rho_j$", labelpad=2)
axs[0, 0].set_ylim(0.3, 1.5)
handles = [
    Line2D([], [], marker="o", ls="", color="k", ms=3.2,
           label="Exp. $M_j=1.42$ survey"),
    Line2D([], [], marker="o", ls="", color="0.72", mfc="none", mew=0.7,
           ms=2.5, label="Exp. $M_j=1.44$ traverse (mirrored)"),
    Line2D([], [], color=C3D, lw=1.1, ls=(0, (4, 4)),
           label="Wall tanh L5, $y$-cut"),
    Line2D([], [], color=C2D, lw=1.1, ls=(4, (4, 4)),
           label="Wall tanh L5, $z$-cut"),
    Patch(facecolor=C3D, alpha=0.15, label="$x\\pm0.05\\,D_e$ envelope"),
]
fig.legend(handles=handles, fontsize=6.5, frameon=False, ncol=5,
           loc="upper center", bbox_to_anchor=(0.5, 1.0),
           handlelength=1.3, columnspacing=0.8, handletextpad=0.4)
fig.tight_layout(w_pad=0.5, h_pad=0.5, rect=[0, 0, 1, 0.955])
fig.savefig(HERE / "paperfig_cmp_radial_walltanh.png", dpi=300)
fig.savefig(HERE / "paperfig_cmp_radial_walltanh.pdf")
print("wrote paperfig_cmp_radial_walltanh")
