#!/usr/bin/env python3
"""Instantaneous structural comparison of the four nozzle inlet profiles
(z=0 streamwise plane): numerical schlieren, vorticity magnitude
|omega| De/Uj, and max(-lambda_{2,d},0) (De/Uj)^2. Rows = profiles,
columns = the three fields. Latest available plotfile of each case.

NOTE: baseline/wall-tanh/thin-shifted are L5 (De/256) at t~10.2-10.4 ms;
the top-hat L5 run failed, so top-hat is shown at L3 (De/64) and
t=5.0 ms -- its smoother small-scale content is a resolution/time
artefact, not a property of the profile. Each panel is labelled with
its level and time."""
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
SHOCK = HERE.parent / "shockcell"
D, UJ = 0.0254, 414.2
S1, S2 = D / UJ, (D / UJ)**2
XW, RW = 8.0, 1.5
# (tag, label, level)
CASES = [
    ("baseline", "Baseline shifted tanh", "L5"),
    ("tophat", "Top hat", "L3"),
    ("walltanh", "Wall tanh", "L5"),
    ("shift100", "Thin shifted tanh", "L5"),
]

data = {}
for tag, _, _ in CASES:
    d = np.load(SHOCK / ("STRUCT4_%s.npz" % tag))
    data[tag] = d
gref = np.nanpercentile(np.concatenate(
    [data[t]["grho"].ravel() for t, _, _ in CASES]), 99.9)

fig, axs = plt.subplots(4, 3, figsize=(8.6, 7.0), sharex=True,
                        sharey=True)
HEADS = ["numerical schlieren", "$|\\omega|\\,D_e/U_j$",
         "$\\max(-\\lambda_{2,d},0)\\,(D_e/U_j)^2$"]
ims = {}
for rI, (tag, lab, lev) in enumerate(CASES):
    d = data[tag]
    x, y = d["x"], d["y"]
    tms = float(d["time"]) * 1e3
    ext = [x[0], x[-1], y[0], y[-1]]
    panels = [
        (np.exp(-8 * d["grho"] / gref), 0, 1, "gray"),
        (d["omag"] * S1, 0, 12, "viridis"),
        (np.maximum(-d["lam2"], 0.0) * S2, 0, 30, "magma"),
    ]
    for cI, (F, v0, v1, cm) in enumerate(panels):
        ax = axs[rI, cI]
        ims[cI] = ax.imshow(F.T, origin="lower", extent=ext, vmin=v0,
                            vmax=v1, cmap=cm, aspect="equal",
                            interpolation="nearest", rasterized=True)
        if rI == 0:
            ax.set_title(HEADS[cI], fontsize=8.5, pad=3)
        tcol = "k" if cm == "gray" else "w"
        ax.text(0.015, 0.94, "(%s)" % chr(97 + 3 * rI + cI),
                transform=ax.transAxes, ha="left", va="top", fontsize=7.5,
                color=tcol,
                bbox=dict(facecolor="white" if tcol == "k" else "black",
                          edgecolor="none", alpha=0.5, pad=1.0))
        ax.set_xlim(0, XW); ax.set_ylim(-RW, RW)
        ax.tick_params(length=2.2, labelsize=7.5)
    axs[rI, 0].set_ylabel("$r/D_e$", labelpad=1)
    axs[rI, 0].text(-0.34, 0.5,
                    "%s\n(%s, $t=%.1f$ ms)" % (lab, lev, tms),
                    transform=axs[rI, 0].transAxes, rotation=90,
                    ha="center", va="center", fontsize=7.5)
for ax in axs[-1]:
    ax.set_xlabel("$x/D_e$", labelpad=2)
fig.subplots_adjust(left=0.10, right=0.985, top=0.955, bottom=0.115,
                    wspace=0.06, hspace=0.06)
for x0, key, lab in ((0.12, 0, "schlieren"),
                     (0.44, 1, "$|\\omega|\\,D_e/U_j$"),
                     (0.74, 2, "$\\max(-\\lambda_{2,d},0)\\,(D_e/U_j)^2$")):
    cax = fig.add_axes([x0, 0.05, 0.22, 0.014])
    cb = fig.colorbar(ims[key], cax=cax, orientation="horizontal")
    cb.set_label(lab, fontsize=7.5, labelpad=1)
    cb.ax.tick_params(length=2, labelsize=7)
fig.savefig(HERE / "paperfig_var4_struct.png", dpi=270)
fig.savefig(HERE / "paperfig_var4_struct.pdf", dpi=270)
print("wrote paperfig_var4_struct")
