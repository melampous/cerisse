#!/usr/bin/env python3
"""Figure 3.17 -- time-averaged centreline density for the four prescribed
nozzle exit profiles at M_j = 1.42 (Baseline / Top hat / Wall-attached tanh / Thin
shifted tanh). Final cumulative mean of each case, centrelines extracted to
8 De at D_e/256, 0.05 D_e display smoothing, rho_j = 1.6413.

Two panels:
 (a) overview of the four profiles + Panda (1999) experiment, distinct line
     styles (not colour alone);
 (b) deviation of each profile from the baseline profile -- removes the
     common shock-cell oscillation so the profile-to-profile phase shift is
     visible instead of four curves overlapping. Panda measures density only,
     so no pressure/Mach panel is shown.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
D_E, RHO_J = 0.0254, 1.6413
XMAX = 8.0
CB = "#aab0b7"                                   # baseline (light reference)
VARIANTS = [                                     # variants (bright)
    ("tophat", "Top hat", "#0b6ef5", (0, (5, 2))),
    ("walltanh", "Wall-attached tanh", "#2ca02c", (0, (3, 1.4, 1, 1.4))),
    ("shift100", "Thin shifted tanh", "#ff6a00", (0, (1, 1.4))),
]


def boxfilt(a, x, wD=0.05):
    w = max(3, int(round(wD / (x[1] - x[0]))))
    if w % 2 == 0:
        w += 1
    out = np.convolve(a, np.ones(w) / w, mode="same")
    out[:w // 2] = a[:w // 2]; out[-(w // 2):] = a[-(w // 2):]
    return out


def load(v):
    d = np.load(HERE / f"AXIS8_{v}.npz")
    x = d["x"]
    return x, boxfilt(d["DensityMEAN"] / RHO_J, x)


xb, rb = load("baseline")
exp = np.load(HERE.parent / "panda_m142den_axis_EPAPS.npz")
rb_on_x = np.interp(exp["x"], xb, rb)

plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})

fig, (axA, axB) = plt.subplots(
    2, 1, figsize=(7.6, 5.4), sharex=True,
    gridspec_kw={"height_ratios": [2.2, 1.15], "hspace": 0.08})
for ax in (axA, axB):
    ax.set_axisbelow(True); ax.grid(True, ls="--", lw=0.5, color="0.88")
    ax.tick_params(length=2.5)

mb = (xb > 0.05) & (xb <= XMAX)
axA.plot(xb[mb], rb[mb], color=CB, lw=2.0, zorder=3)
for v, lab, c, ls in VARIANTS:
    xv, rv = load(v)
    mv = (xv > 0.05) & (xv <= XMAX)
    axA.plot(xv[mv], rv[mv], color=c, lw=1.4, ls=ls, zorder=4)
    axB.plot(xv[mv], (rv - np.interp(xv, xb, rb))[mv], color=c, lw=1.5, ls=ls,
             zorder=4)
axA.plot(exp["x"], exp["rho"], "o", mfc="w", mec="k", mew=0.9, ms=4.2, zorder=5)

axA.set_ylim(0.3, 1.55)
axA.set_ylabel(r"$\langle\rho\rangle/\rho_j$")
handles = [
    Line2D([], [], marker="o", ls="", mfc="w", mec="k", ms=4.2,
           label=r"$M_j=1.42$ experiment"),
    Line2D([], [], color=CB, lw=2.0, label="Baseline"),
] + [Line2D([], [], color=c, lw=1.4, ls=ls, label=lab)
     for _, lab, c, ls in VARIANTS]
axA.legend(handles=handles, fontsize=8, frameon=False, ncol=3,
           loc="upper center", bbox_to_anchor=(0.5, 1.02),
           handlelength=2.0, columnspacing=1.3, handletextpad=0.5)

axB.axhline(0, color=CB, lw=1.5, ls=(0, (6, 3)), zorder=2)
axB.plot(exp["x"], exp["rho"] - rb_on_x, "o", mfc="w", mec="k", mew=0.9,
         ms=4.2, zorder=5)
axB.set_ylim(-0.72, 0.72)
axB.set_xlim(0, XMAX)
axB.set_xlabel(r"$x/D_e$")
axB.set_ylabel(r"$\Delta\langle\rho\rangle/\rho_j$" "\n(rel. baseline profile)",
               fontsize=8.5)
axB.text(0.15, 0.60, "(b)", fontsize=9, va="top")
axA.text(0.15, 0.42, "(a)", fontsize=9, va="top")

fig.subplots_adjust(left=0.10, right=0.975, top=0.94, bottom=0.095)
for ext in ("png", "pdf"):
    fig.savefig(HERE / f"axis_4profile.{ext}", dpi=250)
print("saved axis_4profile.png/pdf (density-only, two-panel)")
