#!/usr/bin/env python3
"""Design figure for the 2D inlet-profile sweep.

(a) every prescribed profile in the proposed matrix, drawn on the real
    lip-corridor mesh (dx = 24.8 um, level 6) so the delivered resolution
    is visible;
(b) the discrimination test: for the three form-matched wall cases, the
    first-cell closure x_fc predicted by three competing laws calibrated
    on the four existing configurations. Where the laws disagree by much
    more than the extraction resolution, one run decides between them.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path

from design_2d_sweep import scalars, fvel, R, D, DX_LIP, EXIST

HERE = Path(__file__).resolve().parent
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})

ARM_A = [0.0, 100.0, 170.0, 254.0, 340.0, 425.0]
ANCHOR = 170.0
WALL = [("B1  $\\delta^*$-matched", 679.5), ("B2  $\\delta_\\omega$-matched", 340.0),
        ("B3  $\\theta$-matched", 270.4), ("C   3D twin", 200.0)]

# laws calibrated on the four existing configurations
y = np.array([e[3] for e in EXIST])
S = [scalars(e[1] * 1e-6, e[2]) for e in EXIST]
LAWS = []
for key, nm, col in [("dstar_c", r"$\delta^*$ law (aperture)", "#b3251e"),
                     ("dw", r"$\delta_\omega$ law", "#1f4e9c"),
                     ("theta_c", r"$\theta$ law", "#2ca02c")]:
    xv = np.array([s[key] for s in S]) / D
    m, c = np.polyfit(xv, y, 1)
    LAWS.append((nm, key, col, m, c))

fig = plt.figure(figsize=(7.4, 3.35))
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.15], wspace=0.30)

# ---------------------------------------------------------------- panel (a)
ax = fig.add_subplot(gs[0, 0])
RLO, RHI = 0.372, 0.5105
rf = np.arange(0.0, R + DX_LIP, DX_LIP)
for k, v in enumerate(rf / D):
    if RLO <= v <= RHI:
        ax.axhline(v, color="0.87", lw=0.35, zorder=0)
        if k % 2 == 0:
            ax.axhspan(v, min(v + DX_LIP / D, RHI), color="0.965", zorder=-1)
r = np.linspace(0.30 * D, R, 6000)
cm = plt.cm.plasma(np.linspace(0.05, 0.80, len(ARM_A)))
for a, c in zip(ARM_A, cm):
    f = fvel(r, a * 1e-6, "shifted")
    ax.plot(np.r_[f, 0.0], np.r_[r / D, R / D], color=c, lw=1.4, zorder=4,
            label="A  $a$=%g" % a)
for (lab, aw), ls in zip(WALL, ["--", "--", "--", ":"]):
    f = fvel(r, aw * 1e-6, "wall")
    ax.plot(np.r_[f, 0.0], np.r_[r / D, R / D], color="0.25", lw=1.0, ls=ls,
            zorder=3)
ax.axhline(R / D, color="k", lw=1.1, zorder=5)
ax.text(0.50, 0.972, "physical lip  $r=R$", transform=ax.transAxes,
        fontsize=7, va="bottom", ha="center")
ax.set_xlim(-0.03, 1.06); ax.set_ylim(RLO, RHI)
ax.set_xlabel("$u/U_e$"); ax.set_ylabel("$r/D_e$")
leg = ax.legend(fontsize=6.3, frameon=False, loc="lower left",
                handlelength=1.2, labelspacing=0.25, borderpad=0.2,
                title="shifted tanh (Arm A)")
leg.get_title().set_fontsize(6.3)
ax.plot([], [], color="0.25", lw=1.0, ls="--")
ax.text(0.965, 0.035, "(a)", transform=ax.transAxes, ha="right", fontsize=9)
ax.text(0.33, 0.135, "dashed: wall tanh (Arms B, C)", transform=ax.transAxes,
        fontsize=6.3, color="0.25")
ax.text(0.33, 0.030,
        "background = real mesh,\n$\\Delta x$=24.8 $\\mu$m (level 6)",
        transform=ax.transAxes, fontsize=6.3, color="0.35", va="bottom")

# ---------------------------------------------------------------- panel (b)
ax2 = fig.add_subplot(gs[0, 1])
sa = scalars(ANCHOR * 1e-6, "shifted")
xa = LAWS[0][4] + LAWS[0][3] * sa["dstar_c"] / D
ax2.axhline(xa, color="0.35", lw=1.0, ls="-", zorder=1)
ax2.text(len(WALL) - 0.45, xa - 0.005,
         "anchor: shifted $a$=170 $\\mu$m (Arm A)", fontsize=6.6,
         color="0.25", ha="right", va="top")
pos = np.arange(len(WALL))
for j, (lab, aw) in enumerate(WALL):
    s = scalars(aw * 1e-6, "wall")
    for (nm, key, col, m, c) in LAWS:
        ax2.plot(j, c + m * s[key] / D, "o", color=col, ms=7, mec="w",
                 mew=0.8, zorder=4)
    vals = [c + m * s[key] / D for (nm, key, col, m, c) in LAWS]
    ax2.plot([j, j], [min(vals), max(vals)], color="0.6", lw=0.8, zorder=2)
    ax2.text(j, max(vals) + 0.004, "%.3f" % (max(vals) - min(vals)),
             ha="center", va="bottom", fontsize=6.6, color="0.3")
ax2.set_xticks(pos)
ax2.set_xticklabels([w[0] + "\n$a$=%g $\\mu$m" % w[1] for w in WALL],
                    fontsize=6.5)
ax2.set_xlim(-0.5, len(WALL) - 0.5)
ax2.set_ylim(0.836, 0.952)
ax2.set_ylabel("$x_{fc}/D_e$ predicted by each competing law")
ax2.grid(True, axis="y", ls="--", lw=0.5, color="0.88")
ax2.set_axisbelow(True)
h = [Line2D([], [], marker="o", ls="", color=c, ms=6, label=n)
     for (n, k, c, m, cc) in LAWS]
ax2.legend(handles=h, fontsize=6.4, frameon=False, loc="lower right",
           handlelength=1.2, labelspacing=0.25,
           title="spread between laws printed above each case")
ax2.get_legend().get_title().set_fontsize(6.4)
ax2.text(0.025, 0.03, "(b)", transform=ax2.transAxes, fontsize=9)
ax2.text(0.025, 0.965,
         "measurable: extraction resolution $\\approx$0.001 $D_e$",
         transform=ax2.transAxes, va="top", fontsize=6.4, color="0.35")

fig.subplots_adjust(left=0.085, right=0.985, top=0.975, bottom=0.175)
fig.savefig(HERE / "sweep2d_design.png", dpi=260)
fig.savefig(HERE / "sweep2d_design.pdf")
print("wrote sweep2d_design.png/pdf")
for j, (lab, aw) in enumerate(WALL):
    s = scalars(aw * 1e-6, "wall")
    vals = {nm: c + m * s[key] / D for (nm, key, col, m, c) in LAWS}
    print("%-26s a=%6.1f  " % (lab, aw)
          + "  ".join("%s=%.4f" % (k.split()[0], v) for k, v in vals.items())
          + "   spread=%.4f" % (max(vals.values()) - min(vals.values())))
