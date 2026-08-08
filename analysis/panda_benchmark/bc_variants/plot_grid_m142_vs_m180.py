#!/usr/bin/env python3
"""Why the m180 level-5 grid is twelve times the m142 one.

Both cases share the same base mesh (D_e/8 over a 24 x 16 x 16 D_e box) and
the same refinement ratio, so the level spacings are identical. What differs
is how level 5 is created:

  m142  a hand-drawn fixed grid file (fixed_grids_v5_L5.dat) with
        regrid_int = 100000, so level 5 is exactly the patch that was drawn;
  m180  geometric forced zones plus the density sensor. The forced zone
        L5-B was a solid cylinder r <= 0.40 D over 1.1 <= x/D <= 2.2, sized
        with a two-dimensional area estimate. In three dimensions that block
        alone is about 120 M cells - more than the whole level-5 budget.

Counts are measured, taken from the "Advanced N cells" lines of each run.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
DX = {3: 198.4, 4: 99.2, 5: 49.6}          # um, but 3D base differs - see below
DX3D = {3: 396.9, 4: 198.4, 5: 99.2}       # um, 3D: base D_e/8 = 3175 um

M142 = {3: 7_372_800, 4: 19_398_656, 5: 12_582_912}
M180 = {3: 12_976_128, 4: 52_953_088, 5: 159_055_872}
STEP = {"m142": 4.68, "m180": 52.0}        # measured s per coarse step

plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "ytick.right": True})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.6, 3.4))

lv = [3, 4, 5]
w = 0.36
xp = np.arange(len(lv))
c142, c180 = "#1f4e9c", "#b3251e"

# ------------------------------------------------------------- panel (a)
a1.bar(xp - w/2, [M142[l]/1e6 for l in lv], w, color=c142,
       label="m142  fixed grid file")
a1.bar(xp + w/2, [M180[l]/1e6 for l in lv], w, color=c180,
       label="m180  forced zones + sensor")
for i, l in enumerate(lv):
    a1.text(xp[i]-w/2, M142[l]/1e6, "%.1f" % (M142[l]/1e6), ha="center",
            va="bottom", fontsize=7, color=c142)
    a1.text(xp[i]+w/2, M180[l]/1e6, "%.1f" % (M180[l]/1e6), ha="center",
            va="bottom", fontsize=7, color=c180)
    a1.text(xp[i]+w/2, M180[l]/1e6*1.18, "%.1f$\\times$" % (M180[l]/M142[l]),
            ha="center", fontsize=7.5, color="0.25")
a1.set_xticks(xp)
a1.set_xticklabels(["L3\n396.9 $\\mu$m", "L4\n198.4 $\\mu$m", "L5\n99.2 $\\mu$m"])
a1.set_ylabel("cells  [millions]")
a1.set_ylim(0, 210)
a1.grid(True, axis="y", ls="--", lw=0.5, color="0.9"); a1.set_axisbelow(True)
a1.legend(fontsize=7, frameon=False, loc="upper left")
a1.text(0.97, 0.05, "(a) cells per level", transform=a1.transAxes,
        ha="right", fontsize=8)

# ------------------------------------------------------------- panel (b)
u142 = [M142[l]*2**l/1e9 for l in lv]
u180 = [M180[l]*2**l/1e9 for l in lv]
a2.bar(xp - w/2, u142, w, color=c142)
a2.bar(xp + w/2, u180, w, color=c180)
for i, l in enumerate(lv):
    a2.text(xp[i]-w/2, u142[i], "%.2f" % u142[i], ha="center", va="bottom",
            fontsize=7, color=c142)
    a2.text(xp[i]+w/2, u180[i], "%.2f" % u180[i], ha="center", va="bottom",
            fontsize=7, color=c180)
a2.set_xticks(xp)
a2.set_xticklabels(["L3  $\\times 8$", "L4  $\\times 16$", "L5  $\\times 32$"])
a2.set_ylabel("cell-updates per coarse step  [billions]")
a2.set_ylim(0, 6.2)
a2.grid(True, axis="y", ls="--", lw=0.5, color="0.9"); a2.set_axisbelow(True)
a2.text(0.03, 0.94,
        "total   m142  %.2f B  ->  %.2f s/step\n"
        "        m180  %.2f B  ->  %.0f s/step"
        % (sum(u142), STEP["m142"], sum(u180), STEP["m180"]),
        transform=a2.transAxes, va="top", fontsize=7.5, family="monospace")
a2.text(0.97, 0.05, "(b) weighted by subcycling", transform=a2.transAxes,
        ha="right", fontsize=8)
fig.tight_layout(w_pad=1.3)
fig.savefig(HERE / "grid_m142_vs_m180.png", dpi=260)
fig.savefig(HERE / "grid_m142_vs_m180.pdf")
print("wrote grid_m142_vs_m180.png/pdf")
print("cell-updates  m142 = %.3f B   m180 = %.3f B   ratio %.1f"
      % (sum(u142), sum(u180), sum(u180)/sum(u142)))
print("step time     m142 = %.2f s   m180 = %.0f s   ratio %.1f"
      % (STEP["m142"], STEP["m180"], STEP["m180"]/STEP["m142"]))
print("\nto match the m142 level-5 budget, m180 L5 must shrink %.0fx"
      % (M180[5]/M142[5]))
