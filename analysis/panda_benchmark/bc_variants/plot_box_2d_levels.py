#!/usr/bin/env python3
"""2D AMR box layout at three maximum refinement levels, read out of the
plotfiles and checkpoints themselves (Level_N/Cell_H for plotfiles,
Level_N/SD_0_New_MF_H for checkpoints).

L4 and L5 are the same top-hat case run with max_level = 4 and 5 under the
lip-focused tagging. L6 is the earlier production run m142_2d_L6stat, which
used the wide shear-corridor tagging and the a = 254 um profile, so its
footprint reflects a different refinement policy as well as a deeper
hierarchy - that is labelled on the panel rather than glossed over.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
D = 0.0254
N0 = (192, 512)                      # base cells, (r, z)
COL = {0: "#c8d8ee", 1: "#9dbde0", 2: "#6f9dd0", 3: "#2ca02c",
       4: "#ff6a00", 5: "#b3251e", 6: "#7b3294"}
PAN = [("BOX2D_L4.txt", "max_level 4   finest 99.2 $\\mu$m   top hat, lip-focused tagging"),
       ("BOX2D_L5.txt", "max_level 5   finest 49.6 $\\mu$m   top hat, lip-focused tagging"),
       ("BOX2D_L6.txt", "max_level 6   finest 24.8 $\\mu$m   m142_2d_L6stat, $a$=254 $\\mu$m, wide-corridor tagging")]


def read(fn):
    lo = hi = None; B = []
    for ln in Path(fn).read_text().splitlines():
        if ln.startswith("# domain"):
            p = ln.split(); lo = [float(p[3]), float(p[4])]
            hi = [float(p[6]), float(p[7])]
        elif ln and not ln.startswith("#"):
            B.append([int(v) for v in ln.split()])
    return np.array(B), lo, hi


have = [(f, t) for f, t in PAN if (HERE/"sweep_out"/f).exists()]
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in"})
fig, axs = plt.subplots(len(have), 1, figsize=(9.0, 2.5*len(have)),
                        sharex=True, gridspec_kw={"hspace": 0.30})
if len(have) == 1:
    axs = [axs]
alllv = set()
for ax, (f, t) in zip(axs, have):
    B, lo, hi = read(HERE/"sweep_out"/f)
    dr0 = (hi[0]-lo[0])/N0[0]
    dz0 = (hi[1]-lo[1])/N0[1]
    n = {}
    for lev, il, jl, ih, jh in B:
        r = 2**lev
        r0 = lo[0] + il*dr0/r; r1 = lo[0] + (ih+1)*dr0/r
        z0 = lo[1] + jl*dz0/r; z1 = lo[1] + (jh+1)*dz0/r
        if z0/D > 4.0 or r0/D > 1.6:
            continue
        ax.add_patch(Rectangle((z0/D, r0/D), (z1-z0)/D, (r1-r0)/D,
                               fc="none", ec=COL.get(lev, "k"), lw=0.5))
        n[lev] = n.get(lev, 0) + 1
        alllv.add(lev)
    ax.axhline(0.5, color="k", lw=0.8, alpha=0.6)
    ax.set_xlim(0, 4.0); ax.set_ylim(0, 1.6)
    ax.set_ylabel("$r/D_e$")
    ax.set_title(t, fontsize=8.5)
    ax.text(0.985, 0.06, "boxes: " + "  ".join("L%d:%d" % (k, n[k])
                                               for k in sorted(n)),
            transform=ax.transAxes, ha="right", fontsize=7, color="0.3")
axs[-1].set_xlabel("$x/D_e$")
axs[0].legend(handles=[Line2D([], [], color=COL[l], lw=1.5,
                              label="L%d  %.1f $\\mu$m" % (l, 1587.5/2**l))
                       for l in sorted(alllv)],
              fontsize=6.4, frameon=False, ncol=len(alllv),
              loc="upper right", handlelength=1.2, columnspacing=0.9)
fig.subplots_adjust(left=0.07, right=0.99, top=0.94, bottom=0.09)
fig.savefig(HERE/"box_2d_L4_L5_L6.png", dpi=260)
fig.savefig(HERE/"box_2d_L4_L5_L6.pdf")
print("wrote box_2d_L4_L5_L6.png/pdf")
