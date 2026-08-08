#!/usr/bin/env python3
"""The actual AMR box layout of m142 and m180, read out of the plotfiles.

Boxes come from Level_N/Cell_H in index space, converted here to physical
coordinates. Only boxes that intersect the z = 0 mid-plane are drawn, so the
panels are a true slice through the hierarchy rather than a projection.

Both cases use the same base mesh: 192 x 128 x 128 over 24 x 16 x 16 D_e,
i.e. 3.175 mm at level 0, so level 5 is 99.2 um in both.
"""
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
D = 0.0254
N0 = (192, 128, 128)
COL = {0: "#c8d8ee", 1: "#9dbde0", 2: "#6f9dd0", 3: "#2ca02c",
       4: "#ff6a00", 5: "#b3251e", 6: "#7b3294"}


def read(fn):
    lo = hi = None
    boxes = []
    for ln in Path(fn).read_text().splitlines():
        if ln.startswith("# domain"):
            p = ln.split()
            lo = [float(v) for v in p[3:6]]
            hi = [float(v) for v in p[7:10]]
        elif ln and not ln.startswith("#"):
            boxes.append([int(v) for v in ln.split()])
    return np.array(boxes), lo, hi


def panel(ax, fn, title):
    B, lo, hi = read(fn)
    dx0 = (hi[0]-lo[0])/N0[0]
    dy0 = (hi[1]-lo[1])/N0[1]
    dz0 = (hi[2]-lo[2])/N0[2]
    counts = {}
    for lev, il, jl, kl, ih, jh, kh in B:
        r = 2**lev
        kmid = (N0[2]*r)//2                      # cell index straddling z = 0
        if not (kl <= kmid <= kh):
            continue
        x0 = lo[0] + il*dx0/r
        x1 = lo[0] + (ih+1)*dx0/r
        y0 = lo[1] + jl*dy0/r
        y1 = lo[1] + (jh+1)*dy0/r
        ax.add_patch(Rectangle((x0/D, y0/D), (x1-x0)/D, (y1-y0)/D,
                               fc="none", ec=COL.get(lev, "k"), lw=0.55))
        counts[lev] = counts.get(lev, 0) + 1
    ax.axhline(0.5, color="k", lw=0.7, alpha=0.5)
    ax.axhline(-0.5, color="k", lw=0.7, alpha=0.5)
    ax.set_xlim(0, 6.0); ax.set_ylim(-2.2, 2.2)
    ax.set_xlabel("$x/D_e$"); ax.set_ylabel("$y/D_e$")
    ax.set_title(title, fontsize=9)
    return counts


CELLS = {"m142": {0: 3_145_728, 1: 4_194_304, 2: 4_194_304, 3: 7_372_800,
                  4: 19_398_656, 5: 12_582_912},
         "m180": None}

plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in"})
files = [("sweep_out/BOX_m142.txt", "m142  $NPR_0$=3.27,  fixed grid file for L5"),
         ("sweep_out/BOX_m180.txt", "m180  $NPR_0$=5.75,  forced zones + sensor")]
have = [(f, t) for f, t in files if (HERE/f).exists()]
fig, axs = plt.subplots(len(have), 1, figsize=(8.6, 3.1*len(have)),
                        sharex=True, sharey=True,
                        gridspec_kw={"hspace": 0.22})
if len(have) == 1:
    axs = [axs]
for ax, (f, t) in zip(axs, have):
    c = panel(ax, HERE/f, t)
    print("%-28s %s" % (t.split()[0], {k: c[k] for k in sorted(c)}))
axs[-1].set_xlabel("$x/D_e$")
lv = sorted(set().union(*[set(read(HERE/f)[0][:, 0]) for f, _ in have]))
axs[0].legend(handles=[Line2D([], [], color=COL.get(int(l), "k"), lw=1.4,
                              label="L%d  %.1f $\\mu$m"
                                    % (l, 3175.0/2**int(l)))
                       for l in lv],
              fontsize=6.6, frameon=False, ncol=len(lv),
              loc="upper right", handlelength=1.3, columnspacing=1.0)
fig.tight_layout()
fig.savefig(HERE / "box_layout_m142_m180.png", dpi=260)
fig.savefig(HERE / "box_layout_m142_m180.pdf")
print("wrote box_layout_m142_m180.png/pdf")
