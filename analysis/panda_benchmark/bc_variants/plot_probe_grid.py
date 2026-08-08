#!/usr/bin/env python3
"""Draw the realised AMR box layout of the L5-lip-only tagging, from the box
dump produced by dump_boxes.py, and print the per-level budget.

Usage: python3 plot_probe_grid.py boxes_np192.txt [out_tag]
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
FN = Path(sys.argv[1])
TAG = sys.argv[2] if len(sys.argv) > 2 else FN.stem
D = 0.0254
DX0 = 1.5875e-3

lev, boxes, time = None, {}, 0.0
for ln in FN.read_text().splitlines():
    if ln.startswith("# time"):
        time = float(ln.split()[2])
    elif ln.startswith("LEVEL"):
        lev = int(ln.split()[1]); boxes[lev] = []
    elif ln.strip() and lev is not None and not ln.startswith("#"):
        v = [float(x) for x in ln.split()]
        if len(v) == 4:
            boxes[lev].append(v)      # rlo rhi zlo zhi

print("plotfile time = %.4f ms" % (time * 1e3))
print("%6s %10s %8s %12s %14s %7s"
      % ("level", "dx[um]", "grids", "cells", "cell-updates", "share"))
tot = 0
rows = []
for l in sorted(boxes):
    dx = DX0 / 2**l
    n = sum(int(round((b[1] - b[0]) / dx)) * int(round((b[3] - b[2]) / dx))
            for b in boxes[l])
    rows.append((l, dx, len(boxes[l]), n, n * 2**l))
    tot += n * 2**l
for l, dx, ng, n, cu in rows:
    print("%6d %10.1f %8d %12d %12.2f M %6.1f%%"
          % (l, dx * 1e6, ng, n, cu / 1e6, 100 * cu / tot))
print("%6s %10s %8d %12d %12.2f M"
      % ("total", "", sum(r[2] for r in rows), sum(r[3] for r in rows),
         tot / 1e6))
print("\nproduction (max_level 6, old tagging) = 233.22 M  ->  speedup %.1fx"
      % (233.22e6 / tot))

# ------------------------------------------------------------------ figure
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in"})
fig, (ax, axz) = plt.subplots(1, 2, figsize=(7.6, 3.6),
                              gridspec_kw={"width_ratios": [1.55, 1.0]})
cmap = plt.cm.viridis(np.linspace(0.05, 0.92, len(rows)))
for a, (xlim, ylim) in zip((ax, axz), (((0, 8), (0, 2.5)),
                                       ((0, 2.0), (0, 0.85)))):
    for (l, dx, ng, n, cu), c in zip(rows, cmap):
        for b in boxes[l]:
            a.add_patch(Rectangle((b[2] / D, b[0] / D),
                                  (b[3] - b[2]) / D, (b[1] - b[0]) / D,
                                  fc="none", ec=c, lw=0.45))
    a.axhline(0.5, color="k", lw=0.9)
    a.set_xlim(*xlim); a.set_ylim(*ylim)
    a.set_xlabel("$x/D_e$"); a.set_ylabel("$r/D_e$")
for xv, lab in [(0.883, "$x_{fc}$"), (1.179, "$x_1$")]:
    for a in (ax, axz):
        a.axvline(xv, color="#b3251e", lw=0.9, ls="--", zorder=5)
    axz.text(xv, 0.82, lab, color="#b3251e", fontsize=8, ha="center",
             va="top")
ax.legend(handles=[Line2D([], [], color=c, lw=1.4,
                          label="L%d  %.1f $\\mu$m  (%d grids)"
                                % (l, dx * 1e6, ng))
                   for (l, dx, ng, n, cu), c in zip(rows, cmap)],
          fontsize=6.4, frameon=False, loc="upper right", labelspacing=0.25,
          handlelength=1.4)
ax.text(0.015, 0.965, "$t$ = %.3f ms" % (time * 1e3), transform=ax.transAxes,
        va="top", fontsize=8)
axz.text(0.03, 0.05, "zoom: first shock cell", transform=axz.transAxes,
         fontsize=7.5)
fig.tight_layout(w_pad=1.2)
fig.savefig(HERE / ("probe_grid_%s.png" % TAG), dpi=260)
fig.savefig(HERE / ("probe_grid_%s.pdf" % TAG))
print("\nwrote probe_grid_%s.png/pdf" % TAG)
