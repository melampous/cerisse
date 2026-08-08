#!/usr/bin/env python3
"""Refinement-footprint budget for the L5-lip-only strategy.

Current production tagging (max_level 6) refines the whole shear corridor
to L5 out to 6 D_e and lets the density sensor build L6 on every shock cell
out to 4 D_e. For an x_fc study that is almost all wasted: L5+L6 carry 86%
of the cost and only the first cell matters.

This script budgets three lip-focused layouts against the measured
production counts and draws the footprints.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
from pathlib import Path

HERE = Path(__file__).resolve().parent
D = 25.4                                   # mm
R = 12.7
DX = {0: 1.5875, 1: 0.79375, 2: 0.396875, 3: 0.1984375,
      4: 0.09921875, 5: 0.049609375, 6: 0.0248046875}

# measured, m142_gpu_prod/phaseB.log
CUR = {0: 98304, 1: 98304, 2: 183296, 3: 569344,
       4: 1642496, 5: 2649088, 6: 1821696}


def cost(counts):
    return sum(n * 2**l for l, n in counts.items())


def cells(area_mm2, lev):
    return int(round(area_mm2 / DX[lev]**2))


# ------------------------------------------------- lip-band geometry
def band_area(zmax_D, w_in=2.5, w_out=2.5, spread=0.17):
    """Area of the band r in [R-w_in, R+w_out+spread*z], 0 <= z <= zmax."""
    z = zmax_D * D
    return (w_in + w_out) * z + 0.5 * spread * z**2


OPTS = {
    "L  lip band to 1.0 De": dict(zl=1.0, axis=0.0),
    "L+ lip band to 1.5 De + axis strip": dict(zl=1.5, axis=1.5),
    "F  whole first cell (r<=0.85D, x<=1.5D)": dict(zl=None, axis=None),
}

# common to all options
A_CORE = 2.0 * D * 1.0 * D                 # forced L4 block  z<=2D, r<=1D
L4_SENSOR = 180000                         # sensor L4 inside the tighter caps
L3_NEW, L2_NEW, L1_NEW, L0_NEW = 150000, 183296, 98304, 98304
PAD = 1.30                                 # blocking-factor / n_error_buf padding

print("measured production cost = %.1f M cell-updates per coarse step"
      % (cost(CUR) / 1e6))
print("   L5+L6 share = %.1f%%,  L0-L3 share = %.1f%%\n"
      % (100 * (CUR[5] * 32 + CUR[6] * 64) / cost(CUR),
         100 * sum(CUR[l] * 2**l for l in (0, 1, 2, 3)) / cost(CUR)))

print("%-42s %9s %9s %9s %10s %8s"
      % ("layout", "L5 cells", "L4 cells", "total", "M cell-upd", "speedup"))
rows = []
for nm, o in OPTS.items():
    if o["zl"] is None:
        a5 = 0.85 * D * 1.5 * D
    else:
        a5 = band_area(o["zl"])
        if o["axis"]:
            a5 += 2.0 * o["axis"] * D      # axis strip r <= 2 mm
    n5 = int(cells(a5, 5) * PAD)
    n4 = int((cells(A_CORE, 4) + L4_SENSOR) * PAD)
    new = {0: L0_NEW, 1: L1_NEW, 2: L2_NEW, 3: L3_NEW, 4: n4, 5: n5}
    c = cost(new)
    rows.append((nm, new, c))
    print("%-42s %9d %9d %9d %10.2f %7.1fx"
          % (nm, n5, n4, sum(new.values()), c / 1e6, cost(CUR) / c))

print("\nper-level breakdown, option L+ :")
nm, new, c = rows[1]
print("%6s %10s %8s %12s %7s" % ("level", "dx[um]", "cells", "cell-updates", "share"))
for l in sorted(new):
    print("%6d %10.1f %8d %12.2f M %6.1f%%"
          % (l, DX[l] * 1e3, new[l], new[l] * 2**l / 1e6,
             100 * new[l] * 2**l / c))

# throughput of the production run: 233.2 M cell-updates in 3.3 s on 4 GPUs
THR = cost(CUR) / 3.3
print("\nproduction throughput = %.1f M cell-updates/s (4 GPUs)" % (THR / 1e6))
for frac, lab in [(1.0, "if throughput held"), (0.4, "at 40% (launch-bound)"),
                  (0.25, "at 25% (pessimistic)")]:
    dt = c / (THR * frac)
    print("   %-24s %5.2f s/coarse step -> %5.1f min for 0.8 ms (1633 steps)"
          % (lab, dt, 1633 * dt / 60))

# ------------------------------------------------------------- figure
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in"})
fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.4, 4.3), sharex=True,
                             gridspec_kw={"hspace": 0.12})


def corridor(ax, zeta_max, lo0, lo1, hi0, hi1, **kw):
    z = np.linspace(0, zeta_max, 80)
    ax.fill_between(z, (lo0 + lo1 * z), (hi0 + hi1 * z), **kw)


# current
corridor(a1, 8.0, 0, 0, 1.5, 0, color="#c8d8ee", lw=0, label="L4 forced + caps")
corridor(a1, 6.0, 0.35, 0.04, 0.65, 0.12, color="#7fa8d8", lw=0,
         label="L5 shear corridor (to 6 $D_e$)")
corridor(a1, 1.0, 0.35, 0.04, 0.65, 0.12, color="#1f4e9c", lw=0,
         label="L6 lip corridor")
a1.add_patch(Rectangle((0, 0), 4.0, 1.25, fc="none", ec="#1f4e9c", lw=1.0,
                       ls="--"))
a1.text(4.1, 1.27, "L6 sensor cap ($z\\leq 4D$, $r\\leq 1.25D$)",
        fontsize=6.6, color="#1f4e9c", va="bottom")
a1.add_patch(Rectangle((0, 0), 8.0, 2.0, fc="none", ec="#4b7ec4", lw=1.0,
                       ls=":"))
a1.text(8.1, 1.72, "L5 sensor cap", fontsize=6.6, color="#4b7ec4",
        va="bottom")
a1.set_ylabel("$r/D_e$")
a1.set_ylim(0, 2.5)
a1.legend(fontsize=6.6, frameon=False, loc="upper right", ncol=3,
          handlelength=1.3, columnspacing=1.0)
a1.text(0.008, 0.90, "current  (max_level 6)   233.2 M cell-updates/step",
        transform=a1.transAxes, fontsize=8, va="top")

# proposed L+
a2.add_patch(Rectangle((0, 0), 2.0, 1.0, fc="#c8d8ee", ec="none"))
zz = np.linspace(0, 1.5, 80)
a2.fill_between(zz, (R - 2.5) / D, (R + 2.5 + 0.17 * zz * D) / D,
                color="#1f4e9c", lw=0)
a2.add_patch(Rectangle((0, 0), 1.5, 2.0 / D, fc="#1f4e9c", ec="none"))
a2.add_patch(Rectangle((0, 0), 3.0, 1.5, fc="none", ec="#4b7ec4", lw=1.0,
                       ls=":"))
a2.text(3.1, 1.52, "L4 sensor cap ($z\\leq 3D$, $r\\leq 1.5D$)",
        fontsize=6.6, color="#4b7ec4", va="bottom")
a2.axhline(0.5, color="k", lw=0.8)
for xv, lab in [(0.883, "$x_{fc}$"), (1.179, "$x_1$")]:
    a2.axvline(xv, color="#b3251e", lw=0.9, ls="--", ymax=0.78)
    a2.text(xv, 1.98, lab, color="#b3251e", fontsize=7.5, ha="center",
            va="top")
a2.set_ylabel("$r/D_e$"); a2.set_xlabel("$x/D_e$")
a2.set_xlim(0, 9.0); a2.set_ylim(0, 2.5)
a2.text(0.008, 0.90,
        "proposed L+  (max_level 5, L5 forced at the lip only)   "
        "%.1f M cell-updates/step" % (rows[1][2] / 1e6),
        transform=a2.transAxes, fontsize=8, va="top")
from matplotlib.patches import Patch
a2.legend(handles=[Patch(fc="#c8d8ee", label="L4 forced core ($z\\leq 2D$, $r\\leq 1D$)"),
                   Patch(fc="#1f4e9c", label="L5 lip band + axis strip")],
          fontsize=6.6, frameon=False, loc="upper right", ncol=2,
          handlelength=1.3)
fig.subplots_adjust(left=0.075, right=0.99, top=0.985, bottom=0.115)
fig.savefig(HERE / "plan_L5_lip.png", dpi=260)
fig.savefig(HERE / "plan_L5_lip.pdf")
print("\nwrote plan_L5_lip.png/pdf")
