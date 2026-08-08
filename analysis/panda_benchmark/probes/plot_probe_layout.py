#!/usr/bin/env python3
"""Where the 210 pressure probes sit.

The probe list in the input deck gives names only; the positions are in the
data-log header, as the level-1 index box each probe averages over. A box
(i0,j0,k0)-(i1,j1,k1) covers cells i0..i1 inclusive, so its centre in index
space is (i0+i1+1)/2 and

    X = prob_lo + (i0+i1+1)/2 * dx_1,    dx_1 = 3175/2 um = 1587.5 um

with prob_lo = (0, -0.2032, -0.2032) m. Each probe therefore reports the mean
pressure over a 2x2x2 block of level-1 cells, 3.175 mm on a side, which is
0.125 D_e - a spatial average, not a point value. That matters for the highest
frequencies: a wave shorter than about twice the block is attenuated by the
averaging, which puts a soft ceiling near 40 kHz well below the sampling
Nyquist of 85-90 kHz.

The layout is built for azimuthal decomposition. Twelve of the fourteen groups
are 16-point rings at 22.5 degree spacing, which resolves modes m = 0 .. +-7
before aliasing; the remaining groups are a streamwise line on the lip line,
three points on the axis, and two four-point rings that are a coarse subset of
the sixteen.

The rings are not exactly circular: the probes are snapped to the Cartesian
grid, so a nominal r = 0.5 D_e ring has radii between 0.44 and 0.53. The
decomposition below therefore uses the actual angle of each probe rather than
an assumed uniform spacing.
"""
from pathlib import Path
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
D = 0.0254
PLO = np.array([0.0, -0.2032, -0.2032])
DX1 = 3175e-6/2

names = re.search(r"^cns\.time_probes\s*=\s*(.*)$",
                  (HERE/"campaign_inputs").read_text(), re.M).group(1).split()
hdr = (HERE/"s3d_skewjst_probes.log").open().readline()
box = np.array(re.findall(
    r'"pressure\(\((-?\d+),(-?\d+),(-?\d+)\)\s*\((-?\d+),(-?\d+),(-?\d+)\)\)"',
    hdr), int)
XYZ = PLO + 0.5*(box[:, :3] + box[:, 3:] + 1)*DX1
P = XYZ/D
R = np.hypot(P[:, 1], P[:, 2])
PHI = np.degrees(np.arctan2(P[:, 2], P[:, 1]))


def grp(n):
    """s1_07 -> s1 ; r1a -> r1 ; q3 -> q ; a2 -> a"""
    if "_" in n:
        return n.split("_")[0]
    if n[-1].isalpha() and any(c.isdigit() for c in n):
        return n.rstrip("abcdefgh")          # r1a..r1d -> r1
    return n.rstrip("0123456789")            # q3 -> q


G = {}
for i, n in enumerate(names):
    G.setdefault(grp(n), []).append(i)

STY = {"s05": ("#0072B2", "o", 34), "s1": ("#0072B2", "o", 34),
       "s2": ("#0072B2", "o", 34), "s3": ("#0072B2", "o", 34),
       "s4": ("#0072B2", "o", 34), "s6": ("#0072B2", "o", 34),
       "o1": ("#D55E00", "s", 30), "o2": ("#D55E00", "s", 30),
       "o3": ("#D55E00", "s", 30), "o4": ("#D55E00", "s", 30),
       "n15": ("#009E73", "^", 34), "n30": ("#CC79A7", "v", 34),
       "q": ("#333333", "D", 26), "a": ("#E69F00", "*", 90),
       "m": ("#56B4E9", "P", 40), "r1": ("0.45", "x", 34),
       "r2": ("0.45", "x", 34)}
LAB = {"s1": "lip-line rings  $r\\simeq0.5$  ($x$ = 0.5, 1, 2, 3, 4, 6)",
       "o1": "outer rings  $r\\simeq1.0$  ($x$ = 1, 2, 3, 4)",
       "n15": "ring  $r$ = 1.5, $x$ = 2", "n30": "ring  $r$ = 3.0, $x$ = 2",
       "q": "lip line, $\\varphi$ = 0", "a": "axis  ($x$ = 1, 2, 3)",
       "m": "$r$ = 0.94  ($x$ = 1, 2)", "r1": "4-point rings  $r$ = 0.5"}

print("%-6s %4s %14s %14s %10s  %s"
      % ("group", "n", "x [D_e]", "r [D_e]", "d_phi", "kind"))
for k in ["q", "a", "m", "r1", "r2", "s05", "s1", "s2", "s3", "s4", "s6",
          "o1", "o2", "o3", "o4", "n15", "n30"]:
    i = G[k]
    ring = len(i) >= 8 and np.ptp(R[i]) < 0.12 and np.ptp(P[i, 0]) < 0.06
    print("%-6s %4d %14s %14s %10s  %s"
          % (k, len(i), "%.2f" % P[i, 0].mean() if np.ptp(P[i, 0]) < 0.02
             else "%.2f-%.2f" % (P[i, 0].min(), P[i, 0].max()),
             "%.3f" % R[i].mean() if np.ptp(R[i]) < 0.02
             else "%.2f-%.2f" % (R[i].min(), R[i].max()),
             "%.1f deg" % (360.0/len(i)) if ring else "-",
             "azimuthal ring, modes to m=%d" % (len(i)//2 - 1) if ring
             else "line/points"))

plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True,
                     "legend.frameon": False})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.0, 4.9),
                             gridspec_kw=dict(width_ratios=[1.75, 1.0]))
fig.subplots_adjust(left=0.058, right=0.985, top=0.925, bottom=0.125,
                    wspace=0.20)

# (a) meridional view: every probe at its own (x, r), sign of y kept so the
# two halves of each ring are visible
seen = set()
for k, i in G.items():
    c, m, s = STY[k]
    lab = LAB.get(k) if k in LAB and LAB[k] not in seen else None
    if lab:
        seen.add(LAB[k])
    a1.scatter(P[i, 0], np.sign(P[i, 1]+1e-9)*R[i], s=s, marker=m, c=c,
               edgecolors="k", linewidths=0.45, zorder=5, label=lab)
a1.axhline(0.5, color="0.75", lw=0.8, ls=":", zorder=2)
a1.axhline(-0.5, color="0.75", lw=0.8, ls=":", zorder=2)
a1.plot([0, 0], [-0.5, 0.5], "-", color="k", lw=2.5, zorder=6)
a1.text(0.06, 0.0, "nozzle", fontsize=8.4, va="center")
a1.set_xlim(-0.25, 6.6)
a1.set_ylim(-3.4, 3.4)
a1.set_xlabel("$x/D_e$")
a1.set_ylabel(r"$\pm r/D_e$")
a1.set_title("(a)  meridional view, all 210 probes", fontsize=10, pad=5)
a1.legend(loc="upper right", fontsize=8.0, handlelength=1.2, ncol=2,
          columnspacing=1.0, scatterpoints=1)
a1.grid(True, ls="--", lw=0.5, color="0.93")
a1.set_axisbelow(True)

# (b) cross-section: the rings seen down the axis
for k in ["s1", "o2", "n15", "n30"]:
    i = G[k]
    c, m, s = STY[k]
    a2.scatter(P[i, 1], P[i, 2], s=s+8, marker=m, c=c, edgecolors="k",
               linewidths=0.45, zorder=5,
               label="%s  ($r\\simeq%.2f$)" % (k, R[i].mean()))
    th = np.linspace(0, 2*np.pi, 200)
    a2.plot(R[i].mean()*np.cos(th), R[i].mean()*np.sin(th), "-", color=c,
            lw=0.7, alpha=0.45, zorder=3)
th = np.linspace(0, 2*np.pi, 200)
a2.plot(0.5*np.cos(th), 0.5*np.sin(th), "-", color="k", lw=1.6, zorder=6)
a2.text(0.0, 0.0, "nozzle", fontsize=8.2, ha="center", va="center")
a2.set_aspect("equal")
a2.set_xlim(-3.4, 3.4)
a2.set_ylim(-3.4, 3.4)
a2.set_xlabel("$y/D_e$")
a2.set_ylabel("$z/D_e$")
a2.set_title(r"(b)  the rings, viewed along the axis", fontsize=10, pad=5)
a2.legend(loc="upper right", fontsize=7.8, handlelength=1.1, scatterpoints=1)
a2.grid(True, ls="--", lw=0.5, color="0.93")
a2.set_axisbelow(True)

np.savez(HERE/"probe_coords.npz", names=np.array(names), xyz=XYZ, box=box,
         r=R, phi=PHI)
fig.savefig(HERE/"probe_layout.png", dpi=250)
fig.savefig(HERE/"probe_layout.pdf")
print("\nwrote probe_layout.png/pdf and probe_coords.npz")
