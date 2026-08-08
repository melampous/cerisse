#!/usr/bin/env python3
"""Two checks on the first shock-cell closure position.

Panel (a): the extracted x_fc marked on each case's own centreline mean
density, so its placement can be judged by eye rather than trusted.
Panel (b): the new reduced-grid case A3_a254 against the earlier full-grid
2D production run m142_2d_L6stat, which uses an identical inlet profile
(shifted tanh, a = 254 um, s = 3). Any difference between them is the cost of
the grid reduction, not of the physics.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
D, RHO_J = 0.0254, 1.6413
T0_NEW = 4.0e-4
T0_OLD = 0.01574728874

CASES = [("A0_tophat", "#404040"), ("A1_a100", "#7b3294"),
         ("A2_a170", "#b3251e"), ("A3_a254", "#e08214"),
         ("S4_a170_s4", "#c2a5cf"), ("S5_a170_s5", "#8073ac"),
         ("B1_wall680", "#1f4e9c"), ("C1_wall200", "#2ca02c")]


def bracket(x, r):
    e = np.where((x >= 0.35) & (x <= 1.0))[0]
    i0 = e[np.nanargmin(r[e])]
    c = np.where((np.arange(x.size) > i0) & (x <= 1.7))[0]
    return i0, c[np.nanargmax(r[c])]


def xfc_E3(x, r):
    i0, i1 = bracket(x, r)
    h = 0.5 * (r[i0] + r[i1])
    k = np.where(r[i0:i1+1] >= h)[0]
    if not len(k):
        return np.nan, i0, i1, h
    k = k[0] + i0
    v = x[i0] if k == i0 else x[k-1] + (h-r[k-1])*(x[k]-x[k-1])/(r[k]-r[k-1])
    return v, i0, i1, h


def final(fn):
    d = np.load(fn)
    ks = sorted([k[:-2] for k in d.files if k.endswith("_t")],
                key=lambda s: int("".join(c for c in s if c.isdigit())))
    dx = float(d[ks[0]+"_t"][1])
    M = d[ks[-1]+"_rho"] / RHO_J
    return (np.arange(M.size)+0.5)*dx/D, M


plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.8, 3.5))

# ---------------------------------------------------------------- panel (a)
for nm, col in CASES:
    f = HERE / "sweep_out" / ("XFC_%s.npz" % nm)
    if not f.exists():
        continue
    x, M = final(f)
    v, i0, i1, h = xfc_E3(x, M)
    m = (x > 0.3) & (x < 1.8)
    a1.plot(x[m], M[m], color=col, lw=1.1)
    a1.plot(v, h, "o", color=col, ms=6, mec="w", mew=0.8, zorder=6)
    a1.plot(x[i0], M[i0], "v", color=col, ms=4, alpha=0.6, zorder=5)
    a1.plot(x[i1], M[i1], "^", color=col, ms=4, alpha=0.6, zorder=5)
a1.set_xlabel("$x / D_e$")
a1.set_ylabel(r"$\langle \rho \rangle / \rho_j$   (centreline)")
a1.set_xlim(0.3, 1.8)
a1.grid(True, ls="--", lw=0.5, color="0.9"); a1.set_axisbelow(True)
a1.legend(handles=[Line2D([], [], marker="v", ls="", color="0.4", ms=5,
                          label="expansion minimum"),
                   Line2D([], [], marker="^", ls="", color="0.4", ms=5,
                          label="compression maximum"),
                   Line2D([], [], marker="o", ls="", color="0.4", ms=6,
                          label="$x_{fc}$ = midpoint crossing")],
           fontsize=6.6, frameon=False, loc="upper left", handlelength=1.0)
a1.text(0.965, 0.05, "(a)", transform=a1.transAxes, ha="right", fontsize=9)

# ---------------------------------------------------------------- panel (b)
xn, Mn = final(HERE / "sweep_out" / "XFC_A3_a254.npz")
vn = xfc_E3(xn, Mn)[0]
a2.plot(xn, Mn, color="#e08214", lw=1.5,
        label="A3_a254  new, lip at level 5")
a2.axvline(vn, color="#e08214", lw=1.0, ls="--")

old = HERE / "XFC_window_L6.npz"
if old.exists():
    xo, Mo = final(old)
    vo = xfc_E3(xo, Mo)[0]
    a2.plot(xo, Mo, color="#1f4e9c", lw=1.2, ls="-", alpha=0.85,
            label="m142_2d_L6stat  earlier, lip at level 6")
    a2.axvline(vo, color="#1f4e9c", lw=1.0, ls=":")
    print("A3_a254 (new, L5 lip, 2.0 ms)        x_fc = %.4f" % vn)
    print("m142_2d_L6stat (earlier, L6, 1.5 ms) x_fc = %.4f" % vo)
    print("difference                                 %+.4f D_e" % (vn - vo))
    a2.text(0.03, 0.06,
            "$x_{fc}$: new %.4f, earlier %.4f\ndifference %+.4f $D_e$"
            % (vn, vo, vn - vo), transform=a2.transAxes, fontsize=7.2)
a2.set_xlabel("$x / D_e$")
a2.set_ylabel(r"$\langle \rho \rangle / \rho_j$   (centreline)")
a2.set_xlim(0.0, 4.0)
a2.grid(True, ls="--", lw=0.5, color="0.9"); a2.set_axisbelow(True)
a2.legend(fontsize=6.8, frameon=False, loc="upper right", handlelength=1.6)
a2.text(0.965, 0.30, "(b)", transform=a2.transAxes, ha="right", fontsize=9)
fig.tight_layout(w_pad=1.3)
fig.savefig(HERE / "xfc_placement.png", dpi=260)
fig.savefig(HERE / "xfc_placement.pdf")
print("\nwrote xfc_placement.png/pdf")
