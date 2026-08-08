#!/usr/bin/env python3
"""NPR 20 top hat: centreline mean density at Level 3 against Level 4.

Same binary (md5 f689f3485f92), same prob.h, same tree, same two-phase
protocol; the only difference between the two runs is amr.max_level, so the
comparison isolates the grid.

If the Level 3 run is still going, its statistics are taken from whatever
frames have landed. That is a legitimate average, since record_stats keeps a
running mean from the start of the statistics phase, but it is a mean over a
shorter window and the panel says so rather than presenting the two curves as
if they covered the same interval.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUT = HERE/"sweep_out"
RHO_J, D = 1.6413, 0.0254


def refine(x, y, i):
    if i <= 0 or i >= len(y)-1:
        return x[i]
    d = y[i-1]-2*y[i]+y[i+1]
    return x[i] if d == 0 else x[i]+0.5*(y[i-1]-y[i+1])/d*(x[1]-x[0])


def prof(f):
    d = np.load(OUT/f)
    ks = sorted([k[:-2] for k in d.files if k.endswith("_t")],
                key=lambda q: int("".join(c for c in q if c.isdigit())))
    dx = float(d[ks[0]+"_t"][1])
    r = d[ks[-1]+"_rho"]/RHO_J
    x = (np.arange(r.size)+0.5)*dx/D
    t = float(np.atleast_1d(d[ks[-1]+"_t"])[0])
    return x, r, dx, t, len(ks)


def xfc(x, r, lo=1.70, hi=3.80):
    m = (x >= lo) & (x <= hi) & np.isfinite(r)
    xx, rr = x[m], r[m]
    i0 = int(np.argmin(rr))
    i1 = i0 + int(np.argmax(rr[np.arange(i0, len(rr))]))
    g = np.gradient(rr, xx)
    s = np.arange(i0, i1+1)
    return refine(xx, g, s[int(np.nanargmax(g[s]))]), xx[i0], rr[i0]


L3F = "XFC_N20_tophat_LO__L3.npz"
if not (OUT/L3F).exists():
    L3F = "XFC_N20_L3_partial.npz"
CASES = [("Level 3", L3F, "#0072B2"),
         ("Level 4", "XFC_N20_tophat_LO__L4.npz", "#D55E00")]

print("%-9s %-30s %8s %7s %9s %9s %9s"
      % ("case", "file", "dx[um]", "frames", "t_end[ms]", "x_fc", "rho_min"))
R = []
for lab, f, c in CASES:
    x, r, dx, t, n = prof(f)
    v, xm, rm = xfc(x, r)
    R.append(dict(lab=lab, x=x, r=r, c=c, dx=dx, t=t, n=n, xfc=v, xm=xm, rm=rm))
    print("%-9s %-30s %8.1f %7d %9.3f %9.4f %9.4f"
          % (lab, f, dx*1e6, n, t*1e3, v, rm))
print("\n   x_fc  L3 %.4f   L4 %.4f   difference %+.4f D_e (%+.1f %%)"
      % (R[0]["xfc"], R[1]["xfc"], R[0]["xfc"]-R[1]["xfc"],
         100*(R[0]["xfc"]/R[1]["xfc"]-1)))
PART = "partial" in L3F

plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True,
                     "legend.frameon": False})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.4, 4.2),
                             gridspec_kw=dict(width_ratios=[1.9, 1.0]))
fig.subplots_adjust(left=0.062, right=0.988, top=0.925, bottom=0.135, wspace=0.20)

for q in R:
    lab = "%s,  %.0f $\\mu$m" % (q["lab"], q["dx"]*1e6)
    if PART and q["lab"] == "Level 3":
        lab += ",  mean to %.2f ms" % (q["t"]*1e3)
    a1.plot(q["x"], q["r"], color=q["c"], lw=1.5, zorder=5, label=lab)
    a1.axvline(q["xfc"], color=q["c"], lw=1.0, ls="--", zorder=3)
    a2.plot(q["x"], q["r"], color=q["c"], lw=1.6, zorder=5)
    a2.axvline(q["xfc"], color=q["c"], lw=1.0, ls="--", zorder=3)
    a2.plot([q["xm"]], [q["rm"]], marker="v", color=q["c"], ms=7, mec="w",
            mew=0.7, zorder=7)

a1.set_xlim(0, 7.0)
a1.set_ylim(0, 5.0)
a1.set_xlabel("$x/D_e$")
a1.set_ylabel(r"$\overline{\rho}/\rho_j$   on the axis")
a1.grid(True, ls="--", lw=0.5, color="0.93")
a1.set_axisbelow(True)
a1.set_title(r"(a)  NPR$_0$ = 20 top hat, centreline mean density", fontsize=10, pad=5)
a1.legend(loc="upper right", fontsize=9.0, handlelength=1.9)

lo = min(q["xfc"] for q in R) - 0.7
hi = max(q["xfc"] for q in R) + 0.7
a2.set_xlim(lo, hi)
a2.set_ylim(0.15, 1.20)
a2.set_xlabel("$x/D_e$")
a2.set_ylabel(r"$\overline{\rho}/\rho_j$")
a2.grid(True, ls="--", lw=0.5, color="0.93")
a2.set_axisbelow(True)
a2.set_title("(b)  the Mach-disk region", fontsize=10, pad=5)
a2.text(0.035, 0.93,
        "$x_{fc}$   L3 %.4f\n        L4 %.4f\n$\\Delta$      %+.4f $D_e$"
        % (R[0]["xfc"], R[1]["xfc"], R[0]["xfc"]-R[1]["xfc"]),
        transform=a2.transAxes, va="top", fontsize=9, color="0.25",
        linespacing=1.5)

fig.savefig(HERE/"n20_grid_L3_L4.png", dpi=250)
fig.savefig(HERE/"n20_grid_L3_L4.pdf")
print("\nwrote n20_grid_L3_L4.png/pdf")
