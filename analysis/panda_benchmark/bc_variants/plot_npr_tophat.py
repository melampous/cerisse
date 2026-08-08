#!/usr/bin/env python3
"""Top-hat inlet at NPR 3.27, 15 and 20: centreline mean density and the
first-cell closure against the Crist scaling.

All three are axisymmetric Level 4 with the same domain, base grid, scheme and
CFL; only p_jet differs, through one literal in prob.h. Statistics are the
running mean over [0.4, 2.4] ms.

The Crist correlation x_s1/D_e = 0.67 sqrt(NPR_0) is an asymptotic
high-pressure-ratio result. Plotting the three points against it shows where
that asymptote is reached, which is the reason for running NPR 15 and 20 at
all: at the Panda condition the measured closure sits well inside the
correlation, and the question is whether that gap is a property of the
correlation or of the simulation.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
OUT = HERE/"sweep_out"
RHO_J, D = 1.6413, 0.0254

CASES = [(3.2734, "XFC_A0_tophat_L4.npz",           0.35, 1.30),
         (15.0,   "XFC_N15_tophat_LO__L4.npz",      1.40, 3.20),
         (20.0,   "XFC_N20_tophat_LO__L4.npz",      1.70, 3.80)]
COL = plt.get_cmap("viridis")(np.linspace(0.10, 0.78, len(CASES)))


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
    return (np.arange(r.size)+0.5)*dx/D, r


def xfc(x, r, lo, hi):
    m = (x >= lo) & (x <= hi) & np.isfinite(r)
    xx, rr = x[m], r[m]
    i0 = int(np.argmin(rr))
    i1 = i0 + int(np.argmax(rr[np.arange(i0, len(rr))]))
    g = np.gradient(rr, xx)
    s = np.arange(i0, i1+1)
    return refine(xx, g, s[int(np.nanargmax(g[s]))]), xx[i0], rr[i0]


R = []
print("%-8s %9s %9s %9s %11s %8s" % ("NPR_0", "x_min", "rho_min", "x_fc",
                                     "0.67 sqrtNPR", "ratio"))
for npr, f, lo, hi in CASES:
    x, r = prof(f)
    v, xm, rm = xfc(x, r, lo, hi)
    cr = 0.67*np.sqrt(npr)
    R.append(dict(npr=npr, x=x, r=r, xfc=v, xm=xm, rm=rm, crist=cr))
    print("%-8.4f %9.4f %9.4f %9.4f %11.4f %8.3f" % (npr, xm, rm, v, cr, v/cr))

plt.rcParams.update({"font.size": 8.5, "axes.linewidth": 0.6,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True,
                     "legend.frameon": False, "axes.titlesize": 9})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.2, 3.3),
                             gridspec_kw=dict(width_ratios=[1.72, 1.0]))
fig.subplots_adjust(left=0.062, right=0.988, top=0.905, bottom=0.145, wspace=0.245)

for q, col in zip(R, COL):
    m = q["x"] <= 8.0
    a1.plot(q["x"][m], q["r"][m], color=col, lw=1.3, zorder=5,
            label=r"NPR$_0$ = %.4g" % q["npr"])
    a1.plot([q["xfc"]], [np.interp(q["xfc"], q["x"], q["r"])], marker="v",
            color=col, ms=7, mec="w", mew=0.7, zorder=7)
a1.set_xlim(0, 8.0)
a1.set_ylim(0, 5.0)
a1.set_xlabel("$x/D_e$")
a1.set_ylabel(r"$\overline{\rho}/\rho_j$   on the axis")
a1.grid(True, ls="--", lw=0.4, color="0.92")
a1.set_axisbelow(True)
a1.set_title("(a)  centreline mean density, top-hat inlet, Level 4", pad=4)
a1.legend(loc="upper right", fontsize=8.2, handlelength=1.8, labelspacing=0.35,
          title=r"$\blacktriangledown$  $x_{fc}$")
a1.get_legend().get_title().set_fontsize(8.2)

nn = np.linspace(2.0, 23.0, 200)
a2.plot(nn, 0.67*np.sqrt(nn), "k--", lw=1.0, zorder=3,
        label=r"Crist, $0.67\sqrt{\mathrm{NPR}_0}$")
for q, col in zip(R, COL):
    a2.plot([q["npr"]], [q["xfc"]], marker="o", ms=7, color=col, mec="k",
            mew=0.6, zorder=6)
a2.set_xlim(0, 23)
a2.set_ylim(0, 3.4)
a2.set_xlabel(r"NPR$_0$")
a2.set_ylabel("$x_{fc}/D_e$")
a2.grid(True, ls="--", lw=0.4, color="0.92")
a2.set_axisbelow(True)
a2.set_title("(b)  first-cell closure against the Crist scaling", pad=4)
for q in R:
    a2.annotate("%+.0f %%" % (100*(q["xfc"]/q["crist"]-1)),
                (q["npr"], q["xfc"]), textcoords="offset points",
                xytext=(7, -9), fontsize=8.0, color="0.3")
a2.legend(loc="upper left", fontsize=8.2, handlelength=2.0)

fig.savefig(HERE/"npr_tophat_xfc.pdf")
fig.savefig(HERE/"npr_tophat_xfc.png", dpi=300)
print("\nwrote npr_tophat_xfc.pdf/png")
