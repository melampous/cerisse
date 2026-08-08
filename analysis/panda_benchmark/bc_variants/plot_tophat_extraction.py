#!/usr/bin/env python3
"""How x_fc is extracted for the top hat, and why the 3D point is the weakest
number in the set.

Top panel: the time-averaged centreline density through the first shock cell,
for the 3D reference case and for the same inlet in 2D on two grids.
Bottom panel: the raw gradient the extraction maximises, over the search
interval defined by the thesis (first expansion minimum to first downstream
principal maximum).

Two properties of the 3D curve are visible directly. Its samples come in
identical pairs, because the axis data is 198.4 um replicated onto the 99.2 um
sampling grid, so the quadratic refinement is fitting a doubled sequence and
its decimals are not meaningful. And the shock is spanned by six cells of
198.4 um, which is what sets the real precision of the peak location.
"""
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUT = HERE/"sweep_out"
SCR = Path("/tmp/claude-1000/-home-qiaoj-testcerisse-cerisse/"
           "efba2e60-081b-46b5-af00-eb8de7ed5170/scratchpad")
D, RHO_J, T0 = 0.0254, 1.6413, 4.0e-4


def refine(x, y, i):
    if i <= 0 or i >= len(y)-1:
        return x[i]
    den = y[i-1]-2*y[i]+y[i+1]
    return x[i] if den == 0 else x[i]+0.5*(y[i-1]-y[i+1])/den*(x[1]-x[0])


def extract(x, r):
    m = np.isfinite(r)
    e = np.where((x >= 0.35) & (x <= 1.0) & m)[0]
    i0 = e[np.nanargmin(r[e])]
    c = np.where((np.arange(x.size) > i0) & (x <= 1.7) & m)[0]
    i1 = c[np.nanargmax(r[c])]
    g = np.gradient(r, x)
    s = np.arange(i0, i1+1)
    ip = s[np.nanargmax(g[s])]
    return dict(i0=i0, i1=i1, ip=ip, g=g, xfc=refine(x, g, ip))


def load2d(fn):
    d = np.load(OUT/fn)
    ks = sorted([k[:-2] for k in d.files if k.endswith("_t")],
                key=lambda q: int("".join(c for c in q if c.isdigit())))
    dx = float(d[ks[0]+"_t"][1])
    M = np.array([d[k+"_rho"] for k in ks])/RHO_J
    return (np.arange(M.shape[1])+0.5)*dx/D, M[-1], dx*1e6


s = np.load(SCR/"STATS3D_tophat.npz")
CASES = [("3D reference,  axis data 198.4 $\\mu$m", s["xx"],
          s["ax_DensityMEAN"][0]/RHO_J, "k", "-", 1.7),
         ("2D L3,  198.4 $\\mu$m everywhere", *load2d("XFC_A0_tophat_LO__L3.npz")[:2],
          "#CC79A7", "-", 1.2),
         ("2D L4,  lip 99.2 $\\mu$m", *load2d("XFC_A0_tophat_L4.npz")[:2],
          "#0072B2", "-", 1.2)]

print("%-34s %8s %8s %8s %8s %10s"
      % ("case", "x_min", "rho_min", "x_max", "x_fc", "cells in shock"))
res = []
for lab, x, r, col, ls, lw in CASES:
    E = extract(x, r)
    dx = x[1]-x[0]
    # distinct samples between the expansion minimum and 90 % of the rise
    seg = r[E["i0"]:E["i1"]]
    tgt = r[E["i0"]] + 0.9*(r[E["i1"]]-r[E["i0"]])
    j = np.argmax(seg >= tgt)
    ncell = len(np.unique(np.round(seg[:j+1], 9)))
    res.append((lab, x, r, col, ls, lw, E, ncell))
    print("%-34s %8.4f %8.4f %8.4f %8.4f %10d"
          % (lab.split(",")[0], x[E["i0"]], r[E["i0"]], x[E["i1"]], E["xfc"], ncell))

plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})
fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.4, 6.6), sharex=True,
                             gridspec_kw={"hspace": 0.08, "height_ratios": [1.15, 1]})

for lab, x, r, col, ls, lw, E, nc in res:
    m = (x >= 0.72) & (x <= 1.55)
    a1.plot(x[m], r[m], color=col, lw=lw, ls=ls, label=lab)
    a1.plot(x[E["i0"]], r[E["i0"]], "v", color=col, ms=6, mec="w", mew=0.6)
    a2.plot(x[m], E["g"][m], color=col, lw=lw, ls=ls)
    a2.axvline(E["xfc"], color=col, lw=0.9, ls=":", alpha=0.8)
    a2.plot(E["xfc"], E["g"][E["ip"]], "o", color=col, ms=6, mec="w", mew=0.7)
    a2.annotate("$x_{fc}$ = %.4f" % E["xfc"], (E["xfc"], E["g"][E["ip"]]),
                textcoords="offset points", xytext=(7, 3), fontsize=7.4,
                color=col)

# the 3D samples, shown individually to make the pairing visible
x3, r3 = res[0][1], res[0][2]
m3 = (x3 >= 0.895) & (x3 <= 0.975)
a1.plot(x3[m3], r3[m3], "o", color="k", ms=3.2, mfc="none", mew=0.7, zorder=6)
a1.annotate("3D samples come in identical pairs:\n"
            "198.4 $\\mu$m data on a 99.2 $\\mu$m grid.\n"
            "The shock spans six 198.4 $\\mu$m cells.",
            (0.935, 0.62), xytext=(1.03, 0.50), textcoords="data", fontsize=7.4,
            color="0.2", arrowprops=dict(arrowstyle="->", lw=0.8, color="0.5"))

a1.set_ylabel(r"$\overline{\rho}/\rho_j$  on the axis")
a1.legend(fontsize=7.6, frameon=False, loc="upper left", handlelength=2.0)
a1.grid(True, ls="--", lw=0.5, color="0.93"); a1.set_axisbelow(True)
a1.text(0.985, 0.05, "$\\blacktriangledown$  first expansion minimum "
        "(start of the search interval)", transform=a1.transAxes, ha="right",
        fontsize=7, color="0.3")

a2.set_xlabel("$x/D_e$")
a2.set_ylabel(r"$\mathrm{d}\overline{\rho}/\mathrm{d}x$   (raw, unsmoothed)")
a2.grid(True, ls="--", lw=0.5, color="0.93"); a2.set_axisbelow(True)
a2.set_xlim(0.72, 1.55)
a2.text(0.985, 0.93, "the extraction takes the maximum of this curve inside "
        "the search interval", transform=a2.transAxes, ha="right", fontsize=7,
        color="0.3")

fig.tight_layout()
fig.savefig(HERE/"tophat_extraction.png", dpi=260)
fig.savefig(HERE/"tophat_extraction.pdf")
print("\nwrote tophat_extraction.png/pdf")
