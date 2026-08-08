#!/usr/bin/env python3
"""Skew-symmetric + JST against Panda & Seasholtz (1999), centreline density.

The measurement is the EPAPS centreline record for the M_j = 1.42 jet,
0.60 < x/D_e < 6.92, 80 points, already normalised by the fully expanded jet
density. The simulation is the level-5 axis extract of the running mean, taken
from the last of nine frames spanning 15.27 to 17.60 ms.

Statistical convergence is shown, not assumed
---------------------------------------------
record_stats keeps a running mean, so the nine frames are nine estimates of the
same quantity with progressively longer windows. The spread between them is the
statistical uncertainty of the curve, and it is plotted as a band. It is small
near the nozzle and large downstream, which is where the flow is unsteady and
the mean converges slowly:

    window          spread (rms)    deviation from experiment (rms)
    0.6 - 3.0            0.017                    0.095
    3.0 - 5.0            0.083                    0.113
    5.0 - 6.92           0.090                    0.063

So the near-field disagreement is five times the statistical scatter and is a
property of the solution; beyond x = 5 the disagreement is smaller than the
scatter and nothing can be concluded there from this averaging window.

AFD-HLLC + TENO5 is drawn alongside because it is the only other scheme with an
axis extract from the same path, at the same time and on the same grid, so the
two can be compared without mixing extraction methods. The LLF pair has no such
extract; adding it would mean comparing an axis extract against a slice.
"""
from pathlib import Path
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

HERE = Path(__file__).resolve().parent
SW = HERE/"bc_variants"/"sweep_out"
D, RHO_J = 0.0254, 1.6413
XLO, XHI = 0.60, 6.92

CASES = [("s3d_skewjst", "skew-symmetric + JST", "#009E73", "-", 2.0),
         ("s3d_afdhllc", "AFD-HLLC + TENO5", "#D55E00", "--", 1.3)]


def series(tag):
    fs = sorted(glob.glob(str(SW/("AX_L5_%s_plt*.npz" % tag))))
    A, T = [], []
    for f in fs:
        d = np.load(f)
        A.append(d["ax_DensityMEAN"]/RHO_J)
        T.append(float(d["time"])*1e3)
    x = np.load(fs[0])["x"]/D
    A = np.array(A)
    return x, A[-1], A.min(0), A.max(0), T


e = np.load(HERE/"panda_m142den_axis_EPAPS.npz")
xe, re = e["x"], e["rho"]

Q = {}
print("%-22s %7s %8s %9s %9s %9s %8s"
      % ("scheme", "frames", "t_end", "rms dev", "bias", "max|dev|", "corr"))
for tag, name, col, ls, lw in CASES:
    x, r, lo, hi, T = series(tag)
    si = np.interp(xe, x, r)
    d = si - re
    Q[tag] = dict(x=x, r=r, lo=lo, hi=hi, T=T, si=si, d=d, name=name,
                  col=col, ls=ls, lw=lw)
    print("%-22s %7d %8.2f %9.4f %+9.4f %9.4f %8.4f"
          % (name, len(T), T[-1], np.sqrt(np.mean(d*d)), d.mean(),
             np.abs(d).max(), np.corrcoef(si, re)[0, 1]))

print("\nby window   (spread = frame-to-frame range of the running mean)")
print("%-22s %-14s %9s %9s %8s" % ("scheme", "window", "spread", "dev rms", "ratio"))
for tag, name, *_ in CASES:
    q = Q[tag]
    for a, b, lab in [(0.6, 3.0, "0.6-3.0"), (3.0, 5.0, "3.0-5.0"),
                      (5.0, 6.92, "5.0-6.92"), (0.6, 6.92, "whole")]:
        w = (q["x"] >= a) & (q["x"] <= b)
        s = np.sqrt(np.mean((q["hi"][w]-q["lo"][w])**2))
        m = (xe >= a) & (xe <= b)
        dr = np.sqrt(np.mean(q["d"][m]**2))
        print("%-22s %-14s %9.4f %9.4f %8.2f" % (name, lab, s, dr, dr/s))

print("\nshock-cell density maxima, prominence > 0.10, separation > 0.5 D_e")
ie, _ = find_peaks(re, prominence=0.10,
                   distance=int(0.5/np.median(np.diff(xe))))
print("  %-3s %9s %9s %9s %9s" % ("n", "x_exp", "x_skew", "dx", "dx [%]"))
for tag, name, *_ in CASES[:1]:
    q = Q[tag]
    w = (q["x"] >= XLO) & (q["x"] <= 7.0)
    xs, rs = q["x"][w], q["r"][w]
    isim, _ = find_peaks(rs, prominence=0.10,
                         distance=int(0.5/(xs[1]-xs[0])))
    for k in range(min(len(ie), len(isim))):
        print("  %-3d %9.3f %9.3f %+9.3f %+9.2f"
              % (k+1, xe[ie[k]], xs[isim[k]], xs[isim[k]]-xe[ie[k]],
                 100*(xs[isim[k]]-xe[ie[k]])/xe[ie[k]]))
    print("  spacing exp   %s"
          % " ".join("%.3f" % v for v in np.diff(xe[ie])))
    print("  spacing skew  %s"
          % " ".join("%.3f" % v for v in np.diff(xs[isim])))

plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True,
                     "legend.frameon": False})
fig, (a1, a2) = plt.subplots(2, 1, figsize=(9.6, 6.4), sharex=True,
                             gridspec_kw=dict(height_ratios=[2.05, 1.0]))
fig.subplots_adjust(left=0.082, right=0.982, top=0.955, bottom=0.098,
                    hspace=0.09)

for tag, name, col, ls, lw in CASES:
    q = Q[tag]
    w = (q["x"] >= 0.0) & (q["x"] <= 7.2)
    a1.fill_between(q["x"][w], q["lo"][w], q["hi"][w], color=col, alpha=0.20,
                    lw=0, zorder=3)
    a1.plot(q["x"][w], q["r"][w], color=col, ls=ls, lw=lw, zorder=5,
            label="%s,  $t$ = %.2f ms" % (name, q["T"][-1]))
    a2.plot(xe, q["d"], color=col, ls=ls, lw=lw, zorder=5, label=name)

a1.plot(xe, re, "o", ms=4.6, mfc="w", mec="k", mew=0.9, ls="", zorder=7,
        label="Panda & Seasholtz (1999), $M_j$ = 1.42")
a1.fill_between([], [], [], color="0.6", alpha=0.20, lw=0,
                label="running-mean spread over the frames")
a1.set_ylim(0.2, 1.6)
a1.set_ylabel(r"$\overline{\rho}/\rho_j$   on the axis")
a1.legend(loc="lower right", fontsize=8.8, handlelength=2.2, ncol=1)
a1.text(0.008, 0.955, "(a)", transform=a1.transAxes, va="top", fontsize=10.5)

# the statistical band of the skew case, mapped onto the deviation panel
q = Q["s3d_skewjst"]
band = np.interp(xe, q["x"], q["hi"]-q["lo"])
a2.fill_between(xe, -0.5*band, 0.5*band, color="0.55", alpha=0.28, lw=0,
                zorder=2, label="statistical spread of the mean")
a2.axhline(0.0, color="0.35", lw=0.9, ls=":", zorder=3)
a2.set_ylim(-0.30, 0.46)
a2.set_xlim(0, 7.2)
a2.set_xlabel("$x/D_e$")
a2.set_ylabel(r"simulation $-$ experiment")
a2.legend(loc="upper left", bbox_to_anchor=(0.048, 1.0), fontsize=8.4,
          handlelength=2.2, ncol=3)
a2.text(0.008, 0.94, "(b)", transform=a2.transAxes, va="top", fontsize=10.5)
for ax in (a1, a2):
    ax.grid(True, ls="--", lw=0.5, color="0.93")
    ax.set_axisbelow(True)

fig.savefig(HERE/"skew_vs_exp_centreline.png", dpi=250)
fig.savefig(HERE/"skew_vs_exp_centreline.pdf")
print("\nwrote skew_vs_exp_centreline.png/pdf")
