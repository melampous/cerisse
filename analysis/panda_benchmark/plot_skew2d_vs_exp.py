#!/usr/bin/env python3
"""2D-RZ L4 skew-symmetric + JST against Panda & Seasholtz (1999), centreline density.

Simulation: s2d_skewjst_L4 (flanged inlet, sharp top-hat, De/256 near field),
online statistics over [2, 6] ms. The curve is the running mean at t = 6 ms;
the band is the frame-to-frame range of the running mean over t >= 3 ms
(window length >= 1 ms), i.e. the statistical uncertainty of the mean.

Reference curves: the 3D Cartesian L5 skew extract (same scheme and constants,
15.27-17.60 ms window) so the dimensionality effect is visible without mixing
extraction methods (both are axis extracts of DensityMEAN).

Experiment: EPAPS centreline record for the M_j = 1.42 jet, 80 points on
0.60 < x/De < 6.92, normalised by the fully expanded jet density rho_j.
"""
from pathlib import Path
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

HERE = Path(__file__).resolve().parent
D, RHO_J = 0.0254, 1.6413
XLO, XHI = 0.60, 6.92

# ---- 2D L4 skew series --------------------------------------------------
fs = sorted(glob.glob(str(HERE/"skew2d_axis"/"AX2D_plt*.npz")))
A, T, X = [], [], None
for f in fs:
    d = np.load(f)
    if X is None:
        X = d["x"]/D
        A.append(d["ax_DensityMEAN"]/RHO_J)
    else:
        # grids can differ slightly frame to frame (regrid): interp onto first
        A.append(np.interp(X, d["x"]/D, d["ax_DensityMEAN"]/RHO_J))
    T.append(float(d["time"])*1e3)
A = np.array(A); T = np.array(T)
band = A[T >= 3.0]
x2, r2, lo2, hi2 = X, A[-1], band.min(0), band.max(0)

# ---- 3D L5 skew reference ----------------------------------------------
f3 = sorted(glob.glob(str(HERE/"bc_variants"/"sweep_out"/"AX_L5_s3d_skewjst_plt*.npz")))
d3 = np.load(f3[-1])
x3, r3 = d3["x"]/D, d3["ax_DensityMEAN"]/RHO_J
t3 = float(d3["time"])*1e3

# ---- experiment ---------------------------------------------------------
e = np.load(HERE/"panda_m142den_axis_EPAPS.npz")
xe, re = e["x"], e["rho"]

# ---- metrics ------------------------------------------------------------
print("%-30s %8s %9s %9s %9s %8s"
      % ("case", "t_end", "rms dev", "bias", "max|dev|", "corr"))
Q = {}
for name, xs, rs, tend in [("2D-RZ L4 skew+JST", x2, r2, T[-1]),
                           ("3D L5 skew+JST",    x3, r3, t3)]:
    si = np.interp(xe, xs, rs)
    dv = si - re
    Q[name] = dv
    print("%-30s %8.2f %9.4f %+9.4f %9.4f %8.4f"
          % (name, tend, np.sqrt(np.mean(dv*dv)), dv.mean(),
             np.abs(dv).max(), np.corrcoef(si, re)[0, 1]))

print("\nshock-cell density maxima, prominence > 0.10, separation > 0.5 De")
ie, _ = find_peaks(re, prominence=0.10, distance=int(0.5/np.median(np.diff(xe))))
w = (x2 >= XLO) & (x2 <= 7.0)
xs, rs = x2[w], r2[w]
isim, _ = find_peaks(rs, prominence=0.10, distance=int(0.5/np.median(np.diff(xs))))
print("  %-3s %9s %9s %9s %9s" % ("n", "x_exp", "x_2Dskew", "dx", "dx [%]"))
for k in range(min(len(ie), len(isim))):
    print("  %-3d %9.3f %9.3f %+9.3f %+9.2f"
          % (k+1, xe[ie[k]], xs[isim[k]], xs[isim[k]]-xe[ie[k]],
             100*(xs[isim[k]]-xe[ie[k]])/xe[ie[k]]))
print("  spacing exp     %s" % " ".join("%.3f" % v for v in np.diff(xe[ie])))
print("  spacing 2Dskew  %s" % " ".join("%.3f" % v for v in np.diff(xs[isim])))

# ---- figure -------------------------------------------------------------
plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True,
                     "legend.frameon": False})
fig, (a1, a2) = plt.subplots(2, 1, figsize=(9.6, 6.4), sharex=True,
                             gridspec_kw=dict(height_ratios=[2.05, 1.0]))
fig.subplots_adjust(left=0.082, right=0.982, top=0.955, bottom=0.098,
                    hspace=0.09)

w2 = (x2 >= 0.0) & (x2 <= 7.2)
a1.fill_between(x2[w2], lo2[w2], hi2[w2], color="#009E73", alpha=0.20, lw=0,
                zorder=3)
a1.plot(x2[w2], r2[w2], color="#009E73", ls="-", lw=2.0, zorder=5,
        label="2D-RZ $L4$ skew+JST,  window [2, 6] ms")
w3 = (x3 >= 0.0) & (x3 <= 7.2)
a1.plot(x3[w3], r3[w3], color="#D55E00", ls="--", lw=1.3, zorder=4,
        label="3D $L5$ skew+JST,  $t$ = %.2f ms" % t3)
a1.plot(xe, re, "o", ms=4.6, mfc="w", mec="k", mew=0.9, ls="", zorder=7,
        label="Panda & Seasholtz (1999), $M_j$ = 1.42")
a1.fill_between([], [], [], color="0.6", alpha=0.20, lw=0,
                label="running-mean spread, $t \\geq 3$ ms")
a1.set_ylim(0.2, 1.6)
a1.set_ylabel(r"$\overline{\rho}/\rho_j$   on the axis")
a1.legend(loc="lower right", fontsize=8.8, handlelength=2.2)
a1.text(0.008, 0.955, "(a)", transform=a1.transAxes, va="top", fontsize=10.5)

bandw = np.interp(xe, x2, hi2-lo2)
a2.fill_between(xe, -0.5*bandw, 0.5*bandw, color="0.55", alpha=0.28, lw=0,
                zorder=2, label="statistical spread of the mean")
a2.plot(xe, Q["2D-RZ L4 skew+JST"], color="#009E73", lw=2.0, zorder=5,
        label="2D-RZ $L4$")
a2.plot(xe, Q["3D L5 skew+JST"], color="#D55E00", ls="--", lw=1.3, zorder=4,
        label="3D $L5$")
a2.axhline(0.0, color="0.35", lw=0.9, ls=":", zorder=3)
a2.set_ylim(-0.45, 0.55)
a2.set_xlim(0, 7.2)
a2.set_xlabel("$x/D_e$")
a2.set_ylabel(r"simulation $-$ experiment")
a2.legend(loc="upper left", fontsize=8.4, handlelength=2.2, ncol=3)
a2.text(0.008, 0.94, "(b)", transform=a2.transAxes, va="top", fontsize=10.5)
for ax in (a1, a2):
    ax.grid(True, ls="--", lw=0.5, color="0.93")
    ax.set_axisbelow(True)

fig.savefig(HERE/"skew2d_vs_exp_centreline.png", dpi=250)
fig.savefig(HERE/"skew2d_vs_exp_centreline.pdf")
print("\nwrote skew2d_vs_exp_centreline.png/pdf")
