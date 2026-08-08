#!/usr/bin/env python3
"""New-tree 4th-order skew-symmetric + JST leg (b2d_skew4_newsrc_L4, flanged
tanh-254um baseline inlet, radial metric-weighted flux split G = rF) against
the old-tree 2D-RZ L4 pair (skew+JST and LLF WENO-Z5, both flanged sharp
top-hat inlet) and the Panda & Seasholtz (1999) centreline density record.

All three simulation legs share the identical grid (De/256 near field, 5-level
AMR), domain and statistics protocol (online statistics over [2, 6] ms; curve
= running mean at t = 6 ms; band = frame-to-frame range of the running mean
over t >= 3 ms).  The new leg differs from the old skew leg in BOTH the code
tree (G = rF radial metric split ported to Skew.h) and the inlet velocity
profile (tanh 254 um vs sharp top-hat), so differences are not attributable
to the scheme port alone.

Experiment: EPAPS centreline record, M_j = 1.42, normalised by rho_j.
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

# ---- old-tree per-frame series (AX2D format) ----------------------------
def series_ax2d(subdir):
    fs = sorted(glob.glob(str(HERE/subdir/"AX2D_plt*.npz")))
    A, T, X = [], [], None
    for f in fs:
        d = np.load(f)
        if X is None:
            X = d["x"]/D
            A.append(d["ax_DensityMEAN"]/RHO_J)
        else:
            A.append(np.interp(X, d["x"]/D, d["ax_DensityMEAN"]/RHO_J))
        T.append(float(d["time"])*1e3)
    A = np.array(A); T = np.array(T)
    band = A[T >= 3.0]
    return X, A[-1], band.min(0), band.max(0), T

# ---- new-tree per-frame series (axis2d_lite AXL format) -----------------
# The phase-2 restart did not set cns.stats_start_time, so time_stat_level
# restarts at cur_time (2 ms) and the on-line mean is normalised by t instead
# of the accumulation window (t - 2 ms).  The exact correction is a per-frame
# factor t/(t - 2 ms), verified against frame-averaged instantaneous fields
# (ratio 4/6 = 0.667 at t = 6 ms, law holds within sampling noise at 3-5 ms).
T_STAT0 = 2.0e-3

def series_axl(subdir):
    fs = sorted(glob.glob(str(HERE/subdir/"AXL_plt*.npz")))
    A, T, X = [], [], None
    for f in fs:
        d = np.load(f)
        if "stat_DensityMEAN" not in d.files:
            continue
        t_s = float(d["time_last"])
        t_ms = t_s*1e3
        if t_ms <= 2.05:               # phase-1 frames carry no statistics
            continue
        m = d["stat_DensityMEAN"] * (t_s/(t_s - T_STAT0))
        if not np.isfinite(m).all() or m.max() <= 0.0:
            continue
        if X is None:
            X = d["z"]/D
            A.append(m/RHO_J)
        else:
            A.append(np.interp(X, d["z"]/D, m/RHO_J))
        T.append(t_ms)
    A = np.array(A); T = np.array(T)
    band = A[T >= 3.0]
    return X, A[-1], band.min(0), band.max(0), T

xs4, rs4, los4, his4, Ts4 = series_axl("skew4newsrc2d_axis")
x2, r2, lo2, hi2, T2 = series_ax2d("skew2d_axis")
xw, rw, low, hiw, Tw = series_ax2d("llfweno2d_axis")

# ---- experiment ---------------------------------------------------------
e = np.load(HERE/"panda_m142den_axis_EPAPS.npz")
xe, re = e["x"], e["rho"]

# ---- metrics ------------------------------------------------------------
print("%-38s %8s %9s %9s %9s %8s"
      % ("case", "t_end", "rms dev", "bias", "max|dev|", "corr"))
Q = {}
for name, xx, rr, tend in [
        ("new-tree skew4+JST (tanh inlet)", xs4, rs4, Ts4[-1]),
        ("old-tree skew+JST (top-hat)",     x2,  r2,  T2[-1]),
        ("old-tree LLF WENO-Z5 (top-hat)",  xw,  rw,  Tw[-1])]:
    si = np.interp(xe, xx, rr)
    dv = si - re
    Q[name] = dv
    print("%-38s %8.2f %9.4f %+9.4f %9.4f %8.4f"
          % (name, tend, np.sqrt(np.mean(dv*dv)), dv.mean(),
             np.abs(dv).max(), np.corrcoef(si, re)[0, 1]))

print("\nshock-cell density maxima, prominence > 0.10, separation > 0.5 De")
ie, _ = find_peaks(re, prominence=0.10, distance=int(0.5/np.median(np.diff(xe))))
for lab, xx, rr in [("new4", xs4, rs4), ("skew", x2, r2), ("weno", xw, rw)]:
    w = (xx >= XLO) & (xx <= 7.0)
    xsw, rsw = xx[w], rr[w]
    isim, _ = find_peaks(rsw, prominence=0.10,
                         distance=int(0.5/np.median(np.diff(xsw))))
    print("  %-4s %-3s %9s %9s %9s %9s" % (lab, "n", "x_exp", "x_sim", "dx", "dx [%]"))
    for k in range(min(len(ie), len(isim))):
        print("  %-4s %-3d %9.3f %9.3f %+9.3f %+9.2f"
              % (lab, k+1, xe[ie[k]], xsw[isim[k]], xsw[isim[k]]-xe[ie[k]],
                 100*(xsw[isim[k]]-xe[ie[k]])/xe[ie[k]]))
    print("  %s spacing sim  %s" % (lab, " ".join("%.3f" % v for v in np.diff(xsw[isim]))))
print("  spacing exp     %s" % " ".join("%.3f" % v for v in np.diff(xe[ie])))

# ---- figure -------------------------------------------------------------
plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True,
                     "legend.frameon": False})
fig, (a1, a2) = plt.subplots(2, 1, figsize=(9.6, 6.4), sharex=True,
                             gridspec_kw=dict(height_ratios=[2.05, 1.0]))
fig.subplots_adjust(left=0.082, right=0.982, top=0.955, bottom=0.098,
                    hspace=0.09)

w4 = (xs4 >= 0.0) & (xs4 <= 7.2)
a1.fill_between(xs4[w4], los4[w4], his4[w4], color="#d1352b", alpha=0.18,
                lw=0, zorder=3)
a1.plot(xs4[w4], rs4[w4], color="#d1352b", ls="-", lw=2.0, zorder=6,
        label="2D-RZ $L4$ skew-4 + JST, new tree ($G=rF$), tanh inlet, [2, 6] ms")
w2 = (x2 >= 0.0) & (x2 <= 7.2)
a1.plot(x2[w2], r2[w2], color="#009E73", ls="--", lw=1.5, zorder=5,
        label="2D-RZ $L4$ skew + JST, old tree, top-hat inlet, [2, 6] ms")
ww = (xw >= 0.0) & (xw <= 7.2)
a1.plot(xw[ww], rw[ww], color="#0b6ef5", ls="-.", lw=1.3, zorder=4,
        label="2D-RZ $L4$ LLF WENO-Z5, old tree, top-hat inlet, [2, 6] ms")
a1.plot(xe, re, "o", ms=4.6, mfc="w", mec="k", mew=0.9, ls="", zorder=7,
        label="Panda & Seasholtz (1999), $M_j$ = 1.42")
a1.fill_between([], [], [], color="0.6", alpha=0.20, lw=0,
                label="running-mean spread, $t \\geq 3$ ms (new leg)")
a1.set_ylim(0.2, 2.95)
a1.set_ylabel(r"$\overline{\rho}/\rho_j$   on the axis")
a1.legend(loc="lower right", fontsize=8.2, handlelength=2.2)
a1.text(0.008, 0.955, "(a)", transform=a1.transAxes, va="top", fontsize=10.5)

band4 = np.interp(xe, xs4, his4-los4)
a2.fill_between(xe, -0.5*band4, 0.5*band4, color="0.55", alpha=0.28, lw=0,
                zorder=2, label="statistical spread of the mean (new leg)")
a2.plot(xe, Q["new-tree skew4+JST (tanh inlet)"], color="#d1352b", lw=2.0,
        zorder=6, label="skew-4 new tree")
a2.plot(xe, Q["old-tree skew+JST (top-hat)"], color="#009E73", ls="--",
        lw=1.4, zorder=5, label="skew old tree")
a2.plot(xe, Q["old-tree LLF WENO-Z5 (top-hat)"], color="#0b6ef5", ls="-.",
        lw=1.2, zorder=4, label="LLF WENO-Z5")
a2.axhline(0.0, color="0.35", lw=0.9, ls=":", zorder=3)
a2.set_ylim(-0.5, 1.95)
a2.set_xlim(0, 7.2)
a2.set_xlabel("$x/D_e$")
a2.set_ylabel(r"simulation $-$ experiment")
a2.legend(loc="upper left", fontsize=8.4, handlelength=2.2, ncol=2)
a2.text(0.008, 0.94, "(b)", transform=a2.transAxes, va="top", fontsize=10.5)
for ax in (a1, a2):
    ax.grid(True, ls="--", lw=0.5, color="0.93")
    ax.set_axisbelow(True)

fig.savefig(HERE/"skew4newsrc_vs_exp_centreline.png", dpi=250)
fig.savefig(HERE/"skew4newsrc_vs_exp_centreline.pdf")
print("\nwrote skew4newsrc_vs_exp_centreline.png/pdf")
