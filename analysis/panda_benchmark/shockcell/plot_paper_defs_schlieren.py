#!/usr/bin/env python3
"""Definition/illustration figure for the shock-cell metrics.
(a) Centreline mean pressure <p>/p_inf of the 3D L5 solution with the
    operational definitions marked: the expansion branch, the steepest
    recompression x_fc = argmax d<p>/dx, and the successive principal
    mean-pressure maxima x_n searched DOWNSTREAM of x_fc
    (<p>/p_inf > 1, prominence > 0.03 in normalised pressure, separation
    >= 0.5 D_e; discrete locations estimated to sub-cell precision by a
    three-point parabolic fit).
(b) Time-mean numerical schlieren exp(-k |grad<rho>|/|grad<rho>|_p99.5)
    of the same solution (z = 0 plane, in-plane gradient of the mean
    density), with x_fc, x_n and the peak-to-peak shock-cell spacings
    S_n = x_{n+1} - x_n annotated."""
import numpy as np
from scipy.signal import find_peaks
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams.update({
    "font.size": 9, "axes.linewidth": 0.6,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
})
HERE = Path(__file__).parent
D = 0.0254; P_AMB = 99780.0
XMAX = 5.2
CFC, CPK = "#1f4e9c", "#b3251e"

def boxfilt(a, x, wD):
    w = max(3, int(round(wD / (x[1] - x[0]))))
    if w % 2 == 0:
        w += 1
    out = np.convolve(a, np.ones(w) / w, mode="same")
    out[:w//2] = a[:w//2]; out[-(w//2):] = a[-(w//2):]
    return out

def refine(x, y, i):
    if i <= 0 or i >= len(y) - 1:
        return x[i]
    d = (y[i-1] - 2*y[i] + y[i+1])
    return x[i] if d == 0 else x[i] + 0.5*(y[i-1] - y[i+1])/d * (x[1]-x[0])

# ---- metrics from the 3D L5 axis mean pressure (full window) ----
d = np.load(HERE / "blockax3d_L5.npz")
prof = d["prof"][-1]
x = (np.arange(prof.shape[1]) + 0.5) * float(d["dx"]) / D
p = prof[0] / P_AMB
ps = boxfilt(p, x, 0.05)
m = (x > 0.4) & (x < 1.6)
dp = np.gradient(ps, x)
i0 = np.where(m)[0]
xfc = refine(x, dp, i0[np.nanargmax(dp[i0])])
dist = max(1, int(round(0.5 / (x[1]-x[0]))))
pk, _ = find_peaks(np.where(x > xfc, ps, -np.inf), height=1.0,
                   prominence=0.03, distance=dist)
xn = [refine(x, ps, i) for i in pk[:4]]
print("xfc=%.3f  xn=%s  Sn=%s" % (xfc, ["%.3f" % v for v in xn],
                                  ["%.3f" % s for s in np.diff(xn)]))

# ---- mean schlieren: per-level gradient composite (remote extraction) ----
sm = np.load(HERE / "schl3d_L5_MEAN.npz")
g = sm["grho"]
dxs = float(sm["dx"]) / D
xs = (np.arange(g.shape[0]) + 0.5) * dxs
ys = (np.arange(g.shape[1]) - g.shape[1]//2 + 0.5) * dxs
# axially adaptive shading reference: sliding-window (1 D_e) p99.5,
# equalises contrast between the sharp near-field and the diffuse
# time-mean downstream cells (stated in the caption)
ref = np.empty(g.shape[0])
for ii in range(g.shape[0]):
    a = max(0, ii - int(0.5/dxs)); b = min(g.shape[0], ii + int(0.5/dxs))
    ref[ii] = np.nanpercentile(g[a:b], 99.5)
ref = np.maximum(ref, 0.12*np.nanmax(ref))
SHL = np.exp(-5.0 * g / ref[:, None])

fig, (ax0, ax1) = plt.subplots(
    2, 1, figsize=(7.0, 4.6), sharex=True,
    gridspec_kw=dict(height_ratios=[1.0, 1.05]))
mm = (x > 0.05) & (x <= XMAX)
ax0.plot(x[mm], ps[mm], color="k", lw=1.1)
ax0.axhline(1.0, color="k", lw=0.5, ls=":")
ax0.axvline(xfc, color=CFC, lw=0.9, ls="--")
for k, xv in enumerate(xn):
    yv = np.interp(xv, x, ps)
    ax0.plot(xv, yv, "o", color=CPK, ms=4.5, mec="k", mew=0.4, zorder=5)
    ax0.annotate("$x_%d$" % (k+1), (xv, yv), xytext=(xv+0.04, yv+0.10),
                 fontsize=8.5, color=CPK)
ax0.annotate("$x_{fc}$", (xfc, 0.62), xytext=(xfc-0.52, 0.52),
             fontsize=8.5, color=CFC,
             arrowprops=dict(arrowstyle="-", color=CFC, lw=0.7))
ax0.annotate("expansion", (0.45, 0.95), xytext=(0.12, 0.62), fontsize=8,
             color="0.3",
             arrowprops=dict(arrowstyle="-", color="0.3", lw=0.6))
ax0.annotate("steepest\nrecompression", (xfc-0.02, 0.86),
             xytext=(1.15, 0.50), fontsize=8, color=CFC,
             arrowprops=dict(arrowstyle="-", color=CFC, lw=0.6))
ax0.set_ylabel("$\\langle p\\rangle/p_\\infty$", labelpad=2)
ax0.set_ylim(np.nanmin(ps[mm]) - 0.05, np.nanmax(ps[mm]) + 0.14)
ax0.tick_params(length=2.5)
ax0.text(0.012, 0.94, "(a)", transform=ax0.transAxes, ha="left", va="top",
         fontsize=9)

ax1.imshow(SHL.T, origin="lower", extent=[xs[0], xs[-1], ys[0], ys[-1]],
           vmin=0, vmax=1, cmap="gray", aspect="equal",
           interpolation="nearest", rasterized=True)
ax1.axvline(xfc, color=CFC, lw=0.9, ls="--")
for xv in xn:
    ax1.axvline(xv, color=CPK, lw=0.8, ls=":")
for k in range(len(xn) - 1):
    ax1.annotate("", (xn[k+1], 0.66), (xn[k], 0.66),
                 arrowprops=dict(arrowstyle="<->", color=CPK, lw=0.9))
    ax1.text(0.5*(xn[k]+xn[k+1]), 0.72, "$S_%d$" % (k+1), ha="center",
             va="bottom", fontsize=8.5, color=CPK)
ax1.set_xlim(0, XMAX); ax1.set_ylim(-0.85, 0.85)
ax1.set_xlabel("$x/D_e$", labelpad=2)
ax1.set_ylabel("$y/D_e$", labelpad=2)
ax1.tick_params(length=2.5)
ax1.text(0.012, 0.94, "(b)", transform=ax1.transAxes, ha="left", va="top",
         fontsize=9, color="k",
         bbox=dict(facecolor="white", edgecolor="none", alpha=0.6, pad=1.0))
fig.align_ylabels([ax0, ax1])
fig.tight_layout(h_pad=0.4)
fig.savefig(HERE / "paperfig_metricdefs_schlieren.png", dpi=300)
fig.savefig(HERE / "paperfig_metricdefs_schlieren.pdf")
print("wrote paperfig_metricdefs_schlieren")
