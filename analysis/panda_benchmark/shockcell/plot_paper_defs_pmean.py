#!/usr/bin/env python3
"""Definition/illustration figure for the shock-cell metrics.
(a) Unsmoothed centreline mean density of the 3D L5 solution, with the
    density-based operational closure x_fc marked. The pressure-based
    principal compression locations x_n are overlaid as positional markers.
(b) Time-mean pressure field of the same solution, with x_fc, x_n and the
    peak-to-peak shock-cell spacings S_n = x_{n+1} - x_n annotated."""
import numpy as np
from scipy.signal import find_peaks
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams.update({
    "font.size": 10.5, "axes.linewidth": 0.6,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
})
HERE = Path(__file__).parent
D = 0.0254; P_AMB = 99780.0; RHO_J = 1.6413
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

# ---- density-based x_fc and pressure-based x_n from the 3D L5 mean ----
d = np.load(HERE / "blockax3d_L5.npz")
prof = d["prof"][-1]
x = (np.arange(prof.shape[1]) + 0.5) * float(d["dx"]) / D
p = prof[0] / P_AMB
ps = boxfilt(p, x, 0.05)

dr = np.load(HERE / "AX3T_L5.npz")
xrho = dr["x"] / D
rho = dr["ax_DensityMEAN"] / RHO_J
drhodx = np.gradient(rho, xrho)
expansion = np.where((xrho >= 0.35) & (xrho <= 1.0))[0]
imin = expansion[np.nanargmin(rho[expansion])]
compression = np.where((np.arange(xrho.size) > imin) & (xrho <= 1.65))[0]
imax = compression[np.nanargmax(rho[compression])]
i0 = np.arange(imin, imax + 1)
xfc = refine(xrho, drhodx, i0[np.nanargmax(drhodx[i0])])

dist = max(1, int(round(0.5 / (x[1]-x[0]))))
pk, _ = find_peaks(np.where(x > xfc, ps, -np.inf), height=1.0,
                   prominence=0.03, distance=dist)
xn = [refine(x, ps, i) for i in pk[:4]]
print("xfc=%.3f  xn=%s  Sn=%s" % (xfc, ["%.3f" % v for v in xn],
                                  ["%.3f" % s for s in np.diff(xn)]))

# ---- mean pressure field used to identify x_n ----
sm = np.load(HERE / "pmean3d_L5_z0.npz")
SHL = sm["p"] / P_AMB
dxs = float(sm["dx"]) / D
xs = (np.arange(SHL.shape[0]) + 0.5) * dxs
ys = (np.arange(SHL.shape[1]) - SHL.shape[1]//2 + 0.5) * dxs

fig, (ax0, ax1) = plt.subplots(
    2, 1, figsize=(7.0, 3.9), sharex=True,
    gridspec_kw=dict(height_ratios=[0.52, 1.0]))
mm = (xrho > 0.05) & (xrho <= XMAX)
ax0.plot(xrho[mm], rho[mm], color="k", lw=1.1)
ax0.axvline(xfc, color=CFC, lw=0.9, ls="--")
for k, xv in enumerate(xn):
    yv = np.interp(xv, xrho, rho)
    ax0.plot(xv, yv, "o", color=CPK, ms=4.5, mec="k", mew=0.4, zorder=5)
    ax0.annotate("$x_%d$" % (k+1), (xv, yv), xytext=(xv, yv-0.16),
                 fontsize=10, color=CPK, ha="center", va="top")
ax0.text(xfc + 0.06, np.interp(xfc, xrho, rho) - 0.20, "$x_{fc}$",
         fontsize=10, color=CFC)
ax0.grid(True, which="major", color="0.88", lw=0.5)
ax0.set_axisbelow(True)
ax0.set_ylabel("$\\langle\\rho\\rangle/\\rho_j$", labelpad=2)
ax0.set_ylim(np.nanmin(rho[mm]) - 0.05, np.nanmax(rho[mm]) + 0.14)
ax0.tick_params(length=2.5)
ax0.text(0.012, 0.88, "(a)", transform=ax0.transAxes, ha="left", va="top",
         fontsize=10.5)

im = ax1.imshow(SHL.T, origin="lower", extent=[xs[0], xs[-1], ys[0], ys[-1]],
           vmin=0.4, vmax=1.6, cmap="RdBu_r", aspect="equal",
           interpolation="nearest", rasterized=True)
cb = fig.colorbar(im, ax=ax1, pad=0.012, fraction=0.035)
cb.set_label("$\\langle p\\rangle/p_\\infty$", labelpad=3, fontsize=10)
cb.ax.tick_params(length=2.5, labelsize=9)
ax1.axvline(xfc, color=CFC, lw=0.9, ls="--")
for xv in xn:
    ax1.axvline(xv, color="k", lw=0.8, ls=":")
for k in range(len(xn) - 1):
    ax1.annotate("", (xn[k+1], 0.66), (xn[k], 0.66),
                 arrowprops=dict(arrowstyle="<->", color="k", lw=0.9))
    ax1.text(0.5*(xn[k]+xn[k+1]), 0.60, "$S_%d$" % (k+1), ha="center",
             va="top", fontsize=10, color="k")
for gx in range(1, 6):
    ax1.axvline(gx, color="0.5", lw=0.4, alpha=0.30, zorder=1)
for gy in (-0.5, 0.0, 0.5):
    ax1.axhline(gy, color="0.5", lw=0.4, alpha=0.30, zorder=1)
ax1.set_xlim(0, XMAX); ax1.set_ylim(-0.85, 0.85)
ax1.set_xlabel("$x/D_e$", labelpad=2)
ax1.set_ylabel("$y/D_e$", labelpad=2)
ax1.tick_params(length=2.5)
ax1.text(0.012, 0.94, "(b)", transform=ax1.transAxes, ha="left", va="top",
         fontsize=10.5, color="k",
         bbox=dict(facecolor="white", edgecolor="none", alpha=0.6, pad=1.0))
fig.align_ylabels([ax0, ax1])
fig.tight_layout(h_pad=0.4)
# align panel widths: shrink (a) to the map area of (b) (colorbar excluded)
p1 = ax1.get_position(); p0 = ax0.get_position()
ax0.set_position([p1.x0, p0.y0, p1.width, p0.height])
fig.savefig(HERE / "paperfig_metricdefs_pmean.png", dpi=300)
fig.savefig(HERE / "paperfig_metricdefs_pmean.pdf")
print("wrote paperfig_metricdefs_pmean")
