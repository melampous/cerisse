#!/usr/bin/env python3
"""Centreline mean density: numerical-scheme sensitivity at M_j = 1.42.
Experiment (m142den phase-mean, primary) vs 3D L5 solutions with three
numerics: baseline WENO-Z5 Navier-Stokes (AX3T_L5, window 6.5-10.16 ms),
TENO5 Navier-Stokes and Euler/WENO-Z5 (phase-B window 6.5-9.3 ms).

TENO5/Euler centrelines are taken from each l5stat final-frame DensityMEAN
z=0 slice (finest-wins composite, coverage to x=7.2) because the innermost-
ring axis in STATS3D_{teno,euler} is NaN for x>2 (finest-owned mask);
baseline uses AX3T_L5. 0.05 D_e display smoothing, rho_j = 1.6413.

Two panels:
 (a) overview of the four curves (distinct line styles, not colour alone);
 (b) deviation from the baseline scheme -- the scheme-to-scheme spread hugs
     zero while the experiment departs several times further, i.e. the
     reconstruction scheme (WENO-Z5 -> TENO5) and the viscous terms
     (NS -> Euler) do not change the primary model-experiment mismatch.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path

plt.rcParams.update({
    "font.size": 9.5, "axes.linewidth": 0.7,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
})
HERE = Path(__file__).parent
D, RHO_J = 0.0254, 1.6413
XMAX = 7.25
CB, CT, CE = "#aab0b7", "#0b6ef5", "#ff6a00"     # baseline (light) / TENO5 / Euler (bright)


def boxfilt(a, x, wD=0.05):
    w = max(3, int(round(wD / (x[1] - x[0]))))
    if w % 2 == 0:
        w += 1
    out = np.convolve(a, np.ones(w) / w, mode="same")
    out[:w // 2] = a[:w // 2]; out[-(w // 2):] = a[-(w // 2):]
    return out


def slice_axis(fn):
    d = np.load(HERE / fn)
    x, y, rho = d["x"], d["y"], d["rho"]
    jc = np.argsort(np.abs(y))[:2]
    r = boxfilt(rho[:, jc].mean(1) / RHO_J, x)
    return x, r


xd = np.load(HERE.parent / "panda_m142den_axis_EPAPS.npz")
a3 = np.load(HERE / "AX3T_L5.npz", allow_pickle=True)
xb = a3["x"] / D
rb = boxfilt(a3["ax_DensityMEAN"] / RHO_J, xb)
xt, rt = slice_axis("SIM_teno_rhoMEAN_z0.npz")
xe, re = slice_axis("SIM_euler_rhoMEAN_z0.npz")

# baseline resampled onto each grid for the deviation panel (rel. baseline)
rb_on_t = np.interp(xt, xb, rb)
rb_on_e = np.interp(xe, xb, rb)
rb_on_x = np.interp(xd["x"], xb, rb)

# third shock convergence = 3rd strong peak of the baseline density gradient
# |d(rho)/dx| (the steep compression rises; weak expansion-side peaks excluded)
from scipy.signal import find_peaks
_grad = boxfilt(np.abs(np.gradient(rb, xb)), xb, 0.08)
_mg = (xb > 0.5) & (xb < 7.5)
_pk, _ = find_peaks(_grad[_mg], height=2.5,
                    distance=int(0.4 / (xb[1] - xb[0])))
_xconv = xb[_mg][_pk]
X3 = float(_xconv[2])          # third convergence (~3.35 De)

fig, (axA, axB) = plt.subplots(
    2, 1, figsize=(7.4, 5.2), sharex=True,
    gridspec_kw={"height_ratios": [2.35, 1.0], "hspace": 0.08})

# ---- (a) overview
for ax in (axA, axB):
    ax.set_axisbelow(True); ax.grid(True, ls="--", lw=0.5, color="0.88")
    ax.tick_params(length=2.5)

mb = (xb > 0.05) & (xb <= XMAX)
mt = (xt > 0.05) & (xt <= XMAX)
me = (xe > 0.05) & (xe <= XMAX)
axA.plot(xb[mb], rb[mb], color=CB, lw=1.9, zorder=4)
axA.plot(xt[mt], rt[mt], color=CT, lw=1.3, ls=(0, (5, 2)), zorder=3)
axA.plot(xe[me], re[me], color=CE, lw=1.3, ls=(0, (1, 1.4)), zorder=3)
axA.plot(xd["x"], xd["rho"], "o", mfc="w", mec="k", mew=0.9, ms=4.2, zorder=5)
axA.set_ylim(0.3, 1.55)
axA.set_ylabel(r"$\langle\rho\rangle/\rho_j$")
axA.axvspan(0, X3, color="0.965", zorder=0)

handles = [
    Line2D([], [], marker="o", ls="", mfc="w", mec="k", ms=4.2,
           label=r"$M_j=1.42$ experiment"),
    Line2D([], [], color=CB, lw=1.9, label="WENO-Z5 NS (baseline)"),
    Line2D([], [], color=CT, lw=1.3, ls=(0, (5, 2)), label="TENO5 NS"),
    Line2D([], [], color=CE, lw=1.3, ls=(0, (1, 1.4)), label="Euler (WENO-Z5)"),
]
axA.legend(handles=handles, fontsize=8, frameon=False, ncol=2,
           loc="upper right", handlelength=2.0, columnspacing=1.4)

# ---- (b) deviation from baseline (baseline = dashed light reference)
axB.axhline(0, color=CB, lw=1.5, ls=(0, (6, 3)), zorder=2)
axB.plot(xt[mt], (rt - rb_on_t)[mt], color=CT, lw=1.5, ls=(0, (5, 2)), zorder=4)
axB.plot(xe[me], (re - rb_on_e)[me], color=CE, lw=1.5, ls=(0, (1, 1.4)), zorder=4)
axB.plot(xd["x"], xd["rho"] - rb_on_x, "o", mfc="w", mec="k", mew=0.9,
         ms=4.2, zorder=5)
axB.set_ylim(-0.34, 0.34)
axB.set_xlim(0, XMAX)
axB.set_xlabel(r"$x/D_e$")
axB.set_ylabel(r"$\Delta\langle\rho\rangle/\rho_j$" "\n(rel. baseline)",
               fontsize=8.5)
axB.text(0.15, 0.28, "(b)", fontsize=9, va="top")
axA.text(0.15, 1.47, "(a)", fontsize=9, va="top")

fig.subplots_adjust(left=0.11, right=0.98, top=0.985, bottom=0.10)
fig.savefig(HERE / "paperfig_cmp_axis_rho_schemes.png", dpi=300)
fig.savefig(HERE / "paperfig_cmp_axis_rho_schemes.pdf")
# quick numeric report
print("max |TENO5-baseline|:", np.nanmax(np.abs((rt-rb_on_t)[mt])).round(3))
print("max |Euler-baseline| :", np.nanmax(np.abs((re-rb_on_e)[me])).round(3))
print("max |exp-baseline|   :", np.nanmax(np.abs(xd["rho"]-rb_on_x)).round(3))
print(f"third shock convergence X3 = {X3:.3f} De")
print("wrote paperfig_cmp_axis_rho_schemes")
