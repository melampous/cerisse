#!/usr/bin/env python3
"""One panel per convective scheme: its mean centreline density against the
existing 3D baseline and the Panda & Seasholtz measurement.

The four scheme runs are identical except for the convective flux. They share
the prob.h (inlet a = 254 um shifted tanh), the closure, the viscous term, the
prescribed level-3 hierarchy, the seed at 10.4 ms, the relaxation leg
(10.4 - 11.4 ms) and the statistics window (11.4 - 14.0 ms). Their axis line is
sampled at the level-3 cell, 396.9 um.

Baseline: the existing 3D production run of the same inlet profile with WENO-Z5
(STATS3D_baseline), whose axis data is 198.4 um - the plateau test shows its
samples come in pairs on the 99.2 um sampling grid, so it is one grid level
finer than the four scheme runs, not a like-for-like scheme comparison.

Experiment: Panda & Seasholtz (1999) EPAPS centreline density for M_j = 1.42,
already normalised by rho_j.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

HERE = Path(__file__).resolve().parent
OUT = HERE/"sweep_out"
PB = HERE.parent                                   # analysis/panda_benchmark
RHO_J = 1.6413
CASES = [("s3d_afdhllc",      "AFD HLLC + TENO5",    "#0072B2"),
         ("s3d_afdhllc_weno", "AFD HLLC + WENO-Z5",  "#009E73"),
         ("s3d_skewjst",      "Skew-JST  $C_2$=1.5, $C_4$=0.016", "#D55E00"),
         ("s3d_musclhllc",    "MUSCL HLLC",          "#7B3294")]
EXP = dict(x_fc=0.888, x_1=1.199, S_1=1.359, S_2=1.186, S_3=1.130, S_4=0.956)


def refine(x, y, i):
    if i <= 0 or i >= len(y)-1:
        return x[i]
    d = y[i-1]-2*y[i]+y[i+1]
    return x[i] if d == 0 else x[i]+0.5*(y[i-1]-y[i+1])/d*(x[1]-x[0])


def boxf(y, x):
    w = max(1, int(round(0.05/(x[1]-x[0]))))
    return np.convolve(y, np.ones(w)/w, mode="same")


def feats(x, r):
    m = np.isfinite(r)
    e = np.where((x >= 0.35) & (x <= 1.0) & m)[0]
    i0 = e[np.nanargmin(r[e])]
    c = np.where((np.arange(x.size) > i0) & (x <= 1.7) & m)[0]
    i1 = c[np.nanargmax(r[c])]
    g = np.gradient(r, x)
    s = np.arange(i0, i1+1)
    xfc = refine(x, g, s[np.nanargmax(g[s])])
    sm = boxf(r, x)
    dist = max(1, int(round(0.5/(x[1]-x[0]))))
    pk, _ = find_peaks(np.where(x > xfc, sm, -np.inf), prominence=0.02,
                       distance=dist)
    return xfc, np.array([refine(x, sm, i) for i in pk[:5]])


R = {}
for key, lab, col in CASES:
    d = np.load(OUT/("AX_%s.npz" % key))
    ts = sorted(float(k[1:].split("_")[0]) for k in d.files
                if k.endswith("_DensityMEAN"))
    M = np.array([d["t%.3f_DensityMEAN" % t] for t in ts])/RHO_J
    x = (np.arange(M.shape[1])+0.5)*(24.0/M.shape[1])
    xfc, xn = feats(x, M[-1])
    R[key] = dict(lab=lab, col=col, x=x, rho=M[-1], xfc=xfc, xn=xn,
                  S=np.diff(xn))

b = np.load(PB/"STATS3D_baseline.npz")
bx = b["xx"]
br = b["ax_DensityMEAN"][0]/RHO_J
bfin = np.isfinite(br)
bxfc, bxn = feats(bx, br)
# the baseline axis array is only valid to x = 2 D_e (512 finite samples), so
# its shock-cell spacings are taken from the thesis table rather than re-derived
bS = np.array([1.300, 1.105, 1.004, 0.968])
print("baseline axis data valid to x = %.3f D_e (%d finite samples)"
      % (bx[bfin][-1], bfin.sum()))

e = np.load(PB/"panda_m142den_axis_EPAPS.npz")
ex, er = e["x"], e["rho"]

print("%-22s %8s %8s %8s %8s %8s | %7s"
      % ("", "x_fc", "S_1", "S_2", "S_3", "S_4", "E_S"))
print("%-22s %8.3f %8.3f %8.3f %8.3f %8.3f | %7s"
      % ("experiment", EXP["x_fc"], EXP["S_1"], EXP["S_2"], EXP["S_3"],
         EXP["S_4"], "-"))
S = list(bS[:4])+[np.nan]*(4-len(bS[:4]))
bES = np.nanmean([abs(S[i]-EXP["S_%d" % (i+1)]) for i in range(4)])
print("%-22s %8.4f %8.4f %8.4f %8.4f %8.4f | %7.4f"
      % ("baseline WENO-Z5 3D", bxfc, *S, bES))
for key, lab, col in CASES:
    r = R[key]
    S = list(r["S"][:4])+[np.nan]*(4-len(r["S"][:4]))
    r["ES"] = np.nanmean([abs(S[i]-EXP["S_%d" % (i+1)]) for i in range(4)])
    print("%-22s %8.4f %8.4f %8.4f %8.4f %8.4f | %7.4f"
          % (lab.split("  ")[0], r["xfc"], *S, r["ES"]))

plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})
fig, axs = plt.subplots(2, 2, figsize=(12.2, 7.4), sharex=True, sharey=True,
                        gridspec_kw={"hspace": 0.13, "wspace": 0.06})

for ax, (key, lab, col) in zip(axs.ravel(), CASES):
    r = R[key]
    mb = bfin & (bx <= 7.0)
    ax.plot(bx[mb], br[mb], color="0.55", lw=1.4, ls="--", zorder=3,
            label="3D baseline, WENO-Z5 (198.4 $\\mu$m, axis data to 2 $D_e$)")
    ax.plot(ex, er, "o", color="k", ms=4.0, mfc="none", mew=0.9, zorder=6,
            label="Panda & Seasholtz 1999")
    m = r["x"] <= 7.0
    ax.plot(r["x"][m], r["rho"][m], color=col, lw=1.6, zorder=5, label=lab)
    ax.axvline(r["xfc"], color=col, lw=1.0, ls=":", zorder=2)
    ax.axvline(EXP["x_fc"], color="k", lw=0.9, alpha=0.45, zorder=2)
    ax.set_xlim(0, 7.0)
    ax.set_ylim(0.30, 1.55)
    ax.grid(True, ls="--", lw=0.5, color="0.93")
    ax.set_axisbelow(True)
    ax.text(0.025, 0.955, lab, transform=ax.transAxes, va="top", fontsize=9.4,
            color=col, fontweight="bold")
    ax.text(0.025, 0.885,
            "$x_{fc}$ = %.4f $\\pm$ 0.0078      baseline %.4f   exp %.3f\n"
            "$S_1$ = %.3f   $S_2$ = %.3f   $S_3$ = %.3f   $S_4$ = %.3f\n"
            "$E_S$ = %.4f      (baseline %.4f, thesis table)"
            % (r["xfc"], bxfc, EXP["x_fc"], *(list(r["S"][:4])+[np.nan]*(4-len(r["S"][:4]))),
               r["ES"], bES),
            transform=ax.transAxes, va="top", fontsize=7.0, color="0.25",
            family="monospace")

for ax in axs[1, :]:
    ax.set_xlabel("$x/D_e$")
for ax in axs[:, 0]:
    ax.set_ylabel(r"$\overline{\rho}/\rho_j$  on the axis")
axs[0, 0].legend(fontsize=7.2, frameon=False, loc="lower right",
                 handlelength=1.9, labelspacing=0.35)
axs[0, 1].text(0.985, 0.05,
               "vertical lines: $x_{fc}$ of the scheme (dotted)\n"
               "and of the experiment (solid grey)",
               transform=axs[0, 1].transAxes, ha="right", va="bottom",
               fontsize=6.8, color="0.4")

fig.savefig(HERE/"schemes_L3_2x2.png", dpi=250, bbox_inches="tight")
fig.savefig(HERE/"schemes_L3_2x2.pdf", bbox_inches="tight")
print("\nwrote schemes_L3_2x2.png/pdf")
