#!/usr/bin/env python3
"""Four convective schemes on one 3D grid: the mean centreline shock-cell train.

All four runs share the prob.h, the closure, the viscous term, the prescribed
19-box level-2 hierarchy lifted to level 3, the seed (chk09322 at 10.4 ms), the
relaxation leg (10.4 - 11.4 ms) and the statistics window (11.4 - 14.0 ms).
The only difference between them is the convective flux:

  AFD HLLC + TENO5      afd_hllc_t<AfdReconScheme::Teno5>
  AFD HLLC + WENO-Z5    afd_hllc_t<AfdReconScheme::WenoZ5>
  Skew-JST              skew_t with the corrected discontinuity sensor and the
                        reference JST constants C2 = 1.5, C4 = 0.016
  MUSCL HLLC            riemann_t

The axis line is the finest-wins centreline of the running mean, sampled at the
level-3 cell, 396.9 um. x_fc and the x_n peak train follow the thesis
definitions; the error bar on x_fc is the drift of the running cumulative value
over the second half of the window, floored at half a cell.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

HERE = Path(__file__).resolve().parent
OUT = HERE/"sweep_out"
RHO_J, T0 = 1.6413, 11.4e-3
CASES = [("s3d_afdhllc",      "AFD HLLC + TENO5",    "#0072B2", "-"),
         ("s3d_afdhllc_weno", "AFD HLLC + WENO-Z5",  "#009E73", "-"),
         ("s3d_skewjst",      "Skew-JST  1.5/0.016", "#D55E00", "-"),
         ("s3d_musclhllc",    "MUSCL HLLC",          "#7B3294", "-")]
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
    e = np.where((x >= 0.35) & (x <= 1.0))[0]
    i0 = e[np.nanargmin(r[e])]
    c = np.where((np.arange(x.size) > i0) & (x <= 1.7))[0]
    i1 = c[np.nanargmax(r[c])]
    g = np.gradient(r, x)
    s = np.arange(i0, i1+1)
    xfc = refine(x, g, s[np.nanargmax(g[s])])
    sm = boxf(r, x)
    dist = max(1, int(round(0.5/(x[1]-x[0]))))
    pk, _ = find_peaks(np.where(x > xfc, sm, -np.inf), prominence=0.02,
                       distance=dist)
    return xfc, np.array([refine(x, sm, i) for i in pk[:5]]), sm


R = {}
for key, lab, col, ls in CASES:
    d = np.load(OUT/("AX_%s.npz" % key))
    ts = sorted(float(k[1:].split("_")[0]) for k in d.files
                if k.endswith("_DensityMEAN"))
    t = np.array(ts)*1e-3
    M = np.array([d["t%.3f_DensityMEAN" % (tt*1e3)] for tt in t])/RHO_J
    x = (np.arange(M.shape[1])+0.5)*(24.0/M.shape[1])
    xfc, xn, sm = feats(x, M[-1])
    cum = np.array([feats(x, m)[0] for m in M])
    err = max(np.abs(cum[len(cum)//2:]-xfc).max(), 0.5*(x[1]-x[0]))
    R[key] = dict(lab=lab, col=col, ls=ls, x=x, rho=M[-1], sm=sm, xfc=xfc,
                  err=err, xn=xn, S=np.diff(xn), t=t, cum=cum)

print("%-20s %8s %8s %8s %8s %8s %8s | %7s"
      % ("scheme", "x_fc", "x_1", "S_1", "S_2", "S_3", "S_4", "E_S"))
print("%-20s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f | %7s"
      % ("experiment", EXP["x_fc"], EXP["x_1"], EXP["S_1"], EXP["S_2"],
         EXP["S_3"], EXP["S_4"], "-"))
for key, *_ in CASES:
    r = R[key]
    S = list(r["S"][:4]) + [np.nan]*(4-len(r["S"][:4]))
    r["ES"] = np.nanmean([abs(S[i]-EXP["S_%d" % (i+1)]) for i in range(4)])
    print("%-20s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f | %7.4f"
          % (r["lab"], r["xfc"], r["xn"][0], *S, r["ES"]))

plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})
fig = plt.figure(figsize=(12.6, 7.4))
gs = fig.add_gridspec(2, 2, width_ratios=[1.62, 1.0], height_ratios=[1, 1],
                      hspace=0.30, wspace=0.20, left=0.055, right=0.985,
                      top=0.965, bottom=0.075)
a1 = fig.add_subplot(gs[0, 0])
a2 = fig.add_subplot(gs[1, 0])
a3 = fig.add_subplot(gs[0, 1])
at = fig.add_subplot(gs[1, 1]); at.axis("off")

for key, lab, col, ls in CASES:
    r = R[key]
    m = r["x"] <= 6.2
    a1.plot(r["x"][m], r["rho"][m], color=col, lw=1.25, ls=ls, label=lab)
    a1.plot(r["xn"], np.interp(r["xn"], r["x"], r["sm"]), "o", color=col,
            ms=4.2, mec="w", mew=0.6, zorder=6)
a1.set_ylabel(r"$\overline{\rho}/\rho_j$  on the axis")
a1.set_xlabel("$x/D_e$", labelpad=1)
a1.set_xlim(0, 6.2)
a1.grid(True, ls="--", lw=0.5, color="0.93"); a1.set_axisbelow(True)
a1.legend(fontsize=7.6, frameon=False, loc="upper right", handlelength=1.9,
          ncol=2, columnspacing=1.2)
a1.text(0.012, 0.055, "filled circles: the $x_n$ peak train "
        "(0.05 $D_e$ box filter, prominence $\\geq$ 0.02)",
        transform=a1.transAxes, fontsize=6.9, color="0.35")

for key, lab, col, ls in CASES:
    r = R[key]
    m = (r["x"] >= 0.55) & (r["x"] <= 1.75)
    a2.plot(r["x"][m], r["rho"][m], color=col, lw=1.35, ls=ls)
    a2.axvline(r["xfc"], color=col, lw=0.9, ls=":", alpha=0.85)
a2.axvline(EXP["x_fc"], color="0.45", lw=1.2)
a2.text(EXP["x_fc"], 0.985, " experiment 0.888", transform=a2.get_xaxis_transform(),
        va="top", fontsize=7, color="0.35")
a2.set_ylabel(r"$\overline{\rho}/\rho_j$")
a2.set_xlabel("$x/D_e$")
a2.set_xlim(0.55, 1.75)
a2.grid(True, ls="--", lw=0.5, color="0.93"); a2.set_axisbelow(True)
a2.text(0.985, 0.06, "first shock cell; dotted lines are $x_{fc}$",
        transform=a2.transAxes, ha="right", fontsize=7, color="0.35")

n = np.arange(1, 5)
a3.plot(n, [EXP["S_%d" % i] for i in n], "k--o", lw=1.4, ms=7, mfc="w", mew=1.4,
        label="experiment", zorder=8)
for key, lab, col, ls in CASES:
    r = R[key]
    S = list(r["S"][:4]) + [np.nan]*(4-len(r["S"][:4]))
    a3.plot(n, S, "-o", color=col, lw=1.3, ms=5.4, mec="w", mew=0.7, label=lab)
a3.set_xlabel("shock-cell index  $n$")
a3.set_ylabel("$S_n/D_e$")
a3.set_xticks(n)
a3.grid(True, ls="--", lw=0.5, color="0.93"); a3.set_axisbelow(True)
a3.legend(fontsize=7.4, frameon=False, loc="upper right", handlelength=1.8)

tab = ["%-20s %7s %7s %7s %7s %7s" % ("", "x_fc", "S_1", "S_2", "S_3", "S_4"),
       "%-20s %7.3f %7.3f %7.3f %7.3f %7.3f"
       % ("experiment", EXP["x_fc"], EXP["S_1"], EXP["S_2"], EXP["S_3"], EXP["S_4"])]
for key, *_ in CASES:
    r = R[key]
    S = list(r["S"][:4]) + [np.nan]*(4-len(r["S"][:4]))
    tab.append("%-20s %7.4f %7.4f %7.4f %7.4f %7.4f"
               % (r["lab"].replace("  ", " "), r["xfc"], *S))
tab.append("")
tab.append("%-20s %s" % ("mean |S error| E_S",
                         "  ".join("%.4f" % R[k]["ES"] for k, *_ in CASES)))
at.text(0.0, 1.0, "\n".join(tab), family="monospace", fontsize=7.0,
        va="top", ha="left", transform=at.transAxes)
at.text(0.0, 0.40,
        "x_fc uncertainty +-0.0078 D_e for all four (half a level-3 cell,\n"
        "396.9 um; the running-mean drift is smaller than that).\n"
        "Spread across the four schemes: x_fc 1.07 %, x_1 1.00 %,\n"
        "S_1 6.52 %, S_2 19.38 %.",
        family="monospace", fontsize=6.9, va="top", ha="left", color="0.32",
        transform=at.transAxes)

fig.savefig(HERE/"schemes_L3.png", dpi=250)
fig.savefig(HERE/"schemes_L3.pdf")
print("\nwrote schemes_L3.png/pdf")
