#!/usr/bin/env python3
"""E_n curvature test + shear-layer growth correlation (zero-cost follow-up).

Left:  E_n vs x_n with a pure-scale (through-origin) reference line and the
       grid-structure landmarks: L5 shear-layer corridor end (0.75D),
       axis L4->L3 (2.0D), L3->L2 (4.5D).
Right: vorticity thickness delta_omega(x) (level-artifact-filtered) and the
       local cell-length ratio L_n^CFD/L_n^exp on the same axis.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/home/qiaoj/testcerisse/cerisse/analysis/panda_benchmark/shockcell"
D = 0.0254
GAM, RGAS = 1.4, 287.0

XE = np.array([1.20, 2.60, 3.80, 4.90, 5.90, 6.60, 7.30])
XC = np.array([1.188, 2.478, 3.587, 4.580, 5.517, 6.390, 7.202])
EN = XC - XE
LE = np.diff(XE)
LC = np.diff(XC)
ratio = LC / LE

d = np.load(f"{BASE}/m142_stats_axis_plane.npz")
P = d["plane"]; CPL = list(d["comps_plane"])
NXP, NZP = P.shape[1], P.shape[2]
xp_ = (np.arange(NXP) + 0.5) * 24.0 / NXP
zp_ = -8.0 + (np.arange(NZP) + 0.5) * 16.0 / NZP
up = P[CPL.index(0)]
kc = NZP // 2
zpos = zp_[kc:]

# delta_omega with gradient computed EXCLUDING the level-boundary contaminated
# columns: recompute per x, then median-filter and mask +-0.15D around known
# plane-level transitions (4.5D, 8.0D at plane resolution)
dom = np.zeros(NXP)
for i in range(NXP):
    prof = up[i, kc:]
    g = np.abs(np.gradient(prof, zpos * D))
    dom[i] = (prof.max()) / max(g.max(), 1e-9) / D
from scipy.signal import medfilt
dom_f = medfilt(dom, 9)
mask = np.ones(NXP, bool)
for xb in (4.5, 8.0):
    mask &= ~(np.abs(xp_ - xb) < 0.2)

fig, axs = plt.subplots(1, 2, figsize=(15, 5.4), constrained_layout=True)
ax = axs[0]
ax.axhline(0, color="0.6")
s0 = np.sum(EN[:4] * XE[:4]) / np.sum(XE[:4] ** 2)   # through-origin LSQ on n<=4
xx = np.linspace(0, 8, 50)
ax.plot(xx, s0 * xx, "b--", lw=1.2,
        label=f"pure-scale reference (through origin, slope {s0*100:.1f}%)")
ax.errorbar(XE, EN, yerr=[0.05, 0.05, 0.05, 0.05, 0.08, 0.10, 0.10],
            fmt="ko-", capsize=3, label="measured E_n")
for xb, lab, c in ((0.75, "L5 corridor end", "green"),
                   (2.0, "axis L4-L3", "orange"),
                   (4.5, "axis L3-L2", "red")):
    ax.axvline(xb, color=c, ls=":", lw=1.5)
    ax.text(xb + 0.05, ax.get_ylim()[0] * 0.15, lab, color=c, rotation=90,
            fontsize=9, va="bottom")
ax.set_xlabel("x_exp / D"); ax.set_ylabel("E_n [D]")
ax.set_title("cumulative phase error: curvature vs pure scaling")
ax.legend(loc="lower left", fontsize=9)

ax = axs[1]
ax.plot(xp_[mask], dom_f[mask], "m-", lw=1.4, label="delta_omega(x)/D (filtered)")
ax.set_xlabel("x/D"); ax.set_ylabel("delta_omega / D"); ax.set_xlim(0, 8)
ax.set_ylim(0, 0.6)
for xb, c in ((0.75, "green"), (2.0, "orange"), (4.5, "red")):
    ax.axvline(xb, color=c, ls=":", lw=1.5)
ax2 = ax.twinx()
xmid = 0.5 * (XE[:-1] + XE[1:])
ax2.plot(xmid, ratio, "ks-", ms=6, label="L_n CFD/exp")
ax2.axhline(1.0, color="0.6", lw=0.8)
ax2.set_ylabel("L_n CFD / exp"); ax2.set_ylim(0.7, 1.4)
ax.set_title("shear-layer thickness vs local cell-length ratio")
ax.legend(loc="upper left", fontsize=9); ax2.legend(loc="upper right", fontsize=9)
fig.savefig(f"{BASE}/fig4_curvature_deltaomega.png", dpi=140)

print("through-origin slope (n<=4):", f"{s0*100:.2f}% per D")
print("local deficit E_n/x_n [%]:", np.round(EN / XE * 100, 2))
print("delta_omega at x=0.75/1.5/2.5/3.5/4.4 D:",
      [round(float(np.interp(x, xp_[mask], dom_f[mask])), 3) for x in (0.75, 1.5, 2.5, 3.5, 4.4)])
grow1 = np.polyfit(xp_[(xp_ > 0.2) & (xp_ < 0.75)], dom_f[(xp_ > 0.2) & (xp_ < 0.75)], 1)[0]
grow2 = np.polyfit(xp_[(xp_ > 1.0) & (xp_ < 2.0)], dom_f[(xp_ > 1.0) & (xp_ < 2.0)], 1)[0]
grow3 = np.polyfit(xp_[(xp_ > 2.2) & (xp_ < 4.2)], dom_f[(xp_ > 2.2) & (xp_ < 4.2)], 1)[0]
print(f"d(delta_omega)/dx: [0.2-0.75D]={grow1:.4f}  [1-2D]={grow2:.4f}  [2.2-4.2D]={grow3:.4f}")
