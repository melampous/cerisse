#!/usr/bin/env python3
"""Inflow-profile effect (delta=0 vs delta=254 um) at two resolutions:
centerline time-mean density, L1 (800 um) and L5 (50 um), vs experiment.
Statistics windows: L1 legs [8.0,10.4] ms; L5 baseline [6.5,10.16] ms,
L5 tophat [6.5,10] ms. fig53."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).parent
RHO_J = 1.6413
exp_axis = np.load(HERE.parent / "panda_m142den_axis_EPAPS.npz")
CASES = [
    ("STATS3D_L1sens.npz",   "baseline L1 (δ=254 µm, 1588 µm)", "#4292c6", "-"),
    ("STATS3D_tophatL1.npz", "tophat L1 (δ=0, 1588 µm)",        "#fd8d3c", "-"),
    ("STATS3D_baseline.npz", "baseline L5 (δ=254 µm, 99 µm)",  "#08306b", "--"),
    ("STATS3D_tophat.npz",   "tophat L5 (δ=0, 99 µm)",         "#a63603", "--"),
]
fig, ax = plt.subplots(figsize=(11, 5.5))
for fn, lab, col, ls in CASES:
    f = HERE / fn
    if not f.exists():
        f = HERE.parent / fn
    d = np.load(f, allow_pickle=True)
    xx = d["xx"]
    rho = d["ax_DensityMEAN"]
    cl = rho[0].copy()
    for rk in (1, 2):
        gap = ~np.isfinite(cl)
        cl[gap] = rho[rk][gap]
    m = (xx > 0.3) & (xx <= 10)
    ax.plot(xx[m], cl[m] / RHO_J, color=col, lw=1.6, ls=ls, label=lab)
ax.plot(exp_axis["x"], exp_axis["rho"], "ko", ms=4, mfc="none",
        label="Panda & Seasholtz (1999)")
ax.set_xlabel("$x/D_e$"); ax.set_ylabel("$\\langle\\rho\\rangle/\\rho_j$")
ax.set_xlim(0.3, 10); ax.grid(alpha=0.25, lw=0.4)
ax.legend(fontsize=8)
ax.set_title("m142 3D — inflow-profile effect (δ) at L1 and L5, centerline density")
fig.tight_layout()
fig.savefig(HERE / "fig53_tophat_l1_compare.png", dpi=160)
print("wrote fig53_tophat_l1_compare.png")

# quantitative: delta effect (tophat-baseline) at each resolution, x/D in [1,6]
def cl_of(fn):
    f = HERE / fn
    if not f.exists():
        f = HERE.parent / fn
    d = np.load(f, allow_pickle=True)
    xx = d["xx"]; rho = d["ax_DensityMEAN"]
    cl = rho[0].copy()
    for rk in (1, 2):
        gap = ~np.isfinite(cl)
        cl[gap] = rho[rk][gap]
    return xx, cl / RHO_J
x1, b1 = cl_of("STATS3D_L1sens.npz"); _, t1 = cl_of("STATS3D_tophatL1.npz")
x5, b5 = cl_of("STATS3D_baseline.npz"); _, t5 = cl_of("STATS3D_tophat.npz")
m1 = (x1 >= 1) & (x1 <= 6); m5 = (x5 >= 1) & (x5 <= 6)
print(f"RMS(tophat-baseline) L1: {np.sqrt(np.nanmean((t1[m1]-b1[m1])**2)):.4f}")
print(f"RMS(tophat-baseline) L5: {np.sqrt(np.nanmean((t5[m5]-b5[m5])**2)):.4f}")
