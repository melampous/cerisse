#!/usr/bin/env python3
"""3D grid sensitivity, preliminary: L4 (frozen L4 grid, TRUE stats [8.0,9.56]ms)
vs L5 baseline ([6.5,10.16]ms) centerline density + experiment.
Extends to L3/L2/L1 automatically as STATS3D_L{n}sens.npz appear."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).parent
RHO_J = 1.6413
exp_axis = np.load(HERE.parent / "panda_m142den_axis_EPAPS.npz")
CAND = [
    ("STATS3D_L1sens.npz", "L1 (1588 µm) sens",  "#c6dbef"),
    ("STATS3D_L2sens.npz", "L2 (794 µm) sens",  "#9ecae1"),
    ("STATS3D_L3sens.npz", "L3 (397 µm) sens",  "#6baed6"),
    ("STATS3D_L4sens.npz", "L4 (198 µm) [8.0,9.56]ms", "#2171b5"),
    ("STATS3D_baseline.npz", "L5 (99 µm) baseline [6.5,10.16]ms", "#08306b"),
]
fig, ax = plt.subplots(figsize=(11, 5.5))
for fn, lab, col in CAND:
    f = HERE / fn
    if not f.exists():
        f = HERE.parent / fn
    if not f.exists():
        continue
    d = np.load(f, allow_pickle=True)
    xx = d["xx"]
    rho = d["ax_DensityMEAN"]
    cl = rho[0].copy()
    for rk in (1, 2):
        gap = ~np.isfinite(cl)
        cl[gap] = rho[rk][gap]
    m = (xx > 0.3) & (xx <= 10)
    ax.plot(xx[m], cl[m] / RHO_J, color=col, lw=1.6, label=lab)
ax.plot(exp_axis["x"], exp_axis["rho"], "ko", ms=4, mfc="none",
        label="Panda & Seasholtz (1999)")
ax.set_xlabel("$x/D_e$"); ax.set_ylabel("$\\langle\\rho\\rangle/\\rho_j$")
ax.set_xlim(0.3, 10); ax.grid(alpha=0.25, lw=0.4)
ax.legend(fontsize=8)
ax.set_title("m142 3D — grid-sensitivity centerline density (TRUE statistics, frozen grids)")
fig.tight_layout()
fig.savefig(HERE / "fig47_3dsens_axis.png", dpi=160)
print("wrote fig47_3dsens_axis.png")
