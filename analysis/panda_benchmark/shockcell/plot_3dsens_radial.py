#!/usr/bin/env python3
"""3D grid sensitivity L1-L5 vs experiment: radial density profiles at the
legacy 8 stations (band = ±0.30D axial offsets). fig52."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from pathlib import Path

HERE = Path(__file__).parent
RHO_J = 1.6413
STATIONS = [0.60, 0.75, 0.90, 1.05, 1.25, 1.55, 2.00, 3.05]
BLUES = plt.cm.Blues(np.linspace(0.35, 1.00, 5))
CASES = [
    ("STATS3D_L1sens.npz",   "L1 (1588 µm) [8.0,10.4]ms",  BLUES[0]),
    ("STATS3D_L2sens.npz",   "L2 (794 µm) [8.0,10.4]ms",  BLUES[1]),
    ("STATS3D_L3sens.npz",   "L3 (397 µm) [8.0,10.4]ms",  BLUES[2]),
    ("STATS3D_L4sens.npz",   "L4 (198 µm) [8.0,9.56]ms",  BLUES[3]),
    ("STATS3D_baseline.npz", "L5 (99 µm) [6.5,10.16]ms",  BLUES[4]),
]
exp5 = np.load(HERE.parent / "panda_fig5a_mj144_EPAPS.npz")
def exp_radial(station):
    key = next(k for k in exp5.files
               if k.startswith("x") and abs(float(k[1:]) - station) < 1e-8)
    d = exp5[key]
    r = np.abs(d[:, 0]); rho = d[:, 1] / float(exp5["rhoj"])
    o = np.argsort(r)
    return r[o], rho[o]

STROKE = [pe.withStroke(linewidth=2.4, foreground="#40506080")]
fig, axs = plt.subplots(2, 4, figsize=(16, 7.5), sharey=True)
for i, st in enumerate(STATIONS):
    ax = axs.flat[i]
    for j, (fn, lab, col) in enumerate(CASES):
        f = HERE / fn
        if not f.exists():
            f = HERE.parent / fn
        d = np.load(f, allow_pickle=True)
        rr = d["rr"]
        prof = d["rad_DensityMEAN"][i] / RHO_J   # (noff, NR)
        mid = prof[15]
        m = rr <= 2.0
        eff = STROKE if j == 0 else None
        ax.plot(rr[m], mid[m], color=col, lw=1.5, path_effects=eff,
                label=lab if i == 0 else None)
    er, ed = exp_radial(st)
    ax.plot(er, ed, "ko", ms=3.5, mfc="none",
            label="Panda & Seasholtz (1999) fig5a" if i == 0 else None)
    ax.set_title(f"$x/D_e$ = {st}", fontsize=10)
    ax.grid(alpha=0.25, lw=0.4)
    if i >= 4:
        ax.set_xlabel("$r/D_e$")
    if i % 4 == 0:
        ax.set_ylabel("$\\langle\\rho\\rangle/\\rho_j$")
fig.legend(loc="lower center", ncol=6, fontsize=8.5, frameon=False)
fig.suptitle("m142 3D grid sensitivity L1-L5 — radial time-mean density vs experiment", y=0.99)
fig.tight_layout(rect=[0, 0.05, 1, 1])
fig.savefig(HERE / "fig52_3dsens_radial_vs_exp.png", dpi=160)
print("wrote fig52_3dsens_radial_vs_exp.png")
