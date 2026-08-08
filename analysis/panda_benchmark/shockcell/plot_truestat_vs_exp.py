#!/usr/bin/env python3
"""2D RZ m142 TRUE statistics L0-L6 vs Panda & Seasholtz (1999) experiment.
fig40: centerline mean density (exp overlay) + density rms.
fig41: radial mean-density profiles at the legacy 8 stations vs exp fig5a
       (band = station +-0.30D axial offsets), needs RADT_*.npz."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from pathlib import Path

HERE = Path(__file__).parent
RHO_J_SIM = 1.6413
BLUES = plt.cm.Blues(np.linspace(0.30, 1.00, 7))
CASES = [
    ("L0lo", "L0 (1584 µm)", BLUES[0]),
    ("L1lo", "L1 (792 µm)",  BLUES[1]),
    ("L2lo", "L2 (396 µm)",  BLUES[2]),
    ("L3",   "L3 (198 µm)",  BLUES[3]),
    ("L4",   "L4 (99 µm)",   BLUES[4]),
    ("L5",   "L5 (50 µm)",   BLUES[5]),
    ("L6",   "L6 (25 µm)",   BLUES[6]),
]
STROKE = [pe.withStroke(linewidth=2.4, foreground="#40506080")]
exp_axis = np.load(HERE.parent / "panda_m142den_axis_EPAPS.npz")

# ---------------- fig40: centerline density ----------------
fig, axs = plt.subplots(2, 1, figsize=(12, 8.5), sharex=True)
for i, (key, lab, col) in enumerate(CASES):
    d = np.load(HERE / f"AXT_{key}.npz", allow_pickle=True)
    z = d["z"] / 0.0254
    rm = d["stat_DensityMEAN"] / RHO_J_SIM
    rs = d["stat_DensitySQR"] / RHO_J_SIM**2
    rrms = np.sqrt(np.maximum(rs - rm**2, 0))
    m = (z > 0.2) & (z <= 10)
    eff = STROKE if i < 2 else None
    axs[0].plot(z[m], rm[m], color=col, lw=1.6, label=lab, path_effects=eff)
    axs[1].plot(z[m], rrms[m], color=col, lw=1.4, label=lab, path_effects=eff)
axs[0].plot(exp_axis["x"], exp_axis["rho"], "ko", ms=4, mfc="none",
            label="Panda & Seasholtz (1999)")
axs[0].set_ylabel("$\\langle\\rho\\rangle/\\rho_j$")
axs[0].set_title("2D RZ m142 — TRUE per-step statistics, centerline density vs experiment, L0-L6")
axs[1].set_ylabel("$\\rho_{rms}/\\rho_j$")
axs[1].set_xlabel("$x/D_e$")
for a in axs:
    a.grid(alpha=0.25, lw=0.4); a.legend(fontsize=8, ncol=2)
fig.tight_layout()
fig.savefig(HERE / "fig40_truestat_axis_vs_exp.png", dpi=160)
print("wrote fig40_truestat_axis_vs_exp.png")

# ---------------- fig41: radial profiles ----------------
exp5 = np.load(HERE.parent / "panda_fig5a_mj144_EPAPS.npz")
def exp_radial(station):
    key = next(k for k in exp5.files
               if k.startswith("x") and abs(float(k[1:]) - station) < 1e-8)
    d = exp5[key]
    r = np.abs(d[:, 0]); rho = d[:, 1] / float(exp5["rhoj"])
    o = np.argsort(r)
    return r[o], rho[o]

have = [(key, lab, col) for key, lab, col in CASES
        if (HERE / f"RADT_{key}.npz").exists()]
if have:
    d0 = np.load(HERE / f"RADT_{have[0][0]}.npz", allow_pickle=True)
    STATIONS = d0["stations"]
    fig, axs = plt.subplots(2, 4, figsize=(16, 7.5), sharey=True)
    for i, st in enumerate(STATIONS):
        ax = axs.flat[i]
        for j, (key, lab, col) in enumerate(have):
            d = np.load(HERE / f"RADT_{key}.npz", allow_pickle=True)
            rr = d["rr"]
            prof = d["rad_DensityMEAN"][i] / RHO_J_SIM   # (noff, NR)
            mid = prof[15]
            lo = np.nanmin(prof, axis=0); hi = np.nanmax(prof, axis=0)
            m = rr <= 2.0
            ax.fill_between(rr[m], lo[m], hi[m], color=col, alpha=0.10, lw=0)
            eff = STROKE if j < 2 else None
            ax.plot(rr[m], mid[m], color=col, lw=1.4, path_effects=eff,
                    label=lab if i == 0 else None)
        er, ed = exp_radial(float(st))
        ax.plot(er, ed, "ko", ms=3.5, mfc="none",
                label="exp (Mj=1.44, fig5a)" if i == 0 else None)
        ax.set_title(f"$x/D_e$ = {st}", fontsize=10)
        ax.grid(alpha=0.25, lw=0.4)
        if i >= 4:
            ax.set_xlabel("$r/D_e$")
        if i % 4 == 0:
            ax.set_ylabel("$\\langle\\rho\\rangle/\\rho_j$")
    fig.legend(loc="lower center", ncol=5, fontsize=9, frameon=False)
    fig.suptitle("2D RZ m142 — TRUE-stat radial density profiles vs experiment "
                 "(band: station ±0.30$D_e$ axial offsets)", y=0.99)
    fig.tight_layout(rect=[0, 0.06, 1, 1])
    fig.savefig(HERE / "fig41_truestat_radial_vs_exp.png", dpi=160)
    print(f"wrote fig41_truestat_radial_vs_exp.png ({len(have)} levels)")
else:
    print("no RADT_*.npz yet — fig41 skipped")
