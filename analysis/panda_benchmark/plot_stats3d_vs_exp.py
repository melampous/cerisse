#!/usr/bin/env python3
"""3D time-MEAN vs Panda & Seasholtz (1999): centerline + radial scans.
Legacy methodology (RADIAL_STATIONS_8, ±0.30D offset band, rho/rho_j norm),
applied to the TRUE statistics windows:
  baseline: [6.5, 10.16] ms (H100 production)   tophat: [6.5, 10] ms (A100)
Inputs: STATS3D_<case>.npz from stats_compare3d_extract.py
Outputs: fig27_axis_vs_exp.png, fig28_radial_vs_exp.png"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent
RHO_J_SIM = 1.6413
STATIONS = [0.60, 0.75, 0.90, 1.05, 1.25, 1.55, 2.00, 3.05]

exp_axis = np.load(HERE / "panda_m142den_axis_EPAPS.npz")
exp5 = np.load(HERE / "panda_fig5a_mj144_EPAPS.npz")

def exp_radial(station):
    key = next(k for k in exp5.files
               if k.startswith("x") and abs(float(k[1:]) - station) < 1e-8)
    d = exp5[key]
    r = np.abs(d[:, 0]); rho = d[:, 1] / float(exp5["rhoj"])
    o = np.argsort(r)
    return r[o], rho[o]

CASES = []
for name, fn, col in [("baseline (WENO, δ=254µm) [6.5,10.16]ms", "STATS3D_baseline.npz", "#084594"),
                      ("tophat (δ=0) [6.5,10]ms",               "STATS3D_tophat.npz",   "#d94801"),
                      ("teno (TENO5, δ=254µm) [6.5,9.3]ms",     "STATS3D_teno.npz",     "#238b45"),
                      ("euler (inviscid, δ=254µm) [6.5,9.3]ms", "STATS3D_euler.npz",    "#756bb1")]:
    f = HERE / "shockcell" / fn
    if not f.exists():
        f = HERE / fn
    if f.exists():
        CASES.append((name, np.load(f), col))

# ---------------- fig27: centerline ----------------
fig, ax = plt.subplots(figsize=(10.5, 5))
for name, d, col in CASES:
    xx = d["xx"]
    rho = d["ax_DensityMEAN"]          # (rings, NXF)
    m = (xx > 0.3) & (xx <= 10)
    # centerline = innermost ring with data (coarse levels have no cell at r<0.01D)
    cl = rho[0].copy()
    for rk in (1, 2):
        gap = ~np.isfinite(cl)
        cl[gap] = rho[rk][gap]
    ax.plot(xx[m], cl[m] / RHO_J_SIM, color=col, lw=1.6, label=name)
    band = rho[1:6] / RHO_J_SIM        # rings r=0.02..0.10 D
    ax.fill_between(xx[m], np.nanmin(band, axis=0)[m], np.nanmax(band, axis=0)[m],
                    color=col, alpha=0.18, lw=0,
                    label=f"{name.split(' ')[0]} near-axis rings r≤0.10D")
ax.plot(exp_axis["x"], exp_axis["rho"], "ko", ms=4, mfc="none",
        label="Panda & Seasholtz (1999) centerline")
ax.set_xlabel("$x/D_e$"); ax.set_ylabel("$\\langle\\rho\\rangle/\\rho_j$")
ax.set_xlim(0.3, 10); ax.grid(alpha=0.25, lw=0.4)
ax.legend(fontsize=8, loc="upper right")
ax.set_title("m142 3D L5 — time-mean centerline density vs experiment")
fig.tight_layout()
fig.savefig(HERE / "fig27_axis_vs_exp.png", dpi=160)
print("wrote fig27_axis_vs_exp.png")

# ---------------- fig28: radial scans ----------------
fig, axs = plt.subplots(2, 4, figsize=(16, 7.5), sharey=True)
for i, st in enumerate(STATIONS):
    ax = axs.flat[i]
    for name, d, col in CASES:
        rr = d["rr"]
        prof = d["rad_DensityMEAN"][i] / RHO_J_SIM   # (noff, NR)
        mid = prof[15]                                # o+0.00
        lo = np.nanmin(prof, axis=0); hi = np.nanmax(prof, axis=0)
        m = rr <= 2.0
        ax.fill_between(rr[m], lo[m], hi[m], color=col, alpha=0.18, lw=0)
        ax.plot(rr[m], mid[m], color=col, lw=1.5,
                label=name.split(" [")[0] if i == 0 else None)
    er, ed = exp_radial(st)
    ax.plot(er, ed, "ko", ms=3.5, mfc="none",
            label="exp (Mj=1.44, fig5a)" if i == 0 else None)
    ax.set_title(f"$x/D_e$ = {st}", fontsize=10)
    ax.grid(alpha=0.25, lw=0.4)
    if i >= 4:
        ax.set_xlabel("$r/D_e$")
    if i % 4 == 0:
        ax.set_ylabel("$\\langle\\rho\\rangle/\\rho_j$")
fig.legend(loc="lower center", ncol=3, fontsize=9, frameon=False)
fig.suptitle("m142 3D L5 — time-mean radial density profiles vs experiment "
             "(band: station ±0.30$D_e$ axial offsets)", y=0.99)
fig.tight_layout(rect=[0, 0.05, 1, 1])
fig.savefig(HERE / "fig28_radial_vs_exp.png", dpi=160)
print("wrote fig28_radial_vs_exp.png")
