#!/usr/bin/env python3
"""Axis true-MEAN with ±RMS band, per refinement level (2D RZ m142).
Acceptance logic: if inter-level mean differences sit inside the turbulent RMS
band, the coarse level is statistically acceptable for position metrics.
Data: stat_* rows of AX2_*.npz (true running means at each case's window)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).parent
P_AMB = 99780.0
CASES = [
    ("L3dev",  "L3 (198 µm) [6.5,8.0]ms",  "#c6dbef"),
    ("L4",     "L4 (99 µm) [0,6.5]ms",     "#9ecae1"),
    ("L5full", "L5 (50 µm) [3.0,6.5]ms",   "#4292c6"),
    ("L6weno", "L6 (25 µm) [~,15.7]ms",    "#084594"),
]
fig, axs = plt.subplots(2, 1, figsize=(12, 8.5), sharex=True)
for key, lab, col in CASES:
    d = np.load(HERE / f"AX2_{key}.npz", allow_pickle=True)
    z = d["z"] / 0.0254
    pm = d["stat_pressureMEAN"] / P_AMB
    ps = d["stat_pressureSQR"] / P_AMB**2
    prms = np.sqrt(np.maximum(ps - pm**2, 0))
    m = (z > 0.2) & (z <= 10)
    axs[0].plot(z[m], pm[m], color=col, lw=1.5, label=lab)
    axs[0].fill_between(z[m], (pm - prms)[m], (pm + prms)[m], color=col, alpha=0.13, lw=0)
    axs[1].plot(z[m], prms[m], color=col, lw=1.4, label=lab)
axs[0].set_ylabel("$\\langle p\\rangle/p_\\infty$  (band: ±$p_{rms}$)")
axs[0].set_title("2D RZ m142 — centerline TRUE mean pressure ± RMS band, by refinement level")
axs[1].set_ylabel("$p_{rms}/p_\\infty$")
axs[1].set_xlabel("$x/D_e$")
for a in axs:
    a.grid(alpha=0.25, lw=0.4); a.legend(fontsize=8)
fig.tight_layout()
fig.savefig(HERE / "fig37_axis_mean_rms_bands.png", dpi=160)
print("wrote fig37_axis_mean_rms_bands.png")

# quantitative: inter-level mean diff vs rms, x in [1,6]D
d6 = np.load(HERE / "AX2_L6weno.npz", allow_pickle=True)
z6 = d6["z"] / 0.0254
pm6 = d6["stat_pressureMEAN"] / P_AMB
print("\nlevel  RMS(mean_L - mean_L6)   median p_rms   ratio (x/D in [1,6])")
for key, lab, _ in CASES[:-1]:
    d = np.load(HERE / f"AX2_{key}.npz", allow_pickle=True)
    z = d["z"] / 0.0254
    pm = d["stat_pressureMEAN"] / P_AMB
    ps = d["stat_pressureSQR"] / P_AMB**2
    prms = np.sqrt(np.maximum(ps - pm**2, 0))
    m = (z >= 1) & (z <= 6)
    pm6i = np.interp(z[m], z6, pm6)
    dm = np.sqrt(np.nanmean((pm[m] - pm6i)**2))
    mr = np.nanmedian(prms[m])
    print(f"{key:8s} {dm:10.4f} {mr:16.4f} {dm/mr:10.2f}")
