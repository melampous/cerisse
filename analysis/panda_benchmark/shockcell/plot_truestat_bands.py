#!/usr/bin/env python3
"""2D RZ m142 — GENUINE per-timestep-accumulated statistics (frozen-grid reruns).
Axis <p> ± p_rms per refinement level + acceptance metric (inter-level mean
difference vs turbulent rms). Windows: L0/L1/L2/L3 [8.0,9.5], L4/L5 [6.5,8.0],
L6 [15.75,17.25] ms — all developed-flow, zero regrid during accumulation."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from pathlib import Path

HERE = Path(__file__).parent
P_AMB = 99780.0
BLUES = plt.cm.Blues(np.linspace(0.30, 1.00, 7))
CASES = [
    ("L0lo", "L0 (1584 µm) [8.0,9.5]ms",  BLUES[0]),
    ("L1lo", "L1 (792 µm) [8.0,9.5]ms",   BLUES[1]),
    ("L2lo", "L2 (396 µm) [8.0,9.5]ms",   BLUES[2]),
    ("L3",   "L3 (198 µm) [8.0,9.5]ms",   BLUES[3]),
    ("L4",   "L4 (99 µm) [6.5,8.0]ms",    BLUES[4]),
    ("L5",   "L5 (50 µm) [6.5,8.0]ms",    BLUES[5]),
    ("L6",   "L6 (25 µm) [15.75,17.25]ms", BLUES[6]),
]
STROKE = [pe.withStroke(linewidth=2.4, foreground="#40506080")]
fig, axs = plt.subplots(2, 1, figsize=(12, 8.5), sharex=True)
for i, (key, lab, col) in enumerate(CASES):
    d = np.load(HERE / f"AXT_{key}.npz", allow_pickle=True)
    z = d["z"] / 0.0254
    pm = d["stat_pressureMEAN"] / P_AMB
    ps = d["stat_pressureSQR"] / P_AMB**2
    prms = np.sqrt(np.maximum(ps - pm**2, 0))
    m = (z > 0.2) & (z <= 10)
    eff = STROKE if i < 2 else None
    axs[0].plot(z[m], pm[m], color=col, lw=1.6, label=lab, path_effects=eff)
    axs[0].fill_between(z[m], (pm - prms)[m], (pm + prms)[m], color=col, alpha=0.10, lw=0)
    axs[1].plot(z[m], prms[m], color=col, lw=1.4, label=lab, path_effects=eff)
axs[0].set_ylabel("$\\langle p\\rangle/p_\\infty$ (band: ±$p_{rms}$)")
axs[0].set_title("2D RZ m142 — TRUE per-step statistics (frozen grid), centerline pressure, L0–L6")
axs[1].set_ylabel("$p_{rms}/p_\\infty$")
axs[1].set_xlabel("$x/D_e$")
for a in axs:
    a.grid(alpha=0.25, lw=0.4); a.legend(fontsize=8)
fig.tight_layout()
fig.savefig(HERE / "fig39_truestat_bands.png", dpi=160)
print("wrote fig39_truestat_bands.png")

d6 = np.load(HERE / "AXT_L6.npz", allow_pickle=True)
z6 = d6["z"] / 0.0254
pm6 = d6["stat_pressureMEAN"] / P_AMB
print("\nlevel  RMS(mean_L-mean_L6)  median p_rms   ratio  (x/D in [1,6])")
for key, lab, _ in CASES[:-1]:
    d = np.load(HERE / f"AXT_{key}.npz", allow_pickle=True)
    z = d["z"] / 0.0254
    pm = d["stat_pressureMEAN"] / P_AMB
    ps = d["stat_pressureSQR"] / P_AMB**2
    prms = np.sqrt(np.maximum(ps - pm**2, 0))
    m = (z >= 1) & (z <= 6)
    dm = np.sqrt(np.nanmean((pm[m] - np.interp(z[m], z6, pm6))**2))
    mr = np.nanmedian(prms[m])
    print(f"{key:6s} {dm:12.4f} {mr:14.4f} {dm/mr:8.2f}")
# 检验: 真时均残余脉动应远低于之前动态网格版
for key in ["L2lo", "L4", "L5", "L6"]:
    d = np.load(HERE / f"AXT_{key}.npz", allow_pickle=True)
    z = d["z"] / 0.0254
    pm = d["stat_pressureMEAN"] / P_AMB
    li = d["last_p"] / P_AMB
    m = (z > 2) & (z < 10)
    k = 51
    tr = lambda a: a - np.convolve(a, np.ones(k)/k, mode="same")
    a, b = tr(pm[m]), tr(li[m])
    ok = np.isfinite(a) & np.isfinite(b)
    print(f"{key}: statMEAN残余脉动/单帧 = {np.std(a[ok])/max(np.std(b[ok]),1e-12):.3f} (动态网格版曾为0.30-0.38)")
