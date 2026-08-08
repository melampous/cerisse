#!/usr/bin/env python3
"""m142 2D: does the quasi-steady L0 behavior survive doubling the base grid?
fig46: centerline <rho> and rho_rms — L0 (1584 um), L0hi (792 um uniform),
L1 (792 um via AMR), L6 (25 um reference) + experiment."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).parent
RHO_J = 1.6413
CASES = [
    ("L0lo",          "L0 (1584 µm, base 192×512)",        "#9ecae1", "-"),
    ("m142_2d_L0hi",  "L0hi (792 µm, base 384×1024)",      "#d94801", "-"),
    ("L1lo",          "L1 (792 µm = 1584 base + 1 lvl)",   "#4292c6", "--"),
    ("L6",            "L6 (25 µm) reference",              "#08306b", "-"),
]
exp_axis = np.load(HERE.parent / "panda_m142den_axis_EPAPS.npz")
fig, axs = plt.subplots(2, 1, figsize=(12, 8.5), sharex=True)
for key, lab, col, ls in CASES:
    d = np.load(HERE / f"AXT_{key}.npz", allow_pickle=True)
    z = d["z"] / 0.0254
    rm = d["stat_DensityMEAN"] / RHO_J
    rs = d["stat_DensitySQR"] / RHO_J**2
    rrms = np.sqrt(np.maximum(rs - rm**2, 0))
    m = (z > 0.2) & (z <= 10)
    axs[0].plot(z[m], rm[m], color=col, lw=1.6, ls=ls, label=lab)
    axs[1].plot(z[m], rrms[m], color=col, lw=1.4, ls=ls, label=lab)
axs[0].plot(exp_axis["x"], exp_axis["rho"], "ko", ms=4, mfc="none",
            label="Panda & Seasholtz (1999)")
axs[0].set_ylabel("$\\langle\\rho\\rangle/\\rho_j$")
axs[0].set_title("m142 2D RZ — base-grid doubling test (TRUE statistics, window 1.5 ms)")
axs[1].set_ylabel("$\\rho_{rms}/\\rho_j$")
axs[1].set_xlabel("$x/D_e$")
for a in axs:
    a.grid(alpha=0.25, lw=0.4); a.legend(fontsize=8)
fig.tight_layout()
fig.savefig(HERE / "fig46_l0hi_compare.png", dpi=160)
print("wrote fig46_l0hi_compare.png")
for key, lab, _, _ in CASES[:3]:
    d = np.load(HERE / f"AXT_{key}.npz", allow_pickle=True)
    z = d["z"] / 0.0254
    rm = d["stat_DensityMEAN"] / RHO_J
    rs = d["stat_DensitySQR"] / RHO_J**2
    prms = np.sqrt(np.maximum(rs - rm**2, 0))
    m = (z >= 1) & (z <= 6)
    print(f"{key:14s} median rho_rms [1,6]D = {np.nanmedian(prms[m]):.4f}")
