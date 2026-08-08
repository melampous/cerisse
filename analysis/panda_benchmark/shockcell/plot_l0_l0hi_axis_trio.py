#!/usr/bin/env python3
"""Three jets, centerline density: 2D L0 / L0hi TRUE statistics vs Panda &
Seasholtz (1999) fig4 EPAPS raw data (each block normalized by its own rhoj).
Output: fig51_axis_trio_l0_l0hi_vs_exp.png"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).parent
exp = np.load(HERE.parent / "panda_fig4_axis_allM_EPAPS.npz")
JOBS = [
    ("m119 (Mj=1.19)", "M1p19", 1.4873,
     [("AXT_m119_2d_L0.npz",  "L0 (1584 µm) [8.0,9.5]ms",       "#4292c6"),
      ("AXT_m119_2d_L0hi.npz","L0hi (792 µm, base×2) [8.0,9.5]ms","#d94801")]),
    ("m142 (Mj=1.44)", "M1p43", 1.6413,
     [("AXT_L0lo.npz",         "L0 (1584 µm) [8.0,9.5]ms",       "#4292c6"),
      ("AXT_m142_2d_L0hi.npz", "L0hi (792 µm, base×2) [8.0,9.5]ms","#d94801")]),
    ("m180 (Mj=1.80)", "M1p80", 1.9100,
     [("AXT_m180_2d_L0.npz",  "L0 (1584 µm) [11.5,13.0]ms",     "#4292c6"),
      ("AXT_m180_2d_L0hi.npz","L0hi (792 µm, base×2) [8.0,9.5]ms","#d94801")]),
]
fig, axs = plt.subplots(3, 1, figsize=(12, 12), sharex=True)
for ax, (title, ekey, rhoj, cases) in zip(axs, JOBS):
    for fn, lab, col in cases:
        d = np.load(HERE / fn, allow_pickle=True)
        z = d["z"] / 0.0254
        rm = d["stat_DensityMEAN"] / rhoj
        m = (z > 0.2) & (z <= 10)
        ax.plot(z[m], rm[m], color=col, lw=1.6, label=lab)
    ax.plot(exp[f"{ekey}_x"], exp[f"{ekey}_rho"], "ko", ms=3.5, mfc="none",
            label=f"Panda & Seasholtz (1999) fig4, $\\rho_j$={float(exp[f'{ekey}_rhoj']):.2f}")
    ax.set_title(title, fontsize=11)
    ax.set_ylabel("$\\langle\\rho\\rangle/\\rho_j$")
    ax.grid(alpha=0.25, lw=0.4)
    ax.legend(fontsize=8, loc="upper right")
axs[-1].set_xlabel("$x/D_e$")
fig.suptitle("2D RZ quasi-steady L0 vs L0hi — centerline density vs experiment", y=0.995)
fig.tight_layout(rect=[0, 0, 1, 0.985])
fig.savefig(HERE / "fig51_axis_trio_l0_l0hi_vs_exp.png", dpi=160)
print("wrote fig51_axis_trio_l0_l0hi_vs_exp.png")
