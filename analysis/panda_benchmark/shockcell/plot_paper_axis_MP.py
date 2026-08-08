#!/usr/bin/env python3
"""Paper figures: centreline <M> (left) and <p>/p_inf (right).
Fig A: 2D levels l_max = 0..6.
Fig B: 3D levels l_max = 1..5; for each level the centreline plus the
axial lines at r/D = 0.02..0.10 (thin, same colour).
Consistent Blues level colours, panel tags, no titles."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams.update({
    "font.size": 9, "axes.linewidth": 0.6,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
})
HERE = Path(__file__).parent
D = 0.0254
P_AMB = 99780.0
GAM, RGAS = 1.4, 287.0
XLIM = (0.0, 10.0)

# ---------------- Fig A: 2D ----------------
C2 = plt.cm.Blues(np.linspace(0.30, 1.00, 7))
CASES_2D = [
    ("AXT_L0lo.npz", 0, "1588"), ("AXT_L1lo.npz", 1, "794"),
    ("AXT_L2lo.npz", 2, "397"),  ("AXT_L3.npz", 3, "198"),
    ("AXT_L4.npz", 4, "99"),     ("AXT_L5.npz", 5, "50"),
    ("AXT_L6.npz", 6, "25"),
]
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9))
for fn, lev, dxs in CASES_2D:
    d = np.load(HERE / fn, allow_pickle=True)
    x = d["z"] / D
    M = np.abs(d["stat_y_velocityMEAN"]) / np.sqrt(
        GAM * RGAS * np.maximum(d["stat_temperatureMEAN"], 1.0))
    p = d["stat_pressureMEAN"] / P_AMB
    m = (x > 0.05) & (x <= XLIM[1])
    lab = f"L{lev}"
    axs[0].plot(x[m], M[m], color=C2[lev], lw=1.1, label=lab)
    axs[1].plot(x[m], p[m], color=C2[lev], lw=1.1, label=lab)
axs[0].set_ylabel("$\\langle M\\rangle$")
axs[1].set_ylabel("$\\langle p\\rangle/p_\\infty$")
axs[0].axhline(1.0, color="k", lw=0.5, ls=":")
axs[1].axhline(1.0, color="k", lw=0.5, ls=":")
for k, ax in enumerate(axs):
    ax.set_xlim(*XLIM)
    ax.set_xlabel("$x/D_e$", labelpad=1.5)
    ax.tick_params(length=2.5)
    ax.text(0.02, 0.965, "(a)" if k == 0 else "(b)", transform=ax.transAxes,
            ha="left", va="top", fontsize=9)
axs[0].legend(fontsize=6.5, ncol=2, frameon=False, loc="upper right",
              handlelength=1.3, columnspacing=0.9, labelspacing=0.25)
fig.tight_layout(w_pad=1.4)
fig.savefig(HERE / "paperfig_axis2d_MP.png", dpi=300)
fig.savefig(HERE / "paperfig_axis2d_MP.pdf")
print("wrote paperfig_axis2d_MP")

# ---------------- Fig B: 3D with near-axis tube lines ----------------
C3 = plt.cm.Blues(np.linspace(0.35, 1.00, 5))
CASES_3D = [
    ("STATS3D_L1sens.npz", 1, "1588"), ("STATS3D_L2sens.npz", 2, "794"),
    ("STATS3D_L3sens.npz", 3, "397"),  ("STATS3D_L4sens.npz", 4, "198"),
    ("STATS3D_baseline.npz", 5, "99"),
]
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9))
for fn, lev, dxs in CASES_3D:
    f = HERE / fn
    if not f.exists():
        f = HERE.parent / fn
    d = np.load(f, allow_pickle=True)
    x = d["xx"]
    pM = d["ax_pressureMEAN"]          # (rings, NX)
    uM = d["ax_x_velocityMEAN"]
    rM = d["ax_DensityMEAN"]
    col = C3[lev - 1]
    m = (x > 0.05) & (x <= XLIM[1])
    nring = pM.shape[0]
    for rk in range(nring - 1, 0, -1):     # tube lines r = 0.02..0.10 D
        if not np.isfinite(pM[rk][m]).any():
            continue
        Mach = np.abs(uM[rk]) * np.sqrt(np.maximum(rM[rk], 1e-6)
                                        / (GAM * np.maximum(pM[rk], 1.0)))
        axs[0].plot(x[m], Mach[m], color=col, lw=0.45, alpha=0.35)
        axs[1].plot(x[m], (pM[rk] / P_AMB)[m], color=col, lw=0.45, alpha=0.35)
    cl_p = pM[0].copy(); cl_u = uM[0].copy(); cl_r = rM[0].copy()
    for rk in (1, 2):                       # centreline fallback (coarse levels)
        gap = ~np.isfinite(cl_p)
        cl_p[gap] = pM[rk][gap]; cl_u[gap] = uM[rk][gap]; cl_r[gap] = rM[rk][gap]
    Mach = np.abs(cl_u) * np.sqrt(np.maximum(cl_r, 1e-6)
                                  / (GAM * np.maximum(cl_p, 1.0)))
    lab = f"L{lev}"
    axs[0].plot(x[m], Mach[m], color=col, lw=1.2, label=lab)
    axs[1].plot(x[m], (cl_p / P_AMB)[m], color=col, lw=1.2, label=lab)
axs[0].set_ylabel("$\\langle M\\rangle$")
axs[1].set_ylabel("$\\langle p\\rangle/p_\\infty$")
axs[0].axhline(1.0, color="k", lw=0.5, ls=":")
axs[1].axhline(1.0, color="k", lw=0.5, ls=":")
for k, ax in enumerate(axs):
    ax.set_xlim(*XLIM)
    ax.set_xlabel("$x/D_e$", labelpad=1.5)
    ax.tick_params(length=2.5)
    ax.text(0.02, 0.965, "(a)" if k == 0 else "(b)", transform=ax.transAxes,
            ha="left", va="top", fontsize=9)
axs[0].legend(fontsize=6.5, ncol=2, frameon=False, loc="upper right",
              handlelength=1.3, columnspacing=0.9, labelspacing=0.25)
fig.tight_layout(w_pad=1.4)
fig.savefig(HERE / "paperfig_axis3d_MP.png", dpi=300)
fig.savefig(HERE / "paperfig_axis3d_MP.pdf")
print("wrote paperfig_axis3d_MP")
