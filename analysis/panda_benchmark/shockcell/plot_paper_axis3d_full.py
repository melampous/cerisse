#!/usr/bin/env python3
"""Paper figure: 3D centreline statistics, l_max = 1-5, vs experiment.
(a) <rho>/rho_j with Panda & Seasholtz (1999) fig4 data (Rayleigh scattering;
    the only measured centreline quantity),
(b) <M>, (c) <p>/p_inf (simulation only),
(d) centreline p_rms/p_inf, box-filtered over 0.15 D to suppress residual
    sampling noise; markers denote each level's global maximum.
Thin lines in (a-c): axial profiles at r/D = 0.02-0.10 (near-axis tube)."""
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
RHO_J = 1.6413
GAM, RGAS = 1.4, 287.0
XLIM = (0.0, 10.0)
C3 = plt.cm.Blues(np.linspace(0.35, 1.00, 5))
RING = [
    ("STATS3D_L1sens.npz", "AX3T_L1.npz", 1),
    ("STATS3D_L2sens.npz", "AX3T_L2.npz", 2),
    ("STATS3D_L3sens.npz", "AX3T_L3.npz", 3),
    ("STATS3D_L4sens.npz", "AX3T_L4.npz", 4),
    ("STATS3D_baseline.npz", "AX3T_L5.npz", 5),
]
exp_axis = np.load(HERE.parent / "panda_m142den_axis_EPAPS.npz")

def boxfilt(a, x, wD=0.15):
    w = max(3, int(round(wD / (x[1] - x[0]))))
    if w % 2 == 0:
        w += 1
    k = np.ones(w) / w
    out = np.convolve(a, k, mode="same")
    out[:w//2] = a[:w//2]
    out[-(w//2):] = a[-(w//2):]
    return out

fig, axs = plt.subplots(2, 2, figsize=(7.6, 4.4), sharex=True)
axr, axm, axp, axf = axs[0, 0], axs[0, 1], axs[1, 0], axs[1, 1]
for fn_ring, fn_ax, lev in RING:
    col = C3[lev - 1]
    f = HERE / fn_ring
    if not f.exists():
        f = HERE.parent / fn_ring
    dr = np.load(f, allow_pickle=True)
    x = dr["xx"]
    m = (x > 0.05) & (x <= XLIM[1])
    pR, uR, rR = dr["ax_pressureMEAN"], dr["ax_x_velocityMEAN"], dr["ax_DensityMEAN"]
    nring = pR.shape[0]
    for rk in range(nring - 1, 0, -1):
        if not np.isfinite(pR[rk][m]).any():
            continue
        Mach = np.abs(uR[rk]) * np.sqrt(np.maximum(rR[rk], 1e-6)
                                        / (GAM * np.maximum(pR[rk], 1.0)))
        axr.plot(x[m], (rR[rk] / RHO_J)[m], color=col, lw=0.4, alpha=0.3)
        axm.plot(x[m], Mach[m], color=col, lw=0.4, alpha=0.3)
        axp.plot(x[m], (pR[rk] / P_AMB)[m], color=col, lw=0.4, alpha=0.3)
    cp, cu, cr = pR[0].copy(), uR[0].copy(), rR[0].copy()
    for rk in (1, 2):
        g = ~np.isfinite(cp)
        cp[g] = pR[rk][g]; cu[g] = uR[rk][g]; cr[g] = rR[rk][g]
    Mach = np.abs(cu) * np.sqrt(np.maximum(cr, 1e-6) / (GAM * np.maximum(cp, 1.0)))
    lab = f"L{lev}"
    axr.plot(x[m], (cr / RHO_J)[m], color=col, lw=1.1, label=lab)
    axm.plot(x[m], Mach[m], color=col, lw=1.1)
    axp.plot(x[m], (cp / P_AMB)[m], color=col, lw=1.1)
    # (d) filtered centreline p_rms from the axis MEAN+SQR extraction
    da = np.load(HERE / fn_ax, allow_pickle=True)
    xa = da["x"] / D
    pm = da["ax_pressureMEAN"]; ps = da["ax_pressureSQR"]
    prms = np.sqrt(np.maximum(ps - pm**2, 0)) / P_AMB
    ma = (xa > 0.05) & (xa <= XLIM[1])
    pf = boxfilt(prms[ma], xa[ma])
    axf.plot(xa[ma], pf, color=col, lw=1.1)
    j = np.nanargmax(pf)
    axf.plot(xa[ma][j], pf[j], "o", color=col, ms=3.5, mec="k", mew=0.4)
axr.plot(exp_axis["x"], exp_axis["rho"], "ko", ms=3, mfc="none", mew=0.7,
         label="Panda & Seasholtz (1999)")
axr.set_ylabel("$\\langle\\rho\\rangle/\\rho_j$", labelpad=1.5)
axm.set_ylabel("$\\langle M\\rangle$", labelpad=1.5)
axp.set_ylabel("$\\langle p\\rangle/p_\\infty$", labelpad=1.5)
axf.set_ylabel("$p_{\\mathrm{rms}}/p_\\infty$", labelpad=1.5)
axm.axhline(1.0, color="k", lw=0.5, ls=":")
axp.axhline(1.0, color="k", lw=0.5, ls=":")
for k, ax in enumerate(axs.flat):
    ax.set_xlim(*XLIM)
    ax.tick_params(length=2.5)
    ax.text(0.015, 0.94, f"({chr(97+k)})", transform=ax.transAxes,
            ha="left", va="top", fontsize=9)
for ax in axs[1, :]:
    ax.set_xlabel("$x/D_e$", labelpad=1.5)
axr.legend(fontsize=6.3, ncol=2, frameon=False, loc="upper right",
           handlelength=1.2, columnspacing=0.8, labelspacing=0.25)
fig.tight_layout(w_pad=1.2, h_pad=0.6)
fig.savefig(HERE / "paperfig_axis3d_full.png", dpi=300)
fig.savefig(HERE / "paperfig_axis3d_full.pdf")
print("wrote paperfig_axis3d_full")
