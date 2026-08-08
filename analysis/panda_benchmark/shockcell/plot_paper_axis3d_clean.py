#!/usr/bin/env python3
"""Paper figure (clean 3-row stack): 3D centreline statistics, l_max = 1-5.
(a) <rho>/rho_j with Panda & Seasholtz (1999), (b) <M>, (c) filtered p_rms.
Display smoothing 0.05 D on mean curves (ring-quantisation removal) and
0.15 D on p_rms; stated in the caption."""
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
GAM = 1.4
XLIM = (0.0, 8.0)
C3 = plt.cm.Blues(np.linspace(0.45, 1.00, 5))
CASES = [
    ("STATS3D_L1sens.npz", "AX3T_L1.npz", 1),
    ("STATS3D_L2sens.npz", "AX3T_L2.npz", 2),
    ("STATS3D_L3sens.npz", "AX3T_L3.npz", 3),
    ("STATS3D_L4sens.npz", "AX3T_L4.npz", 4),
    ("STATS3D_baseline.npz", "AX3T_L5.npz", 5),
]
exp_axis = np.load(HERE.parent / "panda_m142den_axis_EPAPS.npz")

def boxfilt(a, x, wD):
    w = max(3, int(round(wD / (x[1] - x[0]))))
    if w % 2 == 0:
        w += 1
    out = np.convolve(a, np.ones(w) / w, mode="same")
    out[:w//2] = a[:w//2]; out[-(w//2):] = a[-(w//2):]
    return out

fig, axs = plt.subplots(3, 1, figsize=(7.0, 4.9), sharex=True)
for fn_ring, fn_ax, lev in CASES:
    col = C3[lev - 1]
    lw = 1.5 if lev == 5 else 1.0
    f = HERE / fn_ring
    if not f.exists():
        f = HERE.parent / fn_ring
    dr = np.load(f, allow_pickle=True)
    x = dr["xx"]
    m = (x > 0.05) & (x <= XLIM[1])
    pR, uR, rR = dr["ax_pressureMEAN"], dr["ax_x_velocityMEAN"], dr["ax_DensityMEAN"]
    cp, cu, cr = pR[0].copy(), uR[0].copy(), rR[0].copy()
    for rk in (1, 2):
        g = ~np.isfinite(cp)
        cp[g] = pR[rk][g]; cu[g] = uR[rk][g]; cr[g] = rR[rk][g]
    Mach = np.abs(cu) * np.sqrt(np.maximum(cr, 1e-6) / (GAM * np.maximum(cp, 1.0)))
    rho_s = boxfilt((cr / RHO_J)[m], x[m], 0.05)
    M_s = boxfilt(Mach[m], x[m], 0.05)
    # near-axis tube envelope (rings r/D = 0 .. 0.10) as shaded band
    nring = pR.shape[0]
    rho_band, M_band = [], []
    for rk in range(nring):
        if not np.isfinite(pR[rk][m]).any():
            continue
        Mk = np.abs(uR[rk]) * np.sqrt(np.maximum(rR[rk], 1e-6)
                                      / (GAM * np.maximum(pR[rk], 1.0)))
        rho_band.append(boxfilt((rR[rk] / RHO_J)[m], x[m], 0.05))
        M_band.append(boxfilt(Mk[m], x[m], 0.05))
    if len(rho_band) >= 2:
        rb = np.array(rho_band); Mb = np.array(M_band)
        axs[0].fill_between(x[m], np.nanmin(rb, 0), np.nanmax(rb, 0),
                            color=col, alpha=0.16, lw=0)
        axs[1].fill_between(x[m], np.nanmin(Mb, 0), np.nanmax(Mb, 0),
                            color=col, alpha=0.16, lw=0)
    lab = f"L{lev}"
    axs[0].plot(x[m], rho_s, color=col, lw=lw, label=lab)
    axs[1].plot(x[m], M_s, color=col, lw=lw)
    da = np.load(HERE / fn_ax, allow_pickle=True)
    xa = da["x"] / D
    pm = da["ax_pressureMEAN"]; ps = da["ax_pressureSQR"]
    prms = np.sqrt(np.maximum(ps - pm**2, 0)) / P_AMB
    ma = (xa > 0.05) & (xa <= XLIM[1])
    pf = boxfilt(prms[ma], xa[ma], 0.15)
    axs[2].plot(xa[ma], pf, color=col, lw=lw)
    j = np.nanargmax(pf)
    axs[2].plot(xa[ma][j], pf[j], "o", color=col, ms=3.5, mec="k", mew=0.4,
                zorder=5)
axs[0].plot(exp_axis["x"], exp_axis["rho"], "o", color="k", ms=2.8, mfc="none",
            mew=0.7, label="Panda & Seasholtz (1999)")
axs[0].set_ylabel("$\\langle\\rho\\rangle/\\rho_j$", labelpad=2)
axs[1].set_ylabel("$\\langle M\\rangle$", labelpad=2)
axs[2].set_ylabel("$p_{\\mathrm{rms}}/p_\\infty$", labelpad=2)
axs[1].axhline(1.0, color="k", lw=0.5, ls=":")
axs[0].set_ylim(0.3, 1.55)
axs[1].set_ylim(0.7, 2.5)
axs[2].set_yscale("function",
                  functions=(lambda y: np.square(y),
                             lambda y: np.sqrt(np.maximum(y, 0))))
axs[2].set_ylim(0.0, 0.27)
axs[2].set_yticks([0.05, 0.1, 0.15, 0.2, 0.25])
for k, ax in enumerate(axs):
    ax.set_xlim(*XLIM)
    ax.tick_params(length=2.5)
    ax.text(0.012, 0.92, f"({chr(97+k)})", transform=ax.transAxes,
            ha="left", va="top", fontsize=9)
axs[2].set_xlabel("$x/D_e$", labelpad=2)
axs[0].legend(fontsize=6.8, ncol=6, frameon=False, loc="upper right",
              handlelength=1.2, columnspacing=0.9, handletextpad=0.4,
              borderaxespad=0.2)
fig.align_ylabels(axs)
fig.tight_layout(h_pad=0.35)
fig.savefig(HERE / "paperfig_axis3d_clean.png", dpi=300)
fig.savefig(HERE / "paperfig_axis3d_clean.pdf")
print("wrote paperfig_axis3d_clean")
