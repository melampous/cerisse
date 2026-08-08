#!/usr/bin/env python3
"""Paper figures: centreline grid-sensitivity, 4 stacked variables.
(a) M reconstructed from time-averaged primitives, (b) <rho>/rho_j with
Panda & Seasholtz (1999), (c) <u_x>/U_j, (d) rho_rms/rho_j (square-scaled
axis to emphasise peaks). 2D: l_max = 0-6; 3D: l_max = 1-5.
Level colours: turbo ramp (high adjacent-pair distinguishability),
finest level drawn heavier. Display smoothing 0.05 D (means) / 0.15 D (rms)."""
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
RHO_J = 1.6413
UJ = 414.2
GAM, RGAS = 1.4, 287.0
XLIM = (0.0, 10.0)
exp_axis = np.load(HERE.parent / "panda_m142den_axis_EPAPS.npz")

def boxfilt(a, x, wD):
    w = max(3, int(round(wD / (x[1] - x[0]))))
    if w % 2 == 0:
        w += 1
    out = np.convolve(a, np.ones(w) / w, mode="same")
    out[:w//2] = a[:w//2]; out[-(w//2):] = a[-(w//2):]
    return out

CASES2D = [("0", "AXT_L0lo.npz"), ("1", "AXT_L1lo.npz"), ("2", "AXT_L2lo.npz"),
           ("3", "AXT_L3.npz"), ("4", "AXT_L4.npz"), ("5", "AXT_L5.npz"),
           ("6", "AXT_L6.npz")]
CASES3D = [(t, "AX3T_L%s.npz" % t) for t in "12345"]

def load2d(fn):
    d = np.load(HERE / fn, allow_pickle=True)
    x = d["z"] / D
    T = d["stat_temperatureMEAN"]
    u = d["stat_y_velocityMEAN"]          # axial (2D RZ: y = axial)
    v = d["stat_x_velocityMEAN"]
    M = np.sqrt(u**2 + v**2) / np.sqrt(GAM * RGAS * np.maximum(T, 1.0))
    rho = d["stat_DensityMEAN"]
    rr = np.sqrt(np.maximum(d["stat_DensitySQR"] - rho**2, 0))
    return x, M, rho / RHO_J, np.abs(u) / UJ, rr / RHO_J

def load3d(fn):
    d = np.load(HERE / fn, allow_pickle=True)
    x = d["x"] / D
    T = d["ax_temperatureMEAN"]
    u = d["ax_x_velocityMEAN"]
    M = np.abs(u) / np.sqrt(GAM * RGAS * np.maximum(T, 1.0))
    rho = d["ax_DensityMEAN"]
    rr = np.sqrt(np.maximum(d["ax_DensitySQR"] - rho**2, 0))
    return x, M, rho / RHO_J, np.abs(u) / UJ, rr / RHO_J

def draw(cases, loader, outname):
    N = len(cases)
    cols = plt.cm.turbo(np.linspace(0.06, 0.94, N))
    fig, axs = plt.subplots(4, 1, figsize=(7.0, 6.4), sharex=True)
    for k, (ltag, fn) in enumerate(cases):
        col = cols[k]
        lw = 1.6 if k == N - 1 else 1.0
        x, M, rho, u, rr = loader(fn)
        m = (x > 0.05) & (x <= XLIM[1])
        lab = "L%s" % ltag
        axs[0].plot(x[m], boxfilt(M, x, 0.05)[m], color=col, lw=lw, label=lab)
        axs[1].plot(x[m], boxfilt(rho, x, 0.05)[m], color=col, lw=lw)
        axs[2].plot(x[m], boxfilt(u, x, 0.05)[m], color=col, lw=lw)
        axs[3].plot(x[m], boxfilt(rr, x, 0.15)[m], color=col, lw=lw)
    axs[1].plot(exp_axis["x"], exp_axis["rho"], "o", color="k", ms=2.8,
                mfc="none", mew=0.7, label="Panda & Seasholtz (1999)")
    axs[0].axhline(1.0, color="k", lw=0.5, ls=":")
    axs[0].set_ylabel("$M_{\\bar{q}}$", labelpad=2)
    axs[1].set_ylabel("$\\langle\\rho\\rangle/\\rho_j$", labelpad=2)
    axs[2].set_ylabel("$\\langle u_x\\rangle/U_j$", labelpad=2)
    axs[3].set_ylabel("$\\rho_{\\mathrm{rms}}/\\rho_j$", labelpad=2)
    axs[3].set_yscale("function",
                      functions=(lambda y: np.square(y),
                                 lambda y: np.sqrt(np.maximum(y, 0))))
    ymax = max(0.12, 1.06 * max(np.nanmax(l.get_ydata()) for l in axs[3].lines))
    axs[3].set_ylim(0.0, ymax)
    axs[3].set_yticks([t for t in (0.05, 0.1, 0.15, 0.2, 0.25, 0.3)
                       if t < ymax])
    for k, ax in enumerate(axs):
        ax.set_xlim(*XLIM)
        ax.tick_params(length=2.5)
        ax.text(0.012, 0.92, "(%s)" % chr(97 + k), transform=ax.transAxes,
                ha="left", va="top", fontsize=9)
    axs[3].set_xlabel("$x/D_e$", labelpad=2)
    h0, l0 = axs[0].get_legend_handles_labels()
    h1, l1 = axs[1].get_legend_handles_labels()
    axs[0].legend(h0 + h1[-1:], l0 + l1[-1:], fontsize=6.8, ncol=4,
                  frameon=False, loc="upper right", handlelength=1.2,
                  columnspacing=0.9, handletextpad=0.4, borderaxespad=0.2)
    fig.align_ylabels(axs)
    fig.tight_layout(h_pad=0.35)
    fig.savefig(HERE / (outname + ".png"), dpi=300)
    fig.savefig(HERE / (outname + ".pdf"))
    print("wrote", outname)

draw(CASES2D, load2d, "paperfig_axis2d_4var")
draw(CASES3D, load3d, "paperfig_axis3d_4var")
