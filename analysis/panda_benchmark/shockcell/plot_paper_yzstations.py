#!/usr/bin/env python3
"""Appendix-style yz cross-section sequences at four axial stations of the
3D L5 solution: x = x_fc/4, x_fc/2, x_fc and x_1 (x_fc = 0.876, x_1 = 1.178).
Fig 1 (paperfig_yzstations_mach): M / M_qbar / M_rms, definitions and colour
limits identical to Figure A.2.
Fig 2 (paperfig_yzstations_struct): numerical schlieren / |omega| De/Uj /
max(-lambda_{2,d},0) (De/Uj)^2, identical to Figure A.4; full velocity
gradients (d/dx from adjacent x-planes), per-level assembly."""
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
D = 0.0254; GAM, RGAS = 1.4, 287.0
UJ = 414.2                      # T0 = 297.15 K convention
S1 = D / UJ; S2 = (D / UJ)**2
STATIONS = [0.219, 0.438, 0.876, 1.178]
ROWLAB = ["$x_{fc}/4$", "$x_{fc}/2$", "$x_{fc}$", "$x_1$"]
R = 0.8

dm = np.load(HERE / "SIM_L5_yzstations_multi.npz", allow_pickle=True)
WANT = list(dm["want"])
def gv(st, name):
    return dm["x" + str(st).replace(".", "p")][WANT.index(name)].astype(np.float64)
ds = np.load(HERE / "SIM_L5_yzstations_struct.npz")
def sv(st, pre):
    return ds[pre + "_" + str(st).replace(".", "p")]

def draw(rows_data, outname, col_heads, cb_specs):
    """rows_data[row][col] = (field, vmin, vmax, cmap); cb_specs = list of
    (col_indices, label) for shared bottom colorbars."""
    NR, NC = len(rows_data), 3
    PW = 2.05; PH = PW
    ML, MR_, MT, MB = 0.50, 0.10, 0.26, 0.95
    WSP, HSP = 0.12, 0.12
    FW = ML + NC*PW + (NC-1)*WSP + MR_
    FH = MB + NR*PH + (NR-1)*HSP + MT
    fig = plt.figure(figsize=(FW, FH))
    ims = {}
    ny = rows_data[0][0][0].shape[0]
    ext = [-R, R, -R, R]
    for row in range(NR):
        for col in range(NC):
            F, v0, v1, cm = rows_data[row][col]
            left = (ML + col*(PW+WSP))/FW
            bottom = (MB + (NR-1-row)*(PH+HSP))/FH
            ax = fig.add_axes([left, bottom, PW/FW, PH/FH])
            n2 = F.shape[0]
            w = min(n2, int(round(2*R*D / (2*0.85*D/n2) )))
            c0 = (n2 - w)//2
            Fw = F[c0:c0+w, c0:c0+w]
            im = ax.imshow(Fw, origin="lower", extent=ext, vmin=v0, vmax=v1,
                           cmap=cm, aspect="equal", interpolation="nearest",
                           rasterized=True)
            ims[col] = im
            th = np.linspace(0, 2*np.pi, 200)
            lc = "w" if cm in ("magma", "viridis") else "k"
            ax.plot(0.5*np.cos(th), 0.5*np.sin(th), color=lc, lw=0.5, ls=":")
            tag = "(%s)" % chr(97 + NC*row + col)
            lab = tag + (" " + ROWLAB[row] if col == 0 else "")
            tcol = "k" if cm in ("gray", "gray_r") else "w"
            ax.text(0.035, 0.95, lab, transform=ax.transAxes, ha="left",
                    va="top", fontsize=8.5, color=tcol,
                    bbox=dict(facecolor="white" if tcol == "k" else "black",
                              edgecolor="none", alpha=0.45, pad=1.1))
            if row == 0:
                ax.text(0.5, 1.05, col_heads[col], transform=ax.transAxes,
                        ha="center", va="bottom", fontsize=9)
            ax.set_xlim(-R, R); ax.set_ylim(-R, R)
            ax.tick_params(length=2.5, labelsize=8)
            if row < NR-1:
                ax.set_xticklabels([])
            else:
                ax.set_xlabel("$z/D_e$", labelpad=1.5)
            if col == 0:
                ax.set_ylabel("$y/D_e$", labelpad=1.0)
            else:
                ax.set_yticklabels([])
    for cols, lab in cb_specs:
        x0 = (ML + cols[0]*(PW+WSP))/FW
        wdt = (len(cols)*PW + (len(cols)-1)*WSP)/FW
        cax = fig.add_axes([x0, 0.36/FH, wdt, 0.10/FH])
        cb = fig.colorbar(ims[cols[-1]], cax=cax, orientation="horizontal")
        cb.set_label(lab, fontsize=8.5, labelpad=2)
        cb.ax.tick_params(length=2.5, labelsize=8)
    fig.savefig(HERE / (outname + ".png"), dpi=300)
    fig.savefig(HERE / (outname + ".pdf"), dpi=300)
    print("wrote", outname)

# ---------- Fig 1: Mach ----------
def mach_inst(st):
    u2 = gv(st, "x_velocity")**2 + gv(st, "y_velocity")**2 + gv(st, "z_velocity")**2
    return np.sqrt(u2 / (GAM*RGAS*np.maximum(gv(st, "temperature"), 1.0)))
def mach_mean(st):
    u2 = (gv(st, "x_velocityMEAN")**2 + gv(st, "y_velocityMEAN")**2
          + gv(st, "z_velocityMEAN")**2)
    return np.sqrt(u2 / (GAM*RGAS*np.maximum(gv(st, "temperatureMEAN"), 1.0)))
def mach_rms(st):
    var = sum(np.maximum(gv(st, c + "SQR") - gv(st, c + "MEAN")**2, 0)
              for c in ("x_velocity", "y_velocity", "z_velocity"))
    return np.sqrt(var / (GAM*RGAS*np.maximum(gv(st, "temperatureMEAN"), 1.0)))

rows = []
for st in STATIONS:
    rows.append([(mach_inst(st), 0.0, 2.0, "viridis"),
                 (mach_mean(st), 0.0, 2.0, "viridis"),
                 (mach_rms(st), 0.0, 0.5, "magma")])
draw(rows, "paperfig_yzstations_mach",
     ["$M$", "$M_{\\bar{q}}$", "$M_{\\mathrm{rms}}$"],
     [((0, 1), "$M$, $M_{\\bar{q}}$"), ((2,), "$M_{\\mathrm{rms}}$")])

# ---------- Fig 2: structural ----------
gref = np.nanpercentile(np.concatenate(
    [sv(st, "grho").ravel() for st in STATIONS]), 99.9)
rows = []
for st in STATIONS:
    rows.append([(np.exp(-8.0*sv(st, "grho")/gref), 0.0, 1.0, "gray"),
                 (sv(st, "om")*S1, 0.0, 12.0, "viridis"),
                 (np.maximum(-sv(st, "lam2"), 0.0)*S2, 0.0, 30.0, "magma")])
draw(rows, "paperfig_yzstations_struct",
     ["numerical schlieren", "$|\\omega|$", "$-\\lambda_{2,d}$"],
     [((0,), "$\\exp(-8\\,|\\nabla\\rho|/|\\nabla\\rho|_{99.9})$"),
      ((1,), "$|\\omega|\\,D_e/U_j$"),
      ((2,), "$\\max(-\\lambda_{2,d},0)\\,(D_e/U_j)^2$")])
