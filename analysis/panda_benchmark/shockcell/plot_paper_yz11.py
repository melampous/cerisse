#!/usr/bin/env python3
"""Dense yz cross-section sequences of the 3D L5 solution at eleven axial
stations: k/10 * x_fc for k = 1..10 (x_fc = 0.876) plus x_1 = 1.178.
Figures:
  paperfig_yz11_pressure    : p/p_inf | <p>/p_inf | p_rms/p_inf
  paperfig_yz11_struct      : instantaneous numerical schlieren | |omega| |
                              max(-lambda_{2,d},0)   (as Figure A.4)
  paperfig_yz11_struct_mean : same operators applied to the MEAN fields:
                              |grad<rho>|, |curl<u>| (= |<omega>|), and
                              lambda_{2,d} of the mean-velocity gradient
                              (NOT <lambda_{2,d}>)
  paperfig_yz11_struct_std  : sample std over the N=9 statistics-window
                              snapshots of |grad rho|, |omega|, lambda_{2,d}
Full velocity gradients (d/dx from adjacent x-planes), per-level assembly.
U_j = 414.2 m/s (T0 = 297.15 K convention)."""
import sys
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
P_AMB = 99780.0
UJ = 414.2
S1 = D / UJ; S2 = (D / UJ)**2
STS = ["0.0876", "0.1752", "0.2628", "0.3504", "0.438", "0.5256",
       "0.6132", "0.7008", "0.7884", "0.876", "1.178"]
ROWLAB = ["$0.1\\,x_{fc}$", "$0.2\\,x_{fc}$", "$0.3\\,x_{fc}$",
          "$0.4\\,x_{fc}$", "$0.5\\,x_{fc}$", "$0.6\\,x_{fc}$",
          "$0.7\\,x_{fc}$", "$0.8\\,x_{fc}$", "$0.9\\,x_{fc}$",
          "$x_{fc}$", "$x_1$"]
R = 0.8

def key(st):
    return st.replace(".", "p")

def tag(i):
    return chr(97 + i) if i < 26 else "a" + chr(97 + i - 26)

def draw(rows_data, outname, col_heads, cb_specs):
    NR, NC = len(rows_data), 3
    PW = 1.80; PH = PW
    ML, MR_, MT, MB = 0.50, 0.10, 0.26, 0.92
    WSP, HSP = 0.10, 0.10
    FW = ML + NC*PW + (NC-1)*WSP + MR_
    FH = MB + NR*PH + (NR-1)*HSP + MT
    fig = plt.figure(figsize=(FW, FH))
    ims = {}
    for row in range(NR):
        for col in range(NC):
            F, v0, v1, cm = rows_data[row][col]
            left = (ML + col*(PW+WSP))/FW
            bottom = (MB + (NR-1-row)*(PH+HSP))/FH
            ax = fig.add_axes([left, bottom, PW/FW, PH/FH])
            n2 = F.shape[0]
            w = int(round(n2 * R / 0.85))
            c0 = (n2 - w)//2
            Fw = F[c0:c0+w, c0:c0+w]
            im = ax.imshow(Fw, origin="lower", extent=[-R, R, -R, R],
                           vmin=v0, vmax=v1, cmap=cm, aspect="equal",
                           interpolation="nearest", rasterized=True)
            ims[col] = im
            th = np.linspace(0, 2*np.pi, 200)
            lc = "k" if cm in ("gray", "gray_r") else "w"
            ax.plot(0.5*np.cos(th), 0.5*np.sin(th), color=lc, lw=0.5, ls=":")
            tg = "(%s)" % tag(NC*row + col)
            lab = tg + (" " + ROWLAB[row] if col == 0 else "")
            tcol = "k" if cm in ("gray", "gray_r") else "w"
            ax.text(0.04, 0.95, lab, transform=ax.transAxes, ha="left",
                    va="top", fontsize=8, color=tcol,
                    bbox=dict(facecolor="white" if tcol == "k" else "black",
                              edgecolor="none", alpha=0.45, pad=1.0))
            if row == 0:
                ax.text(0.5, 1.06, col_heads[col], transform=ax.transAxes,
                        ha="center", va="bottom", fontsize=9)
            ax.set_xlim(-R, R); ax.set_ylim(-R, R)
            ax.tick_params(length=2.2, labelsize=7.5)
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
        cax = fig.add_axes([x0, 0.34/FH, wdt, 0.09/FH])
        cb = fig.colorbar(ims[cols[-1]], cax=cax, orientation="horizontal")
        cb.set_label(lab, fontsize=8.5, labelpad=2)
        cb.ax.tick_params(length=2.5, labelsize=8)
    fig.savefig(HERE / (outname + ".png"), dpi=300)
    fig.savefig(HERE / (outname + ".pdf"), dpi=300)
    print("wrote", outname)

what = sys.argv[1:] if len(sys.argv) > 1 else ["pressure", "inst", "mean"]

if "pressure" in what:
    dm = np.load(HERE / "SIM_L5_yz11_multi.npz", allow_pickle=True)
    WANT = list(dm["want"])
    def gv(st, name):
        return dm["x" + key(st)][WANT.index(name)].astype(np.float64)
    prms_all = [np.sqrt(np.maximum(gv(st, "pressureSQR")
                                   - gv(st, "pressureMEAN")**2, 0))/P_AMB
                for st in STS]
    VP = float(np.ceil(np.nanpercentile(np.concatenate(
        [a.ravel() for a in prms_all]), 99.5) * 20) / 20)
    rows = []
    for st, pr in zip(STS, prms_all):
        rows.append([(gv(st, "pressure")/P_AMB, 0.4, 1.6, "RdBu_r"),
                     (gv(st, "pressureMEAN")/P_AMB, 0.4, 1.6, "RdBu_r"),
                     (pr, 0.0, VP, "magma")])
    draw(rows, "paperfig_yz11_pressure",
         ["$p/p_\\infty$", "$\\langle p\\rangle/p_\\infty$",
          "$p_{\\mathrm{rms}}/p_\\infty$"],
         [((0, 1), "$p/p_\\infty$, $\\langle p\\rangle/p_\\infty$"),
          ((2,), "$p_{\\mathrm{rms}}/p_\\infty$")])

def struct_rows(ds, pref=""):
    gref = np.nanpercentile(np.concatenate(
        [ds[pref + "grho_" + key(st)].ravel() for st in STS]), 99.9)
    rows = []
    for st in STS:
        rows.append([(np.exp(-8.0*ds[pref + "grho_" + key(st)]/gref),
                      0.0, 1.0, "gray"),
                     (ds[pref + "om_" + key(st)]*S1, 0.0, 12.0, "viridis"),
                     (np.maximum(-ds[pref + "lam2_" + key(st)], 0.0)*S2,
                      0.0, 30.0, "magma")])
    return rows

if "inst" in what:
    ds = np.load(HERE / "SIM_L5_yz11_struct.npz")
    draw(struct_rows(ds), "paperfig_yz11_struct",
         ["numerical schlieren", "$|\\omega|$", "$-\\lambda_{2,d}$"],
         [((0,), "$\\exp(-8\\,|\\nabla\\rho|/|\\nabla\\rho|_{99.9})$"),
          ((1,), "$|\\omega|\\,D_e/U_j$"),
          ((2,), "$\\max(-\\lambda_{2,d},0)\\,(D_e/U_j)^2$")])

if "mean" in what:
    ds = np.load(HERE / "SIM_L5_yz11_struct_MEAN.npz")
    draw(struct_rows(ds), "paperfig_yz11_struct_mean",
         ["numerical schlieren $\\langle\\rho\\rangle$",
          "$|\\nabla\\times\\langle u\\rangle|$",
          "$-\\lambda_{2,d}(\\langle u\\rangle)$"],
         [((0,), "$\\exp(-8\\,|\\nabla\\langle\\rho\\rangle|/|\\nabla\\langle\\rho\\rangle|_{99.9})$"),
          ((1,), "$|\\nabla\\times\\langle u\\rangle|\\,D_e/U_j$"),
          ((2,), "$\\max(-\\lambda_{2,d}(\\langle u\\rangle),0)\\,(D_e/U_j)^2$")])

if "std" in what:
    ds = np.load(HERE / "SIM_L5_yz11_struct_STD.npz")
    N = int(ds["N"])
    def sk(pre, st):
        return ds["std_" + pre + "_" + key(st)]
    vg = float(np.ceil(np.nanpercentile(np.concatenate(
        [sk("grho", st).ravel() for st in STS]), 99.5)))
    vo = float(np.ceil(np.nanpercentile(np.concatenate(
        [(sk("om", st)*S1).ravel() for st in STS]), 99.5)))
    vl = float(np.ceil(np.nanpercentile(np.concatenate(
        [(sk("lam2", st)*S2).ravel() for st in STS]), 99.5)))
    rows = []
    for st in STS:
        rows.append([(sk("grho", st), 0.0, vg, "magma"),
                     (sk("om", st)*S1, 0.0, vo, "magma"),
                     (sk("lam2", st)*S2, 0.0, vl, "magma")])
    draw(rows, "paperfig_yz11_struct_std",
         ["$\\sigma_{|\\nabla\\rho|}$ ($N=%d$)" % N,
          "$\\sigma_{|\\omega|}$", "$\\sigma_{\\lambda_{2,d}}$"],
         [((0,), "$\\sigma_{|\\nabla\\rho|}$ [kg m$^{-4}$]"),
          ((1,), "$\\sigma_{|\\omega|}\\,D_e/U_j$"),
          ((2,), "$\\sigma_{\\lambda_{2,d}}\\,(D_e/U_j)^2$")])
