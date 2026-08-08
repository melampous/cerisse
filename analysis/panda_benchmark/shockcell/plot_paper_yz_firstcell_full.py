#!/usr/bin/env python3
"""Paper figures: yz cross-sections inside the first shock cell
(x/D = 0.2, 0.4, 0.6, 0.8), baseline 3D l_max = 5.
Three figures, each 3 rows x 4 columns:
  Mach     : M (inst) / <M> / M_rms
  pressure : p (inst) / <p> / p_rms   (normalised by p_inf)
  schlieren: in-plane |grad rho| (inst) / |grad <rho>| / rho_rms
Shared bottom colorbars: rows 1-2 together, row 3 separate."""
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
d = np.load(HERE / "SIM_L5_yz4_multi.npz", allow_pickle=True)
WANT = list(d["want"])
def gv(st, name):
    return d["x" + str(st).replace(".", "p")][WANT.index(name)].astype(np.float64)
y, z = d["y"], d["z"]
dy = float(d["dy"]); dz = float(d["dz"])
STATIONS = [0.2, 0.4, 0.6, 0.8]
R = 0.8

def draw(fig_rows, outname, cb1_label, cb2_label):
    PW = 2.30; PH = PW
    ML, MR_, MT, MB = 0.50, 0.10, 0.22, 0.88
    WSP, HSP = 0.10, 0.10
    FW = ML + 4*PW + 3*WSP + MR_
    FH = MB + 3*PH + 2*HSP + MT
    fig = plt.figure(figsize=(FW, FH))
    ims = [None, None]
    for row, (fields, v0, v1, cm, rowlab) in enumerate(fig_rows):
        for col, st in enumerate(STATIONS):
            left = (ML + col*(PW+WSP))/FW
            bottom = (MB + (2-row)*(PH+HSP))/FH
            ax = fig.add_axes([left, bottom, PW/FW, PH/FH])
            F = fields[col]
            vv1 = v1 if v1 is not None else np.nanpercentile(F, 99.5)
            im = ax.imshow(F, origin="lower", extent=[z[0], z[-1], y[0], y[-1]],
                           vmin=v0, vmax=vv1, cmap=cm, aspect="equal",
                           interpolation="nearest", rasterized=True)
            if row < 2:
                ims[0] = im
            else:
                ims[1] = im
            th = np.linspace(0, 2*np.pi, 200)
            lc = "w" if cm in ("magma", "gray_r") else "k"
            ax.plot(0.5*np.cos(th), 0.5*np.sin(th), color=lc, lw=0.5, ls=":")
            ax.set_xlim(-R, R); ax.set_ylim(-R, R)
            tag = f"({chr(97 + 4*row + col)})"
            lab = tag + (f" {rowlab}" if col == 0 else "")
            tcol = "k" if cm == "gray_r" else "w"
            ax.text(0.03, 0.955, lab, transform=ax.transAxes, ha="left",
                    va="top", fontsize=8, color=tcol,
                    bbox=dict(facecolor="white" if tcol == "k" else "black",
                              edgecolor="none", alpha=0.45, pad=1.1))
            if row == 0:
                ax.text(0.5, 1.05, f"$x/D_e={st}$", transform=ax.transAxes,
                        ha="center", va="bottom", fontsize=9)
            ax.tick_params(length=2.5, labelsize=8)
            if row < 2:
                ax.set_xticklabels([])
            else:
                ax.set_xlabel("$z/D_e$", labelpad=1.5)
            if col == 0:
                ax.set_ylabel("$y/D_e$", labelpad=1.0)
            else:
                ax.set_yticklabels([])
    cax1 = fig.add_axes([ML/FW, 0.34/FH, (2.6*PW)/FW, 0.10/FH])
    cb1 = fig.colorbar(ims[0], cax=cax1, orientation="horizontal")
    cb1.set_label(cb1_label, fontsize=8.5, labelpad=2)
    cb1.ax.tick_params(length=2.5, labelsize=8)
    cax2 = fig.add_axes([(ML + 3*(PW+WSP))/FW, 0.34/FH, PW/FW, 0.10/FH])
    cb2 = fig.colorbar(ims[1], cax=cax2, orientation="horizontal")
    cb2.set_label(cb2_label, fontsize=8.5, labelpad=2)
    cb2.ax.tick_params(length=2.5, labelsize=8)
    fig.savefig(HERE / (outname + ".png"), dpi=300)
    fig.savefig(HERE / (outname + ".pdf"), dpi=300)
    print("wrote", outname)

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
def grad_mag(F):
    gy, gz = np.gradient(F, dy, dz)
    return np.sqrt(gy**2 + gz**2)

draw([
    ([mach_inst(s) for s in STATIONS], 0.0, 2.4, "viridis", "$M$"),
    ([mach_mean(s) for s in STATIONS], 0.0, 2.4, "viridis", "$\\langle M\\rangle$"),
    ([mach_rms(s) for s in STATIONS], 0.0, 0.5, "magma", "$M_{\\mathrm{rms}}$"),
], "paperfig_yz_firstcell_mach", "$M$, $\\langle M\\rangle$", "$M_{\\mathrm{rms}}$")

draw([
    ([gv(s, "pressure")/P_AMB for s in STATIONS], 0.2, 1.8, "RdYlBu_r", "$p/p_\\infty$"),
    ([gv(s, "pressureMEAN")/P_AMB for s in STATIONS], 0.2, 1.8, "RdYlBu_r",
     "$\\langle p\\rangle/p_\\infty$"),
    ([np.sqrt(np.maximum(gv(s, "pressureSQR")-gv(s, "pressureMEAN")**2, 0))/P_AMB
      for s in STATIONS], 0.0, None, "magma", "$p_{\\mathrm{rms}}/p_\\infty$"),
], "paperfig_yz_firstcell_pressure", "$p/p_\\infty$, $\\langle p\\rangle/p_\\infty$",
   "$p_{\\mathrm{rms}}/p_\\infty$")

draw([
    ([grad_mag(gv(s, "Density")) for s in STATIONS], 0.0, None, "gray_r",
     "$|\\nabla_{yz}\\rho|$"),
    ([grad_mag(gv(s, "DensityMEAN")) for s in STATIONS], 0.0, None, "gray_r",
     "$|\\nabla_{yz}\\langle\\rho\\rangle|$"),
    ([np.sqrt(np.maximum(gv(s, "DensitySQR")-gv(s, "DensityMEAN")**2, 0))/RHO_J
      for s in STATIONS], 0.0, None, "magma", "$\\rho_{\\mathrm{rms}}/\\rho_j$"),
], "paperfig_yz_firstcell_schlieren", "$|\\nabla_{yz}\\rho|$ [kg m$^{-4}$]",
   "$\\rho_{\\mathrm{rms}}/\\rho_j$")
