#!/usr/bin/env python3
"""L4 (De/256) family: mean density fields + centreline, tophat WENO4 vs
flanged baseline profile, against the Panda M_j=1.42 experiment.

Cases:
  WENO4    s2d_llfweno_L4  : sharp top-hat inlet, LLF WENO-Z5 + 4th-order
           viscous, statistics window [2, 6] ms, eta0 = 3.2735.
  baseline SWEEP step 07   : flanged baseline profile (tanh 254 um), WENO-Z5
           + 2nd-order viscous, staircase quasi-steady window (t ~ 11.7 ms),
           eta0 = 3.273.
Both are flanged-plate geometries at matched eta0; they differ in the exit
velocity profile (and viscous order). Fields on the identical D/256 grid,
r in [0,3De], z in [0,10De].
"""
from pathlib import Path
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
D, RHO_J = 0.0254, 1.6413

w = np.load(HERE/"npr_sweep_case"/"WENO4_field_rmean.npz")
Fw = w["field"]/RHO_J
b = np.load(HERE/"npr_sweep_case"/"SWEEP_field_07_p3.273_stats.npz")
Fb = b["rmean"]/RHO_J
dx = float(w["dx"])
NR, NZ = Fw.shape
zz = (np.arange(NZ)+0.5)*dx/D
rr = (np.arange(NR)+0.5)*dx/D

# axis curves
fs = sorted(glob.glob(str(HERE/"llfweno2d_axis"/"AX2D_plt*.npz")))
dw = np.load(fs[-1])
xw, aw = dw["x"]/D, dw["ax_DensityMEAN"]/RHO_J
xb, ab = zz, Fb[0, :]
e = np.load(HERE/"panda_m142den_axis_EPAPS.npz")
xe, re = e["x"], e["rho"]

for name, xs, rs in [("tophat WENO4", xw, aw), ("baseline flanged", xb, ab)]:
    m = (xe >= 0.6) & (xe <= min(xs.max(), 6.92))
    si = np.interp(xe[m], xs, rs)
    dv = si - re[m]
    print("%-20s rms=%.4f bias=%+.4f corr=%.4f"
          % (name, np.sqrt(np.mean(dv*dv)), dv.mean(),
             np.corrcoef(si, re[m])[0, 1]))

plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "legend.frameon": False})
fig = plt.figure(figsize=(10.4, 8.0))
gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.35], width_ratios=[1, 0.02],
                      hspace=0.28, wspace=0.03,
                      left=0.075, right=0.93, top=0.95, bottom=0.075)
a1 = fig.add_subplot(gs[0, 0])
a2 = fig.add_subplot(gs[1, 0], sharex=a1, sharey=a1)
cax = fig.add_subplot(gs[0:2, 1])
a3 = fig.add_subplot(gs[2, 0])

vmin, vmax = 0.2, 1.5
for ax, F, lab in [(a1, Fw, "top-hat, LLF WENO-Z5 + 4th-order viscous, window [2, 6] ms"),
                   (a2, Fb, "baseline profile (flanged), WENO-Z5 + 2nd-order viscous, quasi-steady")]:
    im = ax.pcolormesh(zz, rr, F, cmap="viridis", vmin=vmin, vmax=vmax,
                       shading="auto", rasterized=True)
    ax.set_ylabel("$r/D_e$")
    ax.set_ylim(0, 3)
    ax.text(0.008, 0.93, lab, transform=ax.transAxes, va="top",
            color="w", fontsize=9)
a2.set_xlabel("")
a1.tick_params(labelbottom=False)
cb = fig.colorbar(im, cax=cax)
cb.set_label(r"$\overline{\rho}/\rho_j$")

a3.plot(xw, aw, color="#0b6ef5", lw=2.0, zorder=5,
        label="top-hat, WENO-Z5 + 4th-order viscous")
a3.plot(xb, ab, color="#8a2be2", lw=1.7, ls="--", zorder=4,
        label="baseline profile (flanged)")
a3.plot(xe, re, "o", ms=4.2, mfc="w", mec="k", mew=0.9, ls="", zorder=7,
        label="Panda & Seasholtz (1999), $M_j$ = 1.42")
a3.set_xlim(0, 10)
a3.set_ylim(0.2, 1.6)
a3.set_xlabel("$x/D_e$")
a3.set_ylabel(r"$\overline{\rho}/\rho_j$  on the axis")
a3.grid(True, ls="--", lw=0.5, color="0.92")
a3.set_axisbelow(True)
a3.legend(loc="lower right", fontsize=8.8, handlelength=2.2)
for ax, tag in [(a1, "(a)"), (a2, "(b)"), (a3, "(c)")]:
    ax.text(-0.06, 1.02, tag, transform=ax.transAxes, fontsize=10.5)

fig.savefig(HERE/"l4_weno4_vs_baseline.png", dpi=230)
fig.savefig(HERE/"l4_weno4_vs_baseline.pdf")
print("wrote l4_weno4_vs_baseline.png/pdf")
