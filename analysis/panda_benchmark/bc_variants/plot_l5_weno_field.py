#!/usr/bin/env python3
"""AFD HLLC + WENO-Z5 on the L5 grid: centreline and mid-plane mean fields.

Top     the time-averaged centreline density, against the same scheme on L3,
        the existing 3D baseline (WENO-Z5) and the Panda & Seasholtz
        measurement.
Middle  mid-plane mean density: this run above the axis, the 3D baseline
        (WENO-Z5, same domain, same base mesh, also five levels) below it.
Bottom  the same split for the mean axial velocity.

Statistics are the running mean over 15.0 - 17.6 ms, ten frames. The slice is
the k plane through z = 0 of the last frame; every level from 0 to 5 is
composited, so the picture carries the true resolution everywhere.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.signal import find_peaks

HERE = Path(__file__).resolve().parent
OUT = HERE/"sweep_out"
PB = HERE.parent
RHO_J = 1.6413
EXP = dict(x_fc=0.888, S_1=1.359)


def refine(x, y, i):
    if i <= 0 or i >= len(y)-1:
        return x[i]
    d = y[i-1]-2*y[i]+y[i+1]
    return x[i] if d == 0 else x[i]+0.5*(y[i-1]-y[i+1])/d*(x[1]-x[0])


def xfc_of(x, r):
    m = np.isfinite(r)
    e = np.where((x >= 0.35) & (x <= 1.0) & m)[0]
    i0 = e[np.nanargmin(r[e])]
    c = np.where((np.arange(x.size) > i0) & (x <= 1.7) & m)[0]
    i1 = c[np.nanargmax(r[c])]
    g = np.gradient(r, x)
    s = np.arange(i0, i1+1)
    return refine(x, g, s[np.nanargmax(g[s])])


def axis(f):
    d = np.load(OUT/f)
    ts = sorted(float(k[1:].split("_")[0]) for k in d.files
                if k.endswith("_DensityMEAN"))
    M = np.array([d["t%.3f_DensityMEAN" % t] for t in ts])/RHO_J
    x = (np.arange(M.shape[1])+0.5)*(24.0/M.shape[1])
    return x, M[-1]


xL5, rL5 = axis("AX_L5_afdhllc_weno.npz")
xL3, rL3 = axis("AX_s3d_afdhllc_weno.npz")
b = np.load(PB/"STATS3D_baseline.npz")
bx, br = b["xx"], b["ax_DensityMEAN"][0]/RHO_J
bfin = np.isfinite(br)
e = np.load(PB/"panda_m142den_axis_EPAPS.npz")
ex, er = e["x"], e["rho"]

fL5, fL3 = xfc_of(xL5, rL5), xfc_of(xL3, rL3)
print("x_fc   L5 %.4f   L3 %.4f   baseline %.4f   experiment %.3f"
      % (fL5, fL3, xfc_of(bx, br), EXP["x_fc"]))

SL = {v: np.load(OUT/("SL_L5_weno_%s.npz" % v))
      for v in ("DensityMEAN", "x_velocityMEAN")}
SB = {v: np.load(OUT/("SL_BASE_%s.npz" % v))
      for v in ("DensityMEAN", "x_velocityMEAN")}
for tag, S in (("AFD-WENO L5", SL), ("baseline", SB)):
    for v, s in S.items():
        print("%-12s %-16s %s  dx %.1f um  t %.3f ms  %.4g .. %.4g"
              % (tag, v, s["field"].shape, s["dx_um"], s["time"]*1e3,
                 np.nanmin(s["field"]), np.nanmax(s["field"])))

plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True,
                     "legend.frameon": False})
fig = plt.figure(figsize=(11.4, 9.4))
gs = fig.add_gridspec(3, 1, height_ratios=[1.05, 1, 1], hspace=0.30,
                      left=0.075, right=0.905, top=0.965, bottom=0.065)
a1 = fig.add_subplot(gs[0])
a2 = fig.add_subplot(gs[1])
a3 = fig.add_subplot(gs[2], sharex=a2)

# ------------------------------------------------------------- centreline
a1.plot(bx[bfin], br[bfin], color="0.55", lw=1.4, ls="--",
        label="3D baseline WENO-Z5 (198.4 $\\mu$m, axis data to 2 $D_e$)")
a1.plot(ex, er, "o", color="k", ms=4.2, mfc="none", mew=0.9,
        label="Panda & Seasholtz 1999")
a1.plot(xL3[xL3 <= 7], rL3[xL3 <= 7], color="#D55E00", lw=1.3,
        label="L3, 396.9 $\\mu$m   $x_{fc}$ = %.4f" % fL3)
a1.plot(xL5[xL5 <= 7], rL5[xL5 <= 7], color="#0072B2", lw=1.7,
        label="L5, lip 99.2 $\\mu$m   $x_{fc}$ = %.4f" % fL5)
a1.axvline(EXP["x_fc"], color="0.45", lw=1.0)
a1.set_xlim(0, 7)
a1.set_ylim(0.30, 1.55)
a1.set_xlabel("$x/D_e$", labelpad=1)
a1.set_ylabel(r"$\overline{\rho}/\rho_j$  on the axis")
a1.grid(True, ls="--", lw=0.5, color="0.93")
a1.set_axisbelow(True)
a1.legend(fontsize=7.8, loc="upper right", handlelength=1.9, ncol=2,
          columnspacing=1.4)
a1.text(0.010, 0.05, "(a)  centreline", transform=a1.transAxes, fontsize=9.5)

# ------------------------------------------------------------- fields
for ax, v, lab, cmap, nrm, pl in (
        (a2, "DensityMEAN", r"$\overline{\rho}$  [kg m$^{-3}$]", "inferno",
         Normalize(0.4, 2.6), "(b)  mean density"),
        (a3, "x_velocityMEAN", r"$\overline{u}$  [m s$^{-1}$]", "viridis",
         Normalize(-60, 620), "(c)  mean axial velocity")):
    s, sb = SL[v], SB[v]
    F = s["field"].T                                   # (y, x)  this run
    G = sb["field"].T                                  #          baseline
    up = s["y"] >= 0
    dn = sb["y"] < 0
    im = ax.pcolormesh(s["x"], s["y"][up], F[up], cmap=cmap, norm=nrm,
                       shading="nearest", rasterized=True)
    ax.pcolormesh(sb["x"], sb["y"][dn], G[dn], cmap=cmap, norm=nrm,
                  shading="nearest", rasterized=True)
    ax.contour(s["x"], s["y"][up], F[up], levels=12, colors="w",
               linewidths=0.35, alpha=0.45)
    ax.contour(sb["x"], sb["y"][dn], G[dn], levels=12, colors="w",
               linewidths=0.35, alpha=0.45)
    ax.axhline(0.0, color="w", lw=1.1)
    ax.axhline(0.5, color="w", lw=0.7, ls=":", alpha=0.8)
    ax.axhline(-0.5, color="w", lw=0.7, ls=":", alpha=0.8)
    ax.axvline(fL5, color="#4dd0e1", lw=1.0)
    ax.axvline(0.8828, color="#ffb74d", lw=1.0)
    ax.text(0.010, 0.79, "AFD HLLC + WENO-Z5, L5,  $t$ = 17.6 ms",
            transform=ax.transAxes, fontsize=8.0, color="w")
    ax.text(0.010, 0.17, "3D baseline, WENO-Z5,  $t$ = 10.16 ms",
            transform=ax.transAxes, fontsize=8.0, color="w")
    ax.set_xlim(0, 7)
    ax.set_ylim(-1.6, 1.6)
    ax.set_aspect("equal")
    ax.set_ylabel("$y/D_e$")
    ax.set_xlabel("$x/D_e$", labelpad=1)
    ax.text(0.010, 0.90, pl, transform=ax.transAxes, fontsize=9.5, color="w")
    cb = fig.colorbar(im, ax=ax, pad=0.012, fraction=0.030, aspect=13)
    cb.set_label(lab, fontsize=8.6)
    cb.ax.tick_params(labelsize=7.6)

for ax in (a2, a3):
    ax.text(0.985, 0.055, "cyan $x_{fc}$ = %.4f (upper)     amber 0.8828 (lower)"
            "     dotted $r=\\pm R$" % fL5, transform=ax.transAxes,
            ha="right", fontsize=7.4, color="w")

fig.savefig(HERE/"l5_weno_field.png", dpi=220)
fig.savefig(HERE/"l5_weno_field.pdf")
print("\nwrote l5_weno_field.png/pdf")
