#!/usr/bin/env python3
"""The reconstruction x flux-path matrix, one panel per combination, each
against the Panda & Seasholtz (1999) centreline density.

Rows are the flux path: the characteristic-wise flux-vector splitting with a
per-face local Lax-Friedrichs flux (the weno_t template, JET_EULER_SCHEME_ID
0 and 2), and the alternative flux discretisation with an HLLC Riemann solver
(afd_hllc_*_t, IDs 5 and 1). Columns are the point-value reconstruction,
WENO-Z5 and TENO5. The inlet profile is the same shifted tanh, a = 254 um, in
all four.

The two rows do NOT share a statistics window or an axis extent:

  LLF row    earlier production runs, windows [6.5, 10.16] and [6.5, 9.30] ms,
             centreline composited from the three innermost near-axis rings
  HLLC row   the L5 legs of 2026-07-26/27, window [15.0, 17.6] ms, centreline
             over the full 24 D_e domain

Both rows resolve the lip at 99.2 um. The windows differ, so the comparison
that carries weight is the shock-cell structure, not the absolute level of a
late-time transient.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
PB = HERE.parent
OUT = HERE/"sweep_out"
RHO_J = 1.6413
D = 0.0254
EXP_XFC = 0.888


def native(x, r):
    """One sample per real cell. Away from the finest patches the stored axis
    value is a coarser cell replicated across 2, 4 or 8 slots of the level-5
    abscissa (198.4, 396.9, 793.7 um against 99.2 um), which draws as a
    staircase. Collapse each run of equal values to its midpoint."""
    m = np.isfinite(r)
    x, r = x[m], r[m]
    edge = np.r_[0, np.flatnonzero(np.diff(r) != 0) + 1, r.size]
    xc = np.array([0.5*(x[a] + x[b-1]) for a, b in zip(edge[:-1], edge[1:])])
    return xc, r[edge[:-1]]


def refine(x, y, i):
    if i <= 0 or i >= len(y)-1:
        return x[i]
    d = y[i-1]-2*y[i]+y[i+1]
    return x[i] if d == 0 else x[i]+0.5*(y[i-1]-y[i+1])/d*(x[1]-x[0])


def xfc_of(x, r):
    """Location of the steepest density rise across the first shock cell."""
    m = np.isfinite(r)
    e = np.where((x >= 0.35) & (x <= 1.0) & m)[0]
    if e.size == 0:
        return np.nan
    i0 = e[np.nanargmin(r[e])]
    c = np.where((np.arange(x.size) > i0) & (x <= 1.7) & m)[0]
    if c.size == 0:
        return np.nan
    i1 = c[np.nanargmax(r[c])]
    g = np.gradient(r, x)   # x is non-uniform; gradient handles that
    s = np.arange(i0, i1+1)
    return refine(x, g, s[np.nanargmax(g[s])])


def from_stats3d(f):
    """ax_DensityMEAN is (near-axis ring, x), not (level, x). Ring 0 hugs the
    axis but the coarse levels carry no cell at r < 0.01 D, so it runs dry at
    2 D_e. The convention already used in plot_stats3d_vs_exp.py is to lay ring
    0 down and fill its gaps from rings 1 and 2, which reaches 10 D_e."""
    z = np.load(PB/f)
    x = z["xx"]
    rings = z["ax_DensityMEAN"]/RHO_J
    r = rings[0].copy()
    for k in (1, 2):
        gap = ~np.isfinite(r)
        r[gap] = rings[k][gap]
    xc, rc = native(x, r)
    return xc, rc, z["time"]*1e3


def from_axpack(f):
    """Multi-frame pack: keys t<ms>_DensityMEAN. Take the last frame."""
    z = np.load(OUT/f)
    ts = sorted(float(k[1:].split("_")[0]) for k in z.files
                if k.endswith("_DensityMEAN"))
    r = z["t%.3f_DensityMEAN" % ts[-1]]/RHO_J
    xc, rc = native((np.arange(r.size)+0.5)*(24.0/r.size), r)
    return xc, rc, ts[-1]


def from_axframe(f):
    z = np.load(OUT/f)
    xc, rc = native(z["x"]/D, z["ax_DensityMEAN"]/RHO_J)
    return xc, rc, float(z["time"])*1e3


CELLS = [
    dict(row=0, col=0, key="llfweno",  sid=0,
         name="FVS + per-face LLF,  WENO-Z5",
         win="[6.5, 10.16] ms", data=lambda: from_stats3d("STATS3D_baseline.npz")),
    dict(row=0, col=1, key="llfteno",  sid=2,
         name="FVS + per-face LLF,  TENO5",
         win="[6.5, 9.30] ms",  data=lambda: from_stats3d("STATS3D_teno.npz")),
    dict(row=1, col=0, key="hllcweno", sid=5,
         name="AFD + HLLC,  WENO-Z5",
         win="[15.0, 17.6] ms", data=lambda: from_axpack("AX_L5_afdhllc_weno.npz")),
    dict(row=1, col=1, key="hllcteno", sid=1,
         name="AFD + HLLC,  TENO5",
         win="[15.0, 17.6] ms", data=lambda: from_axframe("AX_L5_s3d_afdhllc_plt15761.npz")),
]

e = np.load(PB/"panda_m142den_axis_EPAPS.npz")
ex, er = e["x"], e["rho"]

print("%-10s %-34s %-16s %8s %8s %9s"
      % ("cell", "scheme", "window", "t[ms]", "x_fc", "x_max"))
for c in CELLS:
    c["x"], c["r"], c["t"] = c["data"]()
    c["xfc"] = xfc_of(c["x"], c["r"])
    print("%-10s %-34s %-16s %8.2f %8.4f %9.2f"
          % (c["key"], c["name"], c["win"], c["t"], c["xfc"], c["x"].max()))
print("%-10s %-34s %-16s %8s %8.4f %9.2f"
      % ("experiment", "Panda & Seasholtz 1999", "", "-", EXP_XFC, ex.max()))

# ------------------------------------------------------------------ figure
plt.rcParams.update({"font.size": 10, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True,
                     "legend.frameon": False})
SIM = "#0072B2"
ROWLAB = ["FVS + per-face LLF", "AFD + HLLC"]
COLLAB = ["WENO-Z5", "TENO5"]

fig, AX = plt.subplots(2, 2, figsize=(11.0, 7.6), sharex=True, sharey=True)
fig.subplots_adjust(left=0.098, right=0.988, top=0.848, bottom=0.078,
                    wspace=0.070, hspace=0.105)

for c in CELLS:
    ax = AX[c["row"], c["col"]]
    short = c["x"].max() < 7.0
    if short:                                   # only where the data really stops
        ax.axvspan(c["x"].max(), 7.2, color="0.945", zorder=0)
        ax.text(0.5*(c["x"].max()+7.2), 1.42, "no centreline data",
                ha="center", va="center", fontsize=8.8, color="0.52")
    ax.axvline(EXP_XFC, color="0.62", lw=0.9, zorder=1)
    ax.plot(ex, er, "o", color="k", ms=4.0, mfc="none", mew=0.9, zorder=4,
            label="Panda & Seasholtz 1999")
    ax.plot(c["x"], c["r"], color=SIM, lw=1.6, zorder=5, label="simulation")
    ax.axvline(c["xfc"], color=SIM, lw=1.0, ls="--", zorder=3)

    ax.set_xlim(0, 7.2)
    ax.set_ylim(0.30, 1.82)                     # header strip above the data
    ax.axhline(1.50, color="0.85", lw=0.6, zorder=2)
    ax.grid(True, ls="--", lw=0.5, color="0.93")
    ax.set_axisbelow(True)
    ax.text(0.018, 0.958, "(%s)" % "abcd"[2*c["row"]+c["col"]],
            transform=ax.transAxes, fontsize=10.4, va="top")
    ax.text(0.085, 0.955, "mean over %s" % c["win"], transform=ax.transAxes,
            fontsize=8.8, va="top", color="0.30")
    ax.text(0.982, 0.955,
            "$x_{fc}$ = %.4f   (%+.1f %%)"
            % (c["xfc"], 100*(c["xfc"]-EXP_XFC)/EXP_XFC),
            transform=ax.transAxes, fontsize=8.8, va="top", ha="right",
            color="0.30")
    if c["col"] == 0:
        ax.set_ylabel(r"$\overline{\rho}/\rho_j$   on the axis")
    if c["row"] == 1:
        ax.set_xlabel("$x/D_e$")
    if c["row"] == 0:
        ax.set_title(COLLAB[c["col"]], fontsize=11.5, pad=8)

for i, lab in enumerate(ROWLAB):               # matrix row labels
    b = AX[i, 0].get_position()
    fig.text(0.016, 0.5*(b.y0+b.y1), lab, rotation=90, va="center",
             ha="center", fontsize=11.5)

h, l = AX[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.5, 0.965), ncol=2,
           fontsize=9.4, handlelength=1.9, columnspacing=2.6)
fig.text(0.5, 0.985, "Centreline mean density: reconstruction $\\times$ flux path",
         ha="center", fontsize=12)
fig.text(0.5, 0.902,
         "grey vertical line: measured $x_{fc}$ = %.3f          "
         "dashed: computed $x_{fc}$" % EXP_XFC,
         ha="center", fontsize=8.8, color="0.42")

fig.savefig(HERE/"scheme_2x2_axis.png", dpi=250)
fig.savefig(HERE/"scheme_2x2_axis.pdf")
print("\nwrote scheme_2x2_axis.png/pdf")
