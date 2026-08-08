#!/usr/bin/env python3
"""Numerical schlieren of the top-hat jet at NPR 15 and 20, with the Mach-disk
region enlarged beside each full view.

Mean density on the axisymmetric Level 4 grid, composited finest-wins at
99.2 um and mirrored about the axis.

    S = exp( -k |grad rho| / g_ref )

The reference is NOT the global maximum. The nozzle lip is a discontinuity in
the imposed profile and its gradient reaches 4.6e4, seventy times the 99.9th
percentile of everything downstream, so normalising by the maximum leaves the
entire shock structure inside the first two per cent of the colour range and
the panel comes out blank. g_ref is therefore the 99.5th percentile taken over
x > 0.8 D_e, which is the part of the field the figure is about; the lip
saturates to black and is clipped, which is the honest outcome for a
zero-thickness inlet.

The enlargement sits beside the full view rather than on top of it, so nothing
in the plume is covered.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, ConnectionPatch

HERE = Path(__file__).resolve().parent
F = HERE/"fields2d"
D = 0.0254
K = 5.0
ZW, ZH = 0.62, 0.78            # half-width and half-height of the zoom box

CASES = [("NPR$_0$ = 15", "FLD_N15_DensityMEAN.npz", 2.5156),
         ("NPR$_0$ = 20", "FLD_N20_DensityMEAN.npz", 2.9688)]


def schlieren(f):
    z = np.load(F/f)
    rho = z["F"].astype(np.float64)          # (r, z)
    dr = float(z["dr"])
    r, x = z["r"]/D, z["z"]/D
    gr, gz = np.gradient(rho, dr, dr)
    g = np.sqrt(gr*gr + gz*gz)
    gref = np.nanpercentile(g[:, x > 0.8], 99.5)
    s = np.exp(-K*np.clip(g/gref, 0.0, 1.0))
    rr = np.concatenate([-r[::-1], r])
    ss = np.vstack([s[::-1, :], s])
    return x, rr, ss, float(z["time"])*1e3, gref, np.nanmax(g)


plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})
fig, AX = plt.subplots(2, 2, figsize=(11.6, 5.9),
                       gridspec_kw=dict(width_ratios=[2.45, 1.0]))
fig.subplots_adjust(left=0.055, right=0.988, top=0.955, bottom=0.088,
                    wspace=0.105, hspace=0.225)

print("%-14s %8s %12s %12s" % ("case", "t[ms]", "g_ref(p99.5)", "g_max(lip)"))
for row, (lab, f, xfc) in enumerate(CASES):
    x, rr, ss, t, gref, gmax = schlieren(f)
    print("%-14s %8.3f %12.1f %12.1f"
          % (lab.replace("$", ""), t, gref, gmax))
    a, b = AX[row, 0], AX[row, 1]
    for ax in (a, b):
        ax.pcolormesh(x, rr, ss, cmap="gray", shading="nearest",
                      rasterized=True, vmin=0.0, vmax=1.0)
        ax.axvline(xfc, color="#e8482c", lw=1.0, ls="--", zorder=5)
        ax.set_aspect("equal")
    a.set_xlim(0, 6.0)
    a.set_ylim(-1.6, 1.6)
    a.set_ylabel("$r/D_e$")
    a.text(0.010, 0.90, "%s,   $t$ = %.2f ms" % (lab, t), transform=a.transAxes,
           fontsize=10)
    a.add_patch(Rectangle((xfc-ZW, -ZH), 2*ZW, 2*ZH, fill=False,
                          edgecolor="#e8482c", lw=1.1, zorder=6))
    b.set_xlim(xfc-ZW, xfc+ZW)
    b.set_ylim(-ZH, ZH)
    b.set_ylabel("$r/D_e$", labelpad=1)
    b.text(0.035, 0.90, r"$x_{fc}$ = %.3f" % xfc, transform=b.transAxes,
           fontsize=9.5, color="#e8482c")
    for sp in b.spines.values():
        sp.set_color("#e8482c")
        sp.set_linewidth(1.1)
    for yy in (-ZH, ZH):
        fig.add_artist(ConnectionPatch(
            xyA=(xfc+ZW, yy), coordsA=a.transData,
            xyB=(xfc-ZW, yy), coordsB=b.transData,
            color="#e8482c", lw=0.7, ls=":", zorder=1))

AX[1, 0].set_xlabel("$x/D_e$")
AX[1, 1].set_xlabel("$x/D_e$")
fig.savefig(HERE/"npr_schlieren_machdisk.png", dpi=300)
fig.savefig(HERE/"npr_schlieren_machdisk.pdf")
print("\nwrote npr_schlieren_machdisk.png/pdf")
