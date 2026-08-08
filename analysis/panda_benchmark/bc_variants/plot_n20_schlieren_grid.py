#!/usr/bin/env python3
"""NPR 20 top hat, Level 3 against Level 4: numerical schlieren of the mean
density with the Mach-disk region enlarged.

Same binary, same prob.h, same tree; amr.max_level is the only difference, so
the pair shows what the grid does to the disk and to the reflected structure
leaving its rim.

    S = exp( -k |grad rho| / g_ref ),   g_ref = 99.5th percentile over x > 0.8

The reference is a percentile and not the maximum: the top-hat lip is a
discontinuity in the imposed profile and its gradient is around seventy times
anything downstream, so normalising by the maximum leaves the whole shock
system in the first two per cent of the range. Each level gets its own g_ref,
because the coarse grid cannot produce the same peak gradient as the fine one
and a shared reference would make Level 3 look uniformly weaker for a reason
that is only its cell size.

Both fields are AMR composites, so the gradient operator crosses refinement
boundaries. That shows as a faint rectangular texture away from the jet; it is
an artefact of differencing across a resolution jump, not flow structure. The
enlargement sits inside the finest patch, where no such boundary is crossed.
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
ZW, ZH = 0.62, 0.78

CASES = [("Level 3,  198 $\\mu$m", "FLD_N20L3_DensityMEAN.npz", 2.8125),
         ("Level 4,  99 $\\mu$m",  "FLD_N20_DensityMEAN.npz",   2.9688)]


def schlieren(f):
    z = np.load(F/f)
    rho = z["F"].astype(np.float64)
    dr = float(z["dr"])
    r, x = z["r"]/D, z["z"]/D
    gr, gz = np.gradient(rho, dr, dr)
    g = np.sqrt(gr*gr + gz*gz)
    gref = np.nanpercentile(g[:, x > 0.8], 99.5)
    s = np.exp(-K*np.clip(g/gref, 0.0, 1.0))
    rr = np.concatenate([-r[::-1], r])
    ss = np.vstack([s[::-1, :], s])
    return x, rr, ss, float(z["time"])*1e3, gref, dr*1e6


plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})
fig, AX = plt.subplots(2, 2, figsize=(11.6, 5.9),
                       gridspec_kw=dict(width_ratios=[2.45, 1.0]))
fig.subplots_adjust(left=0.055, right=0.988, top=0.955, bottom=0.088,
                    wspace=0.105, hspace=0.225)

print("%-22s %8s %8s %12s" % ("case", "dx[um]", "t[ms]", "g_ref(p99.5)"))
for row, (lab, f, xfc) in enumerate(CASES):
    x, rr, ss, t, gref, dx = schlieren(f)
    print("%-22s %8.1f %8.3f %12.1f"
          % (lab.replace("$", "").replace("\\mu", "u"), dx, t, gref))
    a, b = AX[row, 0], AX[row, 1]
    for ax in (a, b):
        ax.pcolormesh(x, rr, ss, cmap="gray", shading="nearest",
                      rasterized=True, vmin=0.0, vmax=1.0)
        ax.axvline(xfc, color="#e8482c", lw=1.0, ls="--", zorder=5)
        ax.set_aspect("equal")
    a.set_xlim(0, 6.0)
    a.set_ylim(-1.6, 1.6)
    a.set_ylabel("$r/D_e$")
    a.text(0.010, 0.90, "%s,   $t$ = %.2f ms" % (lab, t),
           transform=a.transAxes, fontsize=10)
    a.add_patch(Rectangle((xfc-ZW, -ZH), 2*ZW, 2*ZH, fill=False,
                          edgecolor="#e8482c", lw=1.1, zorder=6))
    b.set_xlim(xfc-ZW, xfc+ZW)
    b.set_ylim(-ZH, ZH)
    b.set_ylabel("$r/D_e$", labelpad=1)
    b.text(0.035, 0.90, r"$x_{fc}$ = %.4f" % xfc, transform=b.transAxes,
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
fig.savefig(HERE/"n20_schlieren_L3_L4.png", dpi=300)
fig.savefig(HERE/"n20_schlieren_L3_L4.pdf")
print("\nwrote n20_schlieren_L3_L4.png/pdf")
