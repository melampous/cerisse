#!/usr/bin/env python3
"""Mean Mach field (M reconstructed from time-averaged primitives) of the
three numerical schemes over 0 < x/De < 16, z=0 plane. The M=1 contour is
outlined; the end of the supersonic core (last axial location where the
centreline mean Mach crosses below 1 and stays subsonic) is marked with a
white star and its x/De printed."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7})
HERE = Path(__file__).resolve().parent
XLIM, RLIM = 16.0, 1.4
CASES = [("baseline", "WENO-Z5 NS (baseline)"),
         ("teno", "TENO5 NS"),
         ("euler", "Euler (WENO-Z5)")]


def core_end(x, y, M):
    """Last x where the centreline (|y| minimal) mean Mach is >=1, i.e.
    the downstream end of the supersonic core."""
    jc = np.argsort(np.abs(y))[:2]
    Mc = M[:, jc].mean(1)
    sup = np.where(Mc >= 1.0)[0]
    if len(sup) == 0:
        return np.nan
    return float(x[sup[-1]])


fig, axes = plt.subplots(3, 1, figsize=(9.6, 4.1), sharex=True,
                         gridspec_kw={"hspace": 0.13})
im = None
for ax, (v, lab) in zip(axes, CASES):
    d = np.load(HERE / ("MACH16_%s.npz" % v))
    x, y, M = d["x"], d["y"], d["M"]
    im = ax.imshow(M.T, origin="lower", extent=[x[0], x[-1], y[0], y[-1]],
                   vmin=0, vmax=2.0, cmap="turbo", aspect="equal",
                   interpolation="nearest", rasterized=True)
    # M = 1 contour
    ax.contour(x, y, M.T, levels=[1.0], colors="w", linewidths=0.7)
    xc = core_end(x, y, M)
    if np.isfinite(xc):
        ax.plot(xc, 0, "*", color="w", ms=11, mec="k", mew=0.6, zorder=6)
        ax.axvline(xc, color="w", lw=0.7, ls="--", alpha=0.7)
        ax.text(xc + 0.15, 0.95, "core end\n$x/D_e$=%.1f" % xc,
                fontsize=7.5, color="w", va="top",
                bbox=dict(facecolor="0.15", edgecolor="none", alpha=0.5,
                          boxstyle="round,pad=0.2"))
    ax.set_xlim(0, XLIM); ax.set_ylim(-RLIM, RLIM)
    ax.set_ylabel("$r/D_e$")
    ax.text(0.011, 0.90, lab, transform=ax.transAxes, fontsize=8.5,
            va="top", color="w",
            bbox=dict(facecolor="0.15", edgecolor="none", alpha=0.55,
                      boxstyle="round,pad=0.25"))
axes[-1].set_xlabel("$x/D_e$")
fig.subplots_adjust(left=0.065, right=0.90, top=0.985, bottom=0.075)
cax = fig.add_axes([0.915, 0.20, 0.014, 0.6])
cb = fig.colorbar(im, cax=cax)
cb.set_label(r"$M_{\bar q}$ (from time-mean primitives)", fontsize=9)
cb.ax.tick_params(labelsize=8)
fig.savefig(HERE / "mach16_schemes.png", dpi=250)
fig.savefig(HERE / "mach16_schemes.pdf", dpi=250)
print("wrote mach16_schemes.png/pdf")
for v, lab in CASES:
    d = np.load(HERE / ("MACH16_%s.npz" % v))
    print("%-10s supersonic core end x/De = %.2f"
          % (lab, core_end(d["x"], d["y"], d["M"])))
