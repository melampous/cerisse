#!/usr/bin/env python3
"""Centreline mean-density comparison, latest M_j=1.42 runs vs experiment.

Curves:
  - Panda & Seasholtz (1999) M_j=1.42 experiment  (filled black circles,
    already normalised by rho_j, x already in x/D)
  - baseline  : order=2 viscous  (m142_3d_gpu/l5_prod_cb, plt09170, t=10.16 ms)
  - order=4   : 4th-order central viscous  (plt13000, t=14.37 ms)
  - WALE      : order=4 + WALE SGS         (plt13000, t=14.37 ms)

Simulation curves use the accumulated DensityMEAN field extracted along the
jet axis (y=z=0) with AMReX fextract at the finest available level.
Normalisation follows shockcell/plot_paper_cmp_panda.py: D=0.0254 m,
rho_j=1.6413 kg/m^3.  A 0.05 D box filter lightly de-noises the profiles.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path

plt.rcParams.update({
    "font.size": 9, "axes.linewidth": 0.6,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
})

HERE = Path(__file__).parent
D = 0.0254
RHO_J = 1.6413


def boxfilt(a, x, wD=0.05):
    dx = np.median(np.diff(x))
    w = max(3, int(round(wD / dx)))
    if w % 2 == 0:
        w += 1
    out = np.convolve(a, np.ones(w) / w, mode="same")
    out[:w // 2] = a[:w // 2]
    out[-(w // 2):] = a[-(w // 2):]
    return out


def load_axis(fn):
    d = np.loadtxt(HERE / fn)
    x = d[:, 0] / D                 # -> x/D
    rho_mean = d[:, 1] / RHO_J      # DensityMEAN / rho_j
    order = np.argsort(x)
    return x[order], rho_mean[order]


# --- experiment (already normalised) ---
xd = np.load(HERE.parent / "panda_m142den_axis_EPAPS.npz")

# --- simulations ---
xb, rb = load_axis("axis_baseline_order2.txt")
xa, ra = load_axis("axis_A_order4.txt")
xw, rw = load_axis("axis_B_wale.txt")

C_BASE, C_O4, C_WALE = "#6b6b6b", "#b3251e", "#1f8a4c"

fig, ax = plt.subplots(figsize=(7.2, 3.1))
ax.set_axisbelow(True)
ax.grid(True, ls="--", lw=0.5, color="0.87")

ax.plot(xb, boxfilt(rb, xb), color=C_BASE, lw=1.3, ls=(0, (5, 2)),
        zorder=3, label="baseline (order 2), t=10.2 ms")
ax.plot(xa, boxfilt(ra, xa), color=C_O4, lw=1.6,
        zorder=4, label="order 4, t=14.4 ms")
ax.plot(xw, boxfilt(rw, xw), color=C_WALE, lw=1.6,
        zorder=4, label="order 4 + WALE, t=14.4 ms")
ax.plot(xd["x"], xd["rho"], "o", color="k", ms=4.0, mfc="k",
        mec="w", mew=0.4, zorder=6, label="$M_j=1.42$ experiment (Panda 1999)")

ax.set_xlim(0, 8)
ax.set_ylim(0.30, 1.55)
ax.set_xlabel("$x/D_e$", labelpad=2)
ax.set_ylabel(r"$\langle\rho\rangle/\rho_j$", labelpad=2)
ax.tick_params(length=2.5)
ax.legend(fontsize=7.2, frameon=False, loc="upper right",
          handlelength=2.2, labelspacing=0.35)
ax.set_title("Centreline mean density — latest $M_j=1.42$ runs vs experiment",
             fontsize=8.5, pad=4)

fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(HERE / f"axis_rho_cmp_20260806.{ext}", dpi=300)
print("wrote axis_rho_cmp_20260806.png / .pdf")

# quick quantitative note: first two shock-cell minima/maxima vs experiment
def first_extrema(x, r, xmax=6.0):
    m = (x > 0.4) & (x < xmax)
    xx, rr = x[m], boxfilt(r, x)[m]
    return xx, rr

print("\n# sim near-field <rho>/rho_j at x/D = 1,2,3,4,5:")
for lbl, (x, r) in {"base": (xb, rb), "o4": (xa, ra), "wale": (xw, rw)}.items():
    rf = boxfilt(r, x)
    vals = [f"{np.interp(s, x, rf):.3f}" for s in (1, 2, 3, 4, 5)]
    print(f"  {lbl:5s}: " + "  ".join(vals))
xe, re = xd["x"], xd["rho"]
vals = [f"{np.interp(s, xe, re):.3f}" for s in (1, 2, 3, 4, 5)]
print(f"  exp  : " + "  ".join(vals))
