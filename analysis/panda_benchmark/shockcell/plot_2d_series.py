#!/usr/bin/env python3
"""2D RZ m142 (NPR=3.27) series: grid convergence (L4/L5/L6) + scheme comparison
(WENO-Z5 / TENO5 / HLLC) on the jet centerline, vs Panda & Seasholtz (1999) exp.
Data: AX2_*.npz produced by axis2d_lite.py (pseudo-time-mean of last N frames).
Output: fig19_2d_gridconv.png, fig20_2d_scheme.png"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).parent
D = 0.0254                 # nozzle exit diameter [m]
RHO_AMB = 1.1698           # ambient density [kg/m^3] (project convention)
P_AMB = 101325.0

exp = np.genfromtxt(HERE/"exp_m142_centerline.csv", delimiter=",", names=True)

def load(name):
    f = HERE/f"AX2_{name}.npz"
    if not f.exists():
        return None
    d = np.load(f, allow_pickle=True)
    return d

CASES_GRID = [  # (key, label, color)  gradient blues by refinement
    ("L3dev",   "L3  ($\\Delta x$=198 µm)", "#c6dbef"),
    ("L4",       "L4  ($\\Delta x$=99 µm)",  "#9ecae1"),
    ("L5full",  "L5  ($\\Delta x$=50 µm)",  "#4292c6"),
    ("L6weno",   "L6  ($\\Delta x$=25 µm)",  "#084594"),
]
CASES_SCHEME = [  # scheme comparison; distinct hues
    ("L6weno", "WENO-Z5", "#084594"),
    ("teno",   "TENO5",   "#d94801"),
    ("hllc",   "HLLC-MUSCL", "#238b45"),
]

def panel(ax, cases, qty, ylab, xmax=12.0):
    for key, lab, col in cases:
        d = load(key)
        if d is None:
            continue
        z = d["z"]/D
        if qty == "rho":
            y = d["pm_rho"]/RHO_AMB
        else:
            y = d["pm_p"]/P_AMB
        m = (z <= xmax) & np.isfinite(y)
        fin = int(d["finest"]); nf = int(d["nframes"]); tl = float(d["time_last"])*1e3
        ax.plot(z[m], y[m], color=col, lw=1.4,
                label=f"{lab}  [{nf}f, t≤{tl:.1f} ms]", path_effects=None)
    if qty == "rho":
        ax.plot(exp["x_over_D"], exp["rho_over_rhoamb"], "ko", ms=3.5, mfc="none",
                label="Panda & Seasholtz (1999)")
    ax.set_xlim(0, xmax); ax.set_xlabel("$x/D_e$")
    ax.set_ylabel(ylab); ax.grid(alpha=0.25, lw=0.4)
    ax.legend(fontsize=7.5, loc="upper right")

# ---- fig19: grid convergence ----
fig, axs = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
panel(axs[0], CASES_GRID, "rho", "$\\rho/\\rho_\\infty$")
panel(axs[1], CASES_GRID, "p",   "$p/p_\\infty$")
axs[0].set_title("2D RZ m142 (NPR$_0$=3.27) — centerline grid convergence, WENO-Z5, pseudo-time-mean")
fig.tight_layout()
fig.savefig(HERE/"fig19_2d_gridconv.png", dpi=170)
print("wrote fig19_2d_gridconv.png")

# ---- fig20: scheme comparison ----
fig, axs = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
panel(axs[0], CASES_SCHEME, "rho", "$\\rho/\\rho_\\infty$")
panel(axs[1], CASES_SCHEME, "p",   "$p/p_\\infty$")
axs[0].set_title("2D RZ m142 — numerical scheme comparison on centerline (each at its own finest level)")
fig.tight_layout()
fig.savefig(HERE/"fig20_2d_scheme.png", dpi=170)
print("wrote fig20_2d_scheme.png")

# ---- quick metrics: first pressure-cell lengths from p minima ----
from scipy.signal import argrelextrema
print("\nshock-cell landmarks (pressure minima positions, x/D):")
for key, lab, _ in CASES_GRID + CASES_SCHEME:
    d = load(key)
    if d is None:
        continue
    z = d["z"]/D; p = d["pm_p"]
    m = (z > 0.1) & (z < 12) & np.isfinite(p)
    zs, ps = z[m], p[m]
    # smooth lightly to suppress cell-level noise
    k = max(3, int(len(ps)*0.002)|1)
    psm = np.convolve(ps, np.ones(k)/k, mode="same")
    idx = argrelextrema(psm, np.less, order=int(0.25/np.mean(np.diff(zs))))[0]
    print(f"  {key:10s}: " + " ".join(f"{zs[i]:.2f}" for i in idx[:6]))
