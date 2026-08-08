#!/usr/bin/env python3
"""Prescribed nozzle-exit velocity profiles of the eight completed cases, and
the corresponding boundary-layer integral parameters.

Panel (a) shows each profile on the actual computational mesh at the lip
(dx = 49.6 um, level 5), so the number of cells resolving the shear layer is
directly visible. Panel (b) plots the measured first shock-cell closure
position against the displacement thickness, with sampling uncertainty bars
obtained from seven independent 0.25 ms averaging windows per case.
"""
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
from design_2d_sweep import scalars, fvel, R, D, DX_LIP

HERE = Path(__file__).resolve().parent
DX = 49.6e-6                      # realised lip resolution, level 5

# name, a[um], s, form, measured x_fc, sampling sigma
MEAS = [
    ("A0_tophat",     0.0, 3.0, "shifted", 0.9000, 0.0032, "#404040"),
    ("A1_a100",     100.0, 3.0, "shifted", 0.9100, 0.0025, "#7b3294"),
    ("A2_a170",     170.0, 3.0, "shifted", 0.8975, 0.0045, "#b3251e"),
    ("A3_a254",     254.0, 3.0, "shifted", 0.8840, 0.0018, "#e08214"),
    ("S4_a170_s4",  170.0, 4.0, "shifted", 0.8968, 0.0023, "#c2a5cf"),
    ("S5_a170_s5",  170.0, 5.0, "shifted", 0.8802, 0.0038, "#8073ac"),
    ("B1_wall680",  679.5, 3.0, "wall",    0.9015, 0.0020, "#1f4e9c"),
    ("C1_wall200",  200.0, 3.0, "wall",    0.9061, 0.0021, "#2ca02c"),
]

plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.6, 3.5),
                             gridspec_kw={"width_ratios": [1.0, 1.05]})

# --------------------------------------------------------------- panel (a)
RLO, RHI = 0.395, 0.508
for k, v in enumerate(np.arange(0.0, R + DX, DX) / D):
    if RLO <= v <= RHI:
        a1.axhline(v, color="0.86", lw=0.35, zorder=0)
        if k % 2 == 0:
            a1.axhspan(v, min(v + DX / D, RHI), color="0.965", zorder=-1)
r = np.linspace(0.34 * D, R, 6000)
rows = []
for nm, a, s, form, xfc, sig, col in MEAS:
    f = fvel(r, a * 1e-6, form, s)
    ls = "-" if form == "shifted" else "--"
    a1.plot(np.r_[f, 0.0], np.r_[r / D, R / D], color=col, lw=1.3, ls=ls,
            zorder=4)
    sc = scalars(a * 1e-6, form, s)
    rows.append((nm, a, s, form, sc, xfc, sig, col))
a1.axhline(0.5, color="k", lw=1.1, zorder=5)
a1.text(0.5, 0.975, "nozzle lip  $r = R$", transform=a1.transAxes,
        ha="center", fontsize=7.5)
a1.set_xlim(-0.04, 1.08); a1.set_ylim(RLO, RHI)
a1.set_xlabel("$u_z / U_e$")
a1.set_ylabel("$r / D_e$")
a1.text(0.03, 0.05,
        "mesh lines = computational cells\n$\\Delta x$ = 49.6 $\\mu$m (level 5)",
        transform=a1.transAxes, fontsize=6.6, color="0.35")
a1.text(0.965, 0.05, "(a)", transform=a1.transAxes, ha="right", fontsize=9)
a1.legend(handles=[Line2D([], [], color=c, lw=1.3,
                          ls="-" if fm == "shifted" else "--",
                          label="%s  $\\delta_\\omega$=%.0f $\\mu$m (%.1f cells)"
                                % (n.split("_")[0], sc["dw"] * 1e6,
                                   sc["dw"] / DX))
                   for n, a, s, fm, sc, xf, sg, c in rows],
           fontsize=6.0, frameon=False, loc="lower left", labelspacing=0.22,
           handlelength=1.4, bbox_to_anchor=(0.0, 0.16))

# --------------------------------------------------------------- panel (b)
for nm, a, s, form, sc, xfc, sig, col in rows:
    mk = "o" if form == "shifted" else "s"
    a2.errorbar(sc["dstar_c"] * 1e6, xfc, yerr=sig, fmt=mk, color=col,
                ms=6, mec="w", mew=0.7, capsize=3, lw=1.1, zorder=4)
    a2.annotate(nm.split("_")[0], (sc["dstar_c"] * 1e6, xfc),
                textcoords="offset points", xytext=(6, 5), fontsize=6.3,
                color=col)
xs = np.array([r0[4]["dstar_c"] * 1e6 for r0 in rows])
ys = np.array([r0[5] for r0 in rows])
m, c = np.polyfit(xs, ys, 1)
xx = np.linspace(0, 1300, 40)
a2.plot(xx, m * xx + c, color="0.5", lw=1.0, zorder=2,
        label="least squares, all 8")
a2.plot(xx, 0.9443 - 1.9461 * xx / (D * 1e6), color="#b3251e", lw=1.0,
        ls="--", zorder=2, label="3D four-case correlation")
a2.set_xlabel(r"displacement thickness  $\delta^*_c$  [$\mu$m]")
a2.set_ylabel("$x_{fc} / D_e$   (time averaged, 2.0 ms)")
a2.set_xlim(-60, 1320)
a2.grid(True, ls="--", lw=0.5, color="0.9"); a2.set_axisbelow(True)
a2.legend(fontsize=6.6, frameon=False, loc="lower left", handlelength=1.8)
a2.text(0.965, 0.05, "(b)", transform=a2.transAxes, ha="right", fontsize=9)
fig.tight_layout(w_pad=1.3)
fig.savefig(HERE / "profiles_and_xfc.png", dpi=260)
fig.savefig(HERE / "profiles_and_xfc.pdf")
print("wrote profiles_and_xfc.png/pdf\n")

# ------------------------------------------------------------------ table
hdr = ("%-12s %-8s %6s %5s | %7s %7s %8s %7s | %6s | %7s %7s"
       % ("case", "form", "a[um]", "s", "dw[um]", "th[um]", "d*[um]", "C_d",
          "cells", "x_fc", "sigma"))
print(hdr); print("-" * len(hdr))
for nm, a, s, form, sc, xfc, sig, col in sorted(rows, key=lambda q: q[4]["dstar_c"]):
    print("%-12s %-8s %6.1f %5.1f | %7.0f %7.1f %8.1f %7.4f | %6.1f | %7.4f %7.4f"
          % (nm, form, a, s, sc["dw"] * 1e6, sc["theta_c"] * 1e6,
             sc["dstar_c"] * 1e6, sc["Cd"], sc["dw"] / DX, xfc, sig))
print("\ncells = delta_omega / 49.6 um, i.e. cells across the vorticity thickness")
print("sigma = uncertainty of the 2.0 ms mean, from 7 independent 0.25 ms windows")
