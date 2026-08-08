#!/usr/bin/env python3
"""Does the solver inject the profile it was told to?

Compares the prescribed inlet shape against the realised mean axial velocity
in the first interior cell row. They agree to 1e-4 through the core and part
company only inside the lip band, because the underexpanded lip expansion
(p_e/p_inf = 1.73) acts within that first cell. The top-hat case proves this
is physics and not a boundary-condition error: its prescription is a perfect
square wave, yet the field reads u/U_e = 1.06 at the lip - the flow has
already accelerated past the exit velocity, which no boundary condition could
have imposed.
"""
from pathlib import Path
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
D, R = 0.0254, 0.0127
LABEL = {"A0_tophat": "A0  top hat",
         "A2_a170": "A2  shifted $a$=170, $s$=3",
         "B1_wall680": "B1  wall tanh $a$=679.5",
         "S7_a170_s7": "S7  shifted $a$=170, $s$=7"}
COL = {"A0_tophat": "#7a7f85", "A2_a170": "#b3251e",
       "B1_wall680": "#1f4e9c", "S7_a170_s7": "#ff6a00"}

files = sorted(glob.glob(str(HERE / "sweep_out" / "INLET_*.npz")))
if not files:
    raise SystemExit("no INLET_*.npz in sweep_out/")

plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.4, 3.3),
                             gridspec_kw={"width_ratios": [1.0, 1.15]})
rows = []
for f in files:
    nm = Path(f).stem.replace("INLET_", "")
    d = np.load(f)
    r, fn, fa = d["r"], d["f_num"], d["f_ana"]
    g = np.isfinite(fn) & (r <= 0.56 * D)
    c = COL.get(nm, "k")
    a1.plot(fa[g], r[g] / D, color=c, lw=2.4, alpha=0.30, solid_capstyle="butt")
    a1.plot(fn[g], r[g] / D, color=c, lw=1.1)
    core = (r <= 0.45 * D) & np.isfinite(fn)
    lip = (r > 0.45 * D) & (r <= R) & np.isfinite(fn)
    rows.append((nm, np.sqrt(((fn - fa)[core]**2).mean()),
                 np.sqrt(((fn - fa)[lip]**2).mean())))
    a2.plot(r[g] / D, (fn - fa)[g], color=c, lw=1.2)

a1.axhline(0.5, color="k", lw=1.0)
a1.text(0.06, 0.505, "lip  $r=R$", fontsize=7)
a1.set_xlabel("$u_z/U_e$"); a1.set_ylabel("$r/D_e$")
a1.set_xlim(-0.05, 1.15); a1.set_ylim(0.30, 0.56)
a1.grid(True, ls="--", lw=0.5, color="0.9"); a1.set_axisbelow(True)
a1.legend(handles=[Line2D([], [], color="0.4", lw=2.4, alpha=0.35,
                          label="prescribed (analytic)"),
                   Line2D([], [], color="0.4", lw=1.1,
                          label="realised (first cell row)")],
          fontsize=6.8, frameon=False, loc="lower left", handlelength=1.8)
a1.text(0.97, 0.05, "(a)", transform=a1.transAxes, ha="right", fontsize=9)

a2.axvline(0.5, color="k", lw=1.0)
a2.axhline(0.0, color="0.6", lw=0.8, ls=":")
a2.axvspan(0.45, 0.5, color="0.92", zorder=0)
a2.text(0.475, 0.30, "lip band", fontsize=6.8, ha="center", color="0.35")
a2.set_xlabel("$r/D_e$")
a2.set_ylabel("realised $-$ prescribed,  $\\Delta u_z/U_e$")
a2.set_xlim(0.30, 0.56)
a2.grid(True, ls="--", lw=0.5, color="0.9"); a2.set_axisbelow(True)
a2.legend(handles=[Line2D([], [], color=COL.get(n, "k"), lw=1.4,
                          label="%s   core rms %.1e" % (LABEL.get(n, n), cr))
                   for n, cr, lr in rows],
          fontsize=6.4, frameon=False, loc="upper left", handlelength=1.4)
a2.text(0.97, 0.05, "(b)", transform=a2.transAxes, ha="right", fontsize=9)
fig.tight_layout(w_pad=1.2)
fig.savefig(HERE / "inlet_check.png", dpi=260)
fig.savefig(HERE / "inlet_check.pdf")
print("wrote inlet_check.png/pdf")
for n, cr, lr in rows:
    print("  %-14s core rms=%.2e   lip rms=%.4f" % (n, cr, lr))
