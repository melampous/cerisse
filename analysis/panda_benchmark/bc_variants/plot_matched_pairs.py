#!/usr/bin/env python3
"""The two matched pairs that separate the inlet integral scales.

Each row holds one pair. Within a row the two profiles agree exactly in one
scalar and differ in the others, so the row isolates that scalar.

  top     matched vorticity thickness, delta_w = 200 um
          wall tanh a = 200 um   against   shifted tanh a = 100 um
          same peak shear, mass deficit differs by a factor of two

  bottom  matched displacement thickness, delta*_c = 533.0 um
          shifted tanh a = 170 um   against   wall tanh a = 679.5 um
          same mass deficit, peak shear differs by a factor of two

The right column is the shear itself, which is where the asymmetry between the
two families shows: the wall-attached form is a half tanh clipped at the lip,
so its steepest point sits exactly on r = R where f = 0; the shifted form is a
complete tanh placed inside the orifice, so its steepest point sits at its own
centre where f = 0.5. Matching delta_w matches the height of those two peaks
but not where they stand, which is why the momentum thickness of the pair does
not match even when delta_w does.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from design_2d_sweep import scalars, fvel, R, D

HERE = Path(__file__).resolve().parent
C_WALL, C_SHIFT = "#D55E00", "#0072B2"


def solve_wall(target, key):
    lo, hi = 1e-6, 4e-3
    for _ in range(90):
        mid = 0.5*(lo + hi)
        if scalars(mid, "wall")[key] < target:
            lo = mid
        else:
            hi = mid
    return 0.5*(lo + hi)


A_DSTAR = solve_wall(scalars(170e-6, "shifted", 3.0)["dstar_c"], "dstar_c")

PAIRS = [
    dict(title=r"matched vorticity thickness,   $\delta_\omega$ = 200 $\mu$m",
         rlo=11.85, rhi=12.95, glo=0.0, ghi=5.6,
         members=[("wall tanh,  $a$ = 200 $\\mu$m", 200e-6, "wall", 3.0, C_WALL, "-"),
                  ("shifted tanh,  $a$ = 100 $\\mu$m", 100e-6, "shifted", 3.0, C_SHIFT, "-")]),
    dict(title=r"matched displacement thickness,   $\delta^*_c$ = 533.0 $\mu$m",
         rlo=9.9, rhi=12.95, glo=0.0, ghi=3.2,
         members=[("wall tanh,  $a$ = 679.5 $\\mu$m", A_DSTAR, "wall", 3.0, C_WALL, "-"),
                  ("shifted tanh,  $a$ = 170 $\\mu$m", 170e-6, "shifted", 3.0, C_SHIFT, "-")]),
]

r = np.linspace(0.0, R, 800001)
rmm, Rmm = r*1e3, R*1e3

print("%-30s %9s %9s %9s %9s %7s %8s"
      % ("profile", "dw[um]", "th_c[um]", "d*_c[um]", "d*_i[um]", "H", "C_d"))
for pr in PAIRS:
    print("--", pr["title"])
    for lab, a, form, s, col, ls in pr["members"]:
        sc = scalars(a, form, s)
        pr.setdefault("sc", []).append(sc)
        print("   %-27s %9.1f %9.1f %9.2f %9.1f %7.2f %8.4f"
              % (lab.replace("$", "").replace("\\mu", "u").replace("\\", ""),
                 sc["dw"]*1e6, sc["theta_c"]*1e6, sc["dstar_c"]*1e6,
                 sc["dstar_i"]*1e6, sc["H"], sc["Cd"]))

plt.rcParams.update({"font.size": 10, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True,
                     "legend.frameon": False})
fig, AX = plt.subplots(2, 2, figsize=(11.6, 7.6))
fig.subplots_adjust(left=0.070, right=0.988, top=0.930, bottom=0.075,
                    wspace=0.185, hspace=0.325)

for row, pr in enumerate(PAIRS):
    a1, a2 = AX[row, 0], AX[row, 1]
    for (lab, a, form, s, col, ls), sc in zip(pr["members"], pr["sc"]):
        f = fvel(r, a, form, s)
        g = np.abs(np.gradient(f, rmm))
        a1.plot(f, rmm, color=col, ls=ls, lw=1.9, zorder=5, label=lab)
        a2.plot(rmm, g, color=col, ls=ls, lw=1.9, zorder=5, label=lab)
        i = int(np.argmax(g))
        a2.plot([rmm[i]], [g[i]], marker="v", color=col, ms=8, mec="w",
                mew=0.9, zorder=7)
        edge = rmm[i] > Rmm - 0.02
        a2.annotate("$f$ = %.2f" % f[i], (rmm[i], g[i]),
                    textcoords="offset points",
                    xytext=(-8 if edge else 0, 11),
                    ha="right" if edge else "center", fontsize=9, color=col)

    for ax in (a1, a2):
        ax.grid(True, ls="--", lw=0.5, color="0.93")
        ax.set_axisbelow(True)
    a1.axhspan(Rmm, 20.0, color="0.915", zorder=0)
    a1.axhline(Rmm, color="0.4", lw=0.9, zorder=2)
    a1.set_xlim(-0.02, 1.06)
    a1.set_ylim(pr["rlo"], pr["rhi"])
    a1.set_xlabel("profile function  $f = u/U_e$")
    a1.set_ylabel("radius  $r$  [mm]")
    a1.legend(loc="lower left", fontsize=9.6, handlelength=2.0)

    a2.axvspan(Rmm, 20.0, color="0.915", zorder=0)
    a2.axvline(Rmm, color="0.4", lw=0.9, zorder=2)
    a2.set_xlim(pr["rlo"], pr["rhi"])
    a2.set_ylim(pr["glo"], pr["ghi"])
    a2.set_xlabel("radius  $r$  [mm]")
    a2.set_ylabel(r"$|\mathrm{d}f/\mathrm{d}r|$   [1/mm]")

    s0, s1 = pr["sc"]
    a2.text(0.022, 0.95,
            "$\\delta_\\omega$   %.0f  /  %.0f $\\mu$m\n"
            "$\\theta_c$   %.0f  /  %.0f $\\mu$m\n"
            "$\\delta^*_c$   %.0f  /  %.0f $\\mu$m\n"
            "$H$    %.2f  /  %.2f"
            % (s0["dw"]*1e6, s1["dw"]*1e6, s0["theta_c"]*1e6, s1["theta_c"]*1e6,
               s0["dstar_c"]*1e6, s1["dstar_c"]*1e6, s0["H"], s1["H"]),
            transform=a2.transAxes, ha="left", va="top", fontsize=9.2,
            color="0.25", linespacing=1.55)

    fig.text(0.5, 0.955 - row*0.487, pr["title"], ha="center", fontsize=11.5)

fig.savefig(HERE/"matched_pairs.png", dpi=250)
fig.savefig(HERE/"matched_pairs.pdf")
print("\nwrote matched_pairs.png/pdf")
