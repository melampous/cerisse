#!/usr/bin/env python3
"""What the red and the blue inlet profiles actually are.

Same presentation as the profile overview of Section 3.3.2: the prescribed
velocity shape f = u/U_e against radius, zoomed on the nozzle lip, with the
finite radial shear beside it.

Red   wall-attached tanh,  f = tanh(max(R-r,0)/a),  a = 200 to 1150 um. The shape factor
      is 2.26 for every member, so the family walks delta* and theta together.
Blue  shifted tanh, f = 1/2 [1 - tanh((r-r0)/a)] with r0 = R - s a, at a fixed
      a = 170 um and s = 3 to 7. Shifting the same tanh outward leaves its slope
      untouched, so this family has one momentum thickness (75.0 - 75.2 um) and
      one vorticity thickness (340 um) throughout while delta* runs from 533 to
      1213 um. It is the series that separates the mass deficit from the
      momentum deficit.

Layout follows Figure 3.3 of Section 3.3.2: full radial profile, lip region,
finite radial shear, with the region outside the nozzle lip shaded in each
panel. No text is placed inside the axes; the key is below the figure.
"""
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
from design_2d_sweep import scalars, fvel, R, D

HERE = Path(__file__).resolve().parent

WALL = [200.0, 340.0, 470.0, 679.5, 900.0, 1150.0]
SHIFT_S = [3.0, 4.0, 5.0, 6.0, 7.0]
A_SHIFT = 170.0
REDS = plt.get_cmap("Reds")(np.linspace(0.42, 0.92, len(WALL)))
BLUES = plt.get_cmap("Blues")(np.linspace(0.42, 0.92, len(SHIFT_S)))

r = np.linspace(0.0, R, 400001)
rmm = r*1e3
Rmm = R*1e3

print("%-22s %8s %9s %9s %7s" % ("profile", "a[um]", "d*_c[um]", "theta_c", "H"))
CASES = []
for a, col in zip(WALL, REDS):
    sc = scalars(a*1e-6, "wall", 3.0)
    CASES.append(("wall  $\\delta_s$=%d" % a, a, 3.0, "wall", col, sc))
    print("%-22s %8.0f %9.1f %9.1f %7.2f"
          % ("wall a=%g" % a, a, sc["dstar_c"]*1e6, sc["theta_c"]*1e6, sc["H"]))
for s, col in zip(SHIFT_S, BLUES):
    sc = scalars(A_SHIFT*1e-6, "shifted", s)
    CASES.append(("shifted  $s$=%d" % s, A_SHIFT, s, "shifted", col, sc))
    print("%-22s %8.0f %9.1f %9.1f %7.2f"
          % ("shifted a=170 s=%g" % s, A_SHIFT, sc["dstar_c"]*1e6,
             sc["theta_c"]*1e6, sc["H"]))

plt.rcParams.update({"font.size": 11.2, "axes.linewidth": 0.7,
                     "axes.labelsize": 11.5,
                     "xtick.labelsize": 10.6,
                     "ytick.labelsize": 10.6,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True,
                     "legend.frameon": False})
fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(10.4, 3.55))
fig.subplots_adjust(left=0.060, right=0.988, top=0.925, bottom=0.205,
                    wspace=0.235)

F = [(lab, fvel(r, a*1e-6, form, s), col, form)
     for lab, a, s, form, col, sc in CASES]

# ------------------------------------------------- (a) full radial profile
a1.axhspan(Rmm, 20.0, color="0.915", zorder=0)
for lab, f, col, form in F:
    a1.plot(f, rmm, color=col, lw=1.6, zorder=4)
a1.axhline(Rmm, color="0.4", lw=0.9, zorder=2)
a1.set_xlabel("profile function  $f = u/U_e$")
a1.set_ylabel("radius  $r$  [mm]")
a1.set_xlim(-0.02, 1.06)
a1.set_ylim(0.0, 13.9)
a1.grid(True, ls="--", lw=0.5, color="0.93")
a1.set_axisbelow(True)
a1.set_title("(a)  full radial profile", fontsize=11.5, pad=5)

# ------------------------------------------------- (b) lip region
a2.axhspan(Rmm, 20.0, color="0.915", zorder=0)
for lab, f, col, form in F:
    a2.plot(f, rmm, color=col, lw=1.7, zorder=4)
a2.axhline(Rmm, color="0.4", lw=0.9, zorder=2)
a2.set_xlabel("profile function  $f = u/U_e$")
a2.set_ylabel("radius  $r$  [mm]")
a2.set_xlim(-0.02, 1.06)
a2.set_ylim(10.9, 13.15)
a2.grid(True, ls="--", lw=0.5, color="0.93")
a2.set_axisbelow(True)
a2.set_title("(b)  lip region", fontsize=11.5, pad=5)

# ------------------------------------------------- (c) finite radial shear
a3.axvspan(Rmm, 20.0, color="0.915", zorder=0)
for lab, f, col, form in F:
    a3.plot(rmm, np.abs(np.gradient(f, rmm)), color=col, lw=1.7, zorder=4)
a3.axvline(Rmm, color="0.4", lw=0.9, zorder=2)
a3.set_xlabel("radius  $r$  [mm]")
a3.set_ylabel(r"$|\mathrm{d}f/\mathrm{d}r|$   [1/mm]")
a3.set_xlim(10.9, 13.15)
a3.set_ylim(0, 5.6)
a3.grid(True, ls="--", lw=0.5, color="0.93")
a3.set_axisbelow(True)
a3.set_title("(c)  finite radial shear", fontsize=11.5, pad=5)

# keys placed in the empty parts of the panels: the wall family and the top hat
# in the lower left of (a), the shifted family in the upper left of (c)
HW = [Line2D([], [], color=c[4], lw=2.0, label="$\\delta_s$ = %g $\\mu$m" % c[1])
       for c in CASES if c[3] == "wall"]
lw_ = a1.legend(handles=HW, fontsize=10.3, loc="center left",
                bbox_to_anchor=(0.045, 0.40), handlelength=1.9,
                labelspacing=0.32, title="wall-attached tanh")
lw_.get_title().set_fontsize(10.6)

HS = [Line2D([], [], color=c[4], lw=2.0, label="$s$ = %d" % c[2])
      for c in CASES if c[3] == "shifted"]
ls_ = a3.legend(handles=HS, fontsize=9.4, loc="upper left",
                bbox_to_anchor=(0.025, 0.975), ncol=2, handlelength=1.5,
                columnspacing=0.8, labelspacing=0.22,
                title="shifted tanh,  $\\delta_s$ = 170 $\\mu$m")
ls_.get_title().set_fontsize(9.8)

fig.savefig(HERE/"kept_profiles.png", dpi=250)
fig.savefig(HERE/"kept_profiles.pdf")
print("\nwrote kept_profiles.png/pdf")
