#!/usr/bin/env python3
"""What the four inlet parameters are, drawn on the profiles themselves.

Left  : the prescribed velocity shape f = u/U_e and the normalised mass flux
        g = rho u / (rho_e U_e) near the lip, for four representative inlets.
Right : the four integral scalars for every profile in the matrix, so the
        ranges and the correlations between them are visible.

  delta*_c = int (1 - g) dr        mass-flux deficit
  theta_c  = int g (1 - f) dr      momentum-flux deficit
  delta_w  = 1 / max|df/dr|        steepness of the shear layer
  C_d      = (2/R^2) int g r dr    discharge coefficient; D_eff/D_e = sqrt(C_d)

The shaded band on the left is the area whose integral is delta*_c: the mass
the profile fails to carry compared with a uniform sonic exit.
"""
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from design_2d_sweep import scalars, fvel, D, R, U_E, T_AMB, T_E, CP, RGAS

HERE = Path(__file__).resolve().parent

SHOW = [("top hat",              0.0, 3.0, "shifted", "#111111"),
        ("wall tanh  $a$=340",  340.0, 3.0, "wall",    "#b3251e"),
        ("shifted  $a$=170, $s$=3", 170.0, 3.0, "shifted", "#1f4e9c"),
        ("shifted  $a$=170, $s$=7", 170.0, 7.0, "shifted", "#6a3d9a")]

ALL = [("top hat",      0.0, 3.0, "shifted", "hat"),
       ("Shift s=0",  170.0, 0.0, "shifted", "trunc"),
       ("wall 100",   100.0, 3.0, "wall",    "wall"),
       ("wall 200",   200.0, 3.0, "wall",    "wall"),
       ("Shift s=1",  170.0, 1.0, "shifted", "trunc"),
       ("wall 340",   340.0, 3.0, "wall",    "wall"),
       ("shift 100",  100.0, 3.0, "shifted", "shift"),
       ("Shift s=2",  170.0, 2.0, "shifted", "trunc"),
       ("wall 470",   470.0, 3.0, "wall",    "wall"),
       ("wall 680",   679.5, 3.0, "wall",    "wall"),
       ("shift 170",  170.0, 3.0, "shifted", "shift"),
       ("s=4",        170.0, 4.0, "shifted", "shift"),
       ("wall 900",   900.0, 3.0, "wall",    "wall"),
       ("shift 254",  254.0, 3.0, "shifted", "shift"),
       ("s=5",        170.0, 5.0, "shifted", "shift"),
       ("wall 1150", 1150.0, 3.0, "wall",    "wall"),
       ("s=6",        170.0, 6.0, "shifted", "shift"),
       ("s=7",        170.0, 7.0, "shifted", "shift")]
COL = {"wall": "#b3251e", "shift": "#1f4e9c", "trunc": "#6a3d9a",
       "hat": "#111111"}
MRK = {"wall": "s", "shift": "o", "trunc": "^", "hat": "P"}

print("%-12s %8s %9s %9s %9s %7s %9s"
      % ("profile", "a[um]", "d*_c[um]", "theta[um]", "dw[um]", "H", "sqrt(Cd)"))
rec = []
for nm, a, s, form, fam in ALL:
    sc = scalars(a*1e-6, form, s)
    rec.append((nm, fam, sc["dstar_c"]*1e6, sc["theta_c"]*1e6, sc["dw"]*1e6,
                sc["H"], np.sqrt(sc["Cd"])))
    print("%-12s %8.0f %9.1f %9.1f %9.1f %7s %9.4f"
          % (nm, a, rec[-1][2], rec[-1][3], rec[-1][4],
             "-" if not np.isfinite(sc["H"]) else "%.2f" % sc["H"], rec[-1][6]))

plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})
fig = plt.figure(figsize=(13.0, 5.4))
gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.0, 1.0], wspace=0.26,
                      left=0.055, right=0.99, top=0.93, bottom=0.115)
a1 = fig.add_subplot(gs[0, 0])
a2 = fig.add_subplot(gs[0, 1])
a3 = fig.add_subplot(gs[0, 2])

r = np.linspace(0.0, R, 400001)
rmm = (R-r)*1e3                                    # distance inside the lip, mm
for lab, a, s, form, col in SHOW:
    f = fvel(r, a*1e-6, form, s)
    u = U_E*f
    T = T_AMB - 0.5*u*u/CP
    g = (1.0/(RGAS*T))*u/((1.0/(RGAS*T_E))*U_E)
    m = rmm <= 4.0
    a1.plot(rmm[m], f[m], color=col, lw=1.5, label=lab)
    a1.plot(rmm[m], g[m], color=col, lw=1.0, ls=":")
    if a > 0:
        sc = scalars(a*1e-6, form, s)
        a1.fill_between(rmm[m], g[m], 1.0, color=col, alpha=0.07, lw=0)
a1.set_xlabel("distance inside the lip,  $R-r$  [mm]")
a1.set_ylabel("$f=u/U_e$   (solid),    $g=\\rho u/\\rho_e U_e$   (dotted)")
a1.set_xlim(0, 4.0)
a1.set_ylim(0, 1.06)
a1.grid(True, ls="--", lw=0.5, color="0.93")
a1.set_axisbelow(True)
a1.legend(fontsize=7.6, frameon=False, loc="lower right", handlelength=1.9)
a1.text(0.03, 0.955, "shaded area $=\\delta^*_c=\\int(1-g)\\,\\mathrm{d}r$",
        transform=a1.transAxes, va="top", fontsize=8, color="0.25")

for i, (xk, xl, yk, yl) in enumerate(
        [(2, r"$\delta^*_c$  [$\mu$m]", 3, r"$\theta_c$  [$\mu$m]"),
         (2, r"$\delta^*_c$  [$\mu$m]", 6, r"$D_{\rm eff}/D_e=\sqrt{C_d}$")]):
    ax = (a2, a3)[i]
    for nm, fam, ds, th, dw, H, sq in rec:
        ax.plot((ds, ds)[0], (th, sq)[i], MRK[fam], color=COL[fam],
                ms=7 if fam == "hat" else 5.4, mec="w", mew=0.7)
    X = np.array([q[2] for q in rec])
    Y = np.array([q[yk] for q in rec])
    m, c = np.polyfit(X, Y, 1)
    res = Y-(c+m*X)
    r2 = 1-(res**2).sum()/((Y-Y.mean())**2).sum()
    xx = np.linspace(0, X.max()*1.04, 20)
    ax.plot(xx, c+m*xx, color="0.45", lw=1.1, ls="--")
    ax.set_xlabel(xl)
    ax.set_ylabel(yl)
    ax.grid(True, ls="--", lw=0.5, color="0.93")
    ax.set_axisbelow(True)
    ax.text(0.04, 0.94 if i else 0.94, "$R^2$ = %.4f" % r2,
            transform=ax.transAxes, va="top", fontsize=8.4,
            bbox=dict(fc="w", ec="0.85", lw=0.5, pad=2.2))

a2.text(0.97, 0.06,
        "the wall family (red) and the shifted\nfamily (blue) separate: at the "
        "same\n$\\delta^*$ they carry different $\\theta$",
        transform=a2.transAxes, ha="right", va="bottom", fontsize=7.2,
        color="0.35")
a3.text(0.97, 0.94,
        "one line, all families:\n$\\sqrt{C_d}=1.0000-1.9888\\,\\delta^*_c/D_e$",
        transform=a3.transAxes, ha="right", va="top", fontsize=7.2,
        color="0.35")

fig.savefig(HERE/"param_definitions.png", dpi=250)
fig.savefig(HERE/"param_definitions.pdf")
print("\nwrote param_definitions.png/pdf")
