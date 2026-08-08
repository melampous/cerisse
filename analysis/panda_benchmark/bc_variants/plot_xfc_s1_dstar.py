#!/usr/bin/env python3
"""x_fc and S_1 against the prescribed displacement thickness.

Upper panel  x_fc : the four 3D reference cases with the line fitted through
                    them, and the 2D sweep.
Lower panel  S_1  : the 3D cases only. The 2D spacings are not used - on the L4
                    grid S_1 correlates with the inlet parameter at R^2 = 0.50
                    against 0.86 for x_fc, so those points would not measure
                    what this figure is about.

The 2D set is the L4 sweep, except for the top hat, which is taken from L3.
The top hat is the one case whose lip gradient is set by the mesh rather than
by the boundary condition, so it moves 0.047 D_e across the 2D grids; L3
(198.4 um throughout) is the grid that matches the resolution the 3D cases
actually carry where x_fc is measured, and its top hat, 0.9432, is the closest
of the five 2D variants to the 3D value of 0.945.

Colours separate the two inlet families: the wall tanh, whose shape factor is
2.26 for every width, and the shifted (thin) tanh, whose shape factor runs from
1.4 to 14.0. Both are drawn against the same abscissa, so where they separate
vertically the displacement thickness alone is not sufficient.
"""
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
from design_2d_sweep import scalars, D

HERE = Path(__file__).resolve().parent
OUT = HERE/"sweep_out"
RHO_J, T0 = 1.6413, 4.0e-4
CWALL, CTHIN, CHAT = "#b3251e", "#1f4e9c", "#111111"
EXP_XFC, EXP_S1 = 0.888, 1.359
EXP_LO, EXP_HI = 0.848, 0.928

# name, a[um], s, form, family, label, grid to read
PROF = [("A0_tophat",     0.0, 3.0, "shifted", "hat",   "top hat",          "L3"),
        ("Shift_s0",    170.0, 0.0, "shifted", "trunc", "shifted 170, s=0", "L4"),
        ("Wall100",     100.0, 3.0, "wall",    "wall",  "wall 100",         "L4"),
        ("C1_wall200",  200.0, 3.0, "wall",    "wall",  "wall 200",         "L4"),
        ("Shift_s1",    170.0, 1.0, "shifted", "trunc", "shifted 170, s=1", "L4"),
        ("Wall340",     340.0, 3.0, "wall",    "wall",  "wall 340",         "L4"),
        ("A1_a100",     100.0, 3.0, "shifted", "thin",  "shifted 100",      "L4"),
        ("Shift_s2",    170.0, 2.0, "shifted", "trunc", "shifted 170, s=2", "L4"),
        ("Wall470",     470.0, 3.0, "wall",    "wall",  "wall 470",         "L4"),
        ("B1_wall680",  679.5, 3.0, "wall",    "wall",  "wall 680",         "L4"),
        ("A2_a170",     170.0, 3.0, "shifted", "thin",  "shifted 170, s=3", "L4"),
        ("S4_a170_s4",  170.0, 4.0, "shifted", "thin",  "shifted 170, s=4", "L4"),
        ("Wall900",     900.0, 3.0, "wall",    "wall",  "wall 900",         "L4"),
        ("A3_a254",     254.0, 3.0, "shifted", "thin",  "shifted 254",      "L4"),
        ("S5_a170_s5",  170.0, 5.0, "shifted", "thin",  "shifted 170, s=5", "L4"),
        ("Wall1150",   1150.0, 3.0, "wall",    "wall",  "wall 1150",        "L4"),
        ("S6_a170_s6",  170.0, 6.0, "shifted", "thin",  "shifted 170, s=6", "L4"),
        ("S7_a170_s7",  170.0, 7.0, "shifted", "thin",  "shifted 170, s=7", "L4")]
COL = {"wall": CWALL, "thin": CTHIN, "trunc": CTHIN, "hat": CHAT}
MRK = {"wall": "s", "thin": "o", "trunc": "^", "hat": "o"}

THREE = [("top hat",       0.0, 3.0, "shifted", 0.945, 1.244),
         ("wall 200",    200.0, 3.0, "wall",    0.930, 1.415),
         ("shifted 100", 100.0, 3.0, "shifted", 0.922, 1.508),
         ("shifted 254", 254.0, 3.0, "shifted", 0.883, 1.300)]


def refine(x, y, i):
    if i <= 0 or i >= len(y)-1:
        return x[i]
    d = y[i-1]-2*y[i]+y[i+1]
    return x[i] if d == 0 else x[i]+0.5*(y[i-1]-y[i+1])/d*(x[1]-x[0])


def xfc_of(x, rho):
    e = np.where((x >= 0.35) & (x <= 1.0))[0]
    i0 = e[np.nanargmin(rho[e])]
    c = np.where((np.arange(x.size) > i0) & (x <= 1.7))[0]
    i1 = c[np.nanargmax(rho[c])]
    g = np.gradient(rho, x)
    s = np.arange(i0, i1+1)
    return refine(x, g, s[np.nanargmax(g[s])])


def load(nm, grid):
    pats = (("XFC_%s__L4.npz", "XFC_%s_L4.npz") if grid == "L4"
            else ("XFC_%s_LO__L3.npz", "XFC_%s__L3.npz"))
    for pat in pats:
        f = OUT/(pat % nm)
        if f.exists():
            break
    else:
        return None, None
    d = np.load(f)
    ks = sorted([k[:-2] for k in d.files if k.endswith("_t")],
                key=lambda q: int("".join(c for c in q if c.isdigit())))
    dx = float(d[ks[0]+"_t"][1])
    M = np.array([d[k+"_rho"] for k in ks])/RHO_J
    x = (np.arange(M.shape[1])+0.5)*dx/D
    cum = np.array([xfc_of(x, m) for m in M])
    full = cum[-1]
    err = max(np.abs(cum[len(cum)//2:]-full).max(), 0.5*dx/D)
    return full, err


rows = []
for nm, a, s, form, fam, lab, grid in PROF:
    v, e = load(nm, grid)
    if v is None:
        continue
    sc = scalars(a*1e-6, form, s)
    rows.append(dict(nm=nm, fam=fam, lab=lab, grid=grid, ds=sc["dstar_c"]*1e6,
                     H=sc["H"], v=v, e=e))
rows.sort(key=lambda q: q["ds"])

X3 = np.array([scalars(a*1e-6, f, s)["dstar_c"]*1e6 for _, a, s, f, *_ in THREE])
Y3f = np.array([q[4] for q in THREE])
Y3s = np.array([q[5] for q in THREE])


def fit(x, y):
    m, c = np.polyfit(x, y, 1)
    r = y-(c+m*x)
    return m, c, 1-(r**2).sum()/((y-y.mean())**2).sum()


mf, cf, r2f = fit(X3, Y3f)
ms, cs, r2s = fit(X3, Y3s)
X2 = np.array([r["ds"] for r in rows])
Y2 = np.array([r["v"] for r in rows])
m2, c2, r22 = fit(X2, Y2)

print("3D  x_fc = %.4f %+.4f d*/1000um   R2 %.4f" % (cf, mf*1000, r2f))
print("3D  S_1  = %.4f %+.4f d*/1000um   R2 %.4f" % (cs, ms*1000, r2s))
print("2D  x_fc = %.4f %+.4f d*/1000um   R2 %.4f   (n=%d)"
      % (c2, m2*1000, r22, len(rows)))
print("\n%-20s %6s %8s %6s %9s %9s %9s"
      % ("profile", "grid", "d*_c", "H", "x_fc", "+-", "resid 3D"))
for r in rows:
    print("%-20s %6s %8.1f %6s %9.4f %9.4f %+9.4f"
          % (r["lab"], r["grid"], r["ds"],
             "-" if not np.isfinite(r["H"]) else "%.2f" % r["H"],
             r["v"], r["e"], r["v"]-(cf+mf*r["ds"])))

plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})
fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.6, 8.2), sharex=True,
                             gridspec_kw={"hspace": 0.08,
                                          "height_ratios": [1.18, 1]})
xx = np.linspace(-45, 1290, 40)

# ------------------------------------------------------------- x_fc
a1.axhspan(EXP_LO, EXP_HI, color="#f4f4f4", zorder=0)
a1.axhline(EXP_XFC, color="0.5", lw=1.0, zorder=1)
a1.text(1282, EXP_XFC+0.0018, "experiment  %.3f" % EXP_XFC, ha="right",
        fontsize=7.4, color="0.32")
a1.text(1282, EXP_HI-0.0042, "experimental sampling interval", ha="right",
        fontsize=6.8, color="0.52")
a1.plot(xx, cf+mf*xx, color="k", lw=1.4, ls="--", zorder=3)
a1.text(660, 0.8665,
        "3D fit:   $x_{fc}$ = %.4f $-$ %.4f $\\delta^*_c$/1000$\\,\\mu$m\n"
        "          $R^2$ = %.4f" % (cf, -mf*1000, r2f),
        fontsize=7.6, color="0.15", ha="left", va="center",
        bbox=dict(fc="w", ec="0.86", lw=0.5, pad=3.0))
for r in rows:
    hat = r["fam"] == "hat"
    a1.errorbar(r["ds"], r["v"], yerr=r["e"], fmt=MRK[r["fam"]],
                color=COL[r["fam"]], ms=15.0 if hat else 6.2,
                mfc="none" if hat else COL[r["fam"]],
                mec=COL[r["fam"]] if hat else "w", mew=1.8 if hat else 0.8,
                capsize=2.6, lw=1.0, ls="", zorder=9 if hat else 5)
a1.plot(X3, Y3f, "*", color="k", ms=16, zorder=8)
for (nm, a, s, f, v, _), xv in zip(THREE, X3):
    a1.annotate(nm, (xv, v), textcoords="offset points", xytext=(10, 5),
                fontsize=7.4)
a1.annotate("2D top hat on L3 = %.4f, versus 3D 0.945:\n"
            "L3 is the 2D grid whose resolution at $x_{fc}$\n"
            "matches the 3D cases, and it is the closest\n"
            "of the five 2D grids." % rows[0]["v"],
            xy=(24, rows[0]["v"]), xytext=(300, 0.9395), textcoords="data",
            fontsize=7.2, color="0.22", va="center",
            arrowprops=dict(arrowstyle="->", lw=0.8, color="0.5",
                            shrinkA=3, shrinkB=9))
a1.set_ylabel("$x_{fc}/D_e$")
a1.set_ylim(0.845, 0.955)
a1.grid(True, ls="--", lw=0.5, color="0.93")
a1.set_axisbelow(True)
a1.legend(handles=[
    Line2D([], [], ls="", marker="*", color="k", ms=14,
           label="3D reference (4 cases)"),
    Line2D([], [], color="k", lw=1.4, ls="--", label="fit through the 3D cases"),
    Line2D([], [], ls="", marker="o", color=CHAT, ms=11, mfc="none", mew=1.8,
           label="2D top hat (L3)"),
    Line2D([], [], ls="", marker="s", color=CWALL, ms=6.2, mec="w",
           label="2D wall tanh, $H$ = 2.26  (7)"),
    Line2D([], [], ls="", marker="o", color=CTHIN, ms=6.2, mec="w",
           label="2D thin (shifted) tanh, $H$ = 4.1 - 14.0  (7)"),
    Line2D([], [], ls="", marker="^", color=CTHIN, ms=6.2, mec="w",
           label="2D shifted, truncated  $s$ = 0, 1, 2  (3)")],
    fontsize=7.2, frameon=False, loc="lower left", handlelength=1.7,
    labelspacing=0.34, borderaxespad=0.7)

# ------------------------------------------------------------- S_1
a2.axhline(EXP_S1, color="0.5", lw=1.0, zorder=1)
a2.text(1282, EXP_S1+0.004, "experiment  %.3f" % EXP_S1, ha="right",
        fontsize=7.4, color="0.32")
a2.plot(xx, cs+ms*xx, color="k", lw=1.4, ls="--", zorder=3)
a2.text(560, 1.386,
        "3D fit:   $S_1$ = %.4f $-$ %.4f $\\delta^*_c$/1000$\\,\\mu$m\n"
        "          $R^2$ = %.4f" % (cs, -ms*1000, r2s),
        fontsize=7.6, color="0.15", ha="left", va="center",
        bbox=dict(fc="w", ec="0.86", lw=0.5, pad=3.0))
a2.plot(X3, Y3s, "*", color="k", ms=16, zorder=8)
for (nm, a, s, f, _, v), xv in zip(THREE, X3):
    a2.annotate(nm, (xv, v), textcoords="offset points", xytext=(10, 5),
                fontsize=7.4)
a2.set_xlabel(r"prescribed displacement thickness  $\delta^*_c$  [$\mu$m]")
a2.set_ylabel("$S_1/D_e$")
a2.set_xlim(-50, 1290)
a2.grid(True, ls="--", lw=0.5, color="0.93")
a2.set_axisbelow(True)
a2.text(0.985, 0.10, "3D only: the 2D shock-cell spacings are not used\n"
        "($R^2$ = 0.50 against the inlet parameter, versus 0.86 for $x_{fc}$)",
        transform=a2.transAxes, ha="right", va="bottom", fontsize=7.2,
        color="0.35")

fig.savefig(HERE/"xfc_s1_dstar.png", dpi=260, bbox_inches="tight")
fig.savefig(HERE/"xfc_s1_dstar.pdf", bbox_inches="tight")
print("\nwrote xfc_s1_dstar.png/pdf")
