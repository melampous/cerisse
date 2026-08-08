#!/usr/bin/env python3
"""Does x_fc collapse onto the displacement thickness alone, or does the shape
of the inlet profile matter as well?

Eighteen inlet profiles, all run on the same grid (lip cells 99.2 um, the
resolution of the three-dimensional reference cases), form two ladders that
overlap in displacement thickness but differ in shape factor H = delta*/theta:

  wall-attached tanh   H = 2.26 throughout, seven points from 78 to 902 um;
  shifted tanh         H = 4.1 to 14.0, eight points from 314 to 1213 um;
  truncation series    s = 0, 1, 2 at fixed a = 170 um, H = 1.39, 2.41, 4.09.

If the two ladders lie on one curve, delta* is the single controlling scalar.
Where they separate, a second parameter is required - or the points are not
grid-converged, which the vertical markers at one and two lip cells identify.

x_fc uses the thesis extraction: the maximum raw density gradient between the
first expansion minimum and the first downstream principal maximum, refined by
a local quadratic fit, on the time-averaged centreline density of the finest
available level. Error bars are the standard error of seven disjoint
sub-windows of the 2.0 ms averaging interval.
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
RHO_J, T0 = 1.6413, 4.0e-4
C3, M3 = 0.9443, -0.0766 / 1000.0          # correlation fitted to the 3D cases
EXP, EXP_LO, EXP_HI = 0.888, 0.848, 0.928
DX_L4 = 99.2                               # um, lip cell of this grid

WALL = [("Wall100", 100.0), ("C1_wall200", 200.0), ("Wall340", 340.0),
        ("Wall470", 470.0), ("B1_wall680", 679.5), ("Wall900", 900.0),
        ("Wall1150", 1150.0)]
SHIFT = [("A1_a100", 100.0, 3.0), ("A2_a170", 170.0, 3.0), ("A3_a254", 254.0, 3.0),
         ("S4_a170_s4", 170.0, 4.0), ("S5_a170_s5", 170.0, 5.0),
         ("S6_a170_s6", 170.0, 6.0), ("S7_a170_s7", 170.0, 7.0)]
TRUNC = [("Shift_s0", 170.0, 0.0), ("Shift_s1", 170.0, 1.0), ("Shift_s2", 170.0, 2.0)]
NEW = {"Wall100", "Wall340", "Wall470", "Wall900", "Wall1150",
       "Shift_s0", "Shift_s1", "Shift_s2"}
THREE_D = [("top hat", 0.0, 3.0, "shifted", 0.945),
           ("wall $a$=200", 200.0, 3.0, "wall", 0.930),
           ("shifted $a$=100", 100.0, 3.0, "shifted", 0.922),
           ("shifted $a$=254", 254.0, 3.0, "shifted", 0.883)]

CW, CS, CT = "#b3251e", "#1f4e9c", "#6a3d9a"


def refine(x, y, i):
    if i <= 0 or i >= len(y) - 1:
        return x[i]
    den = y[i-1] - 2*y[i] + y[i+1]
    return x[i] if den == 0 else x[i] + 0.5*(y[i-1]-y[i+1])/den*(x[1]-x[0])


def xfc(x, rho):
    e = np.where((x >= 0.35) & (x <= 1.0))[0]
    i0 = e[np.nanargmin(rho[e])]
    c = np.where((np.arange(x.size) > i0) & (x <= 1.7))[0]
    i1 = c[np.nanargmax(rho[c])]
    g = np.gradient(rho, x)
    s = np.arange(i0, i1+1)
    return refine(x, g, s[np.nanargmax(g[s])])


def measure(nm):
    """New cases were written as XFC_<name>__L4.npz, earlier ones as _L4.npz."""
    for fn in ("XFC_%s__L4.npz" % nm, "XFC_%s_L4.npz" % nm):
        f = HERE/"sweep_out"/fn
        if f.exists():
            break
    else:
        return None, None
    d = np.load(f)
    ks = sorted([k[:-2] for k in d.files if k.endswith("_t")],
                key=lambda q: int("".join(c for c in q if c.isdigit())))
    dx = float(d[ks[0]+"_t"][1])
    t = np.array([float(d[k+"_t"][0]) for k in ks])
    M = np.array([d[k+"_rho"] for k in ks])/RHO_J
    x = (np.arange(M.shape[1])+0.5)*dx/D
    W = t - T0
    sub = [(W[i+1]*M[i+1]-W[i]*M[i])/(t[i+1]-t[i])
           for i in range(len(t)-1) if t[i+1]-t[i] > 1.5e-4]
    v = np.array([xfc(x, q) for q in sub])
    return xfc(x, M[-1]), (np.nanstd(v, ddof=1)/np.sqrt(len(v)) if len(v) > 1
                           else np.nan)


def series(spec, form):
    out = []
    for row in spec:
        nm, a = row[0], row[1]
        s = row[2] if len(row) > 2 else 3.0
        sc = scalars(a*1e-6, form, s)
        v, e = measure(nm)
        if v is None:
            continue
        out.append(dict(nm=nm, ds=sc["dstar_c"]*1e6, dw=sc["dw"]*1e6,
                        H=sc["H"], v=v, e=e, new=nm in NEW,
                        resid=v-(C3+M3*sc["dstar_c"]*1e6)))
    return sorted(out, key=lambda q: q["ds"])


wall, shift, trunc = series(WALL, "wall"), series(SHIFT, "shifted"), \
    series(TRUNC, "shifted")
th, eh = measure("A0_tophat")

print("%-12s %8s %6s %8s %8s %9s %8s"
      % ("case", "d*[um]", "H", "x_fc", "sigma", "resid3D", "d*/dx"))
for grp, lab in ((wall, "wall"), (shift, "shifted"), (trunc, "trunc")):
    for q in grp:
        print("%-12s %8.1f %6.2f %8.4f %8.4f %+9.4f %8.2f  %s"
              % (q["nm"], q["ds"], q["H"], q["v"], q["e"], q["resid"],
                 q["ds"]/DX_L4, "NEW" if q["new"] else ""))

for grp, lab in ((wall, "wall  H=2.26"), (shift, "shifted  H>=6")):
    X = np.array([q["ds"] for q in grp]); Y = np.array([q["v"] for q in grp])
    m, c = np.polyfit(X, Y, 1)
    r2 = 1-((Y-(m*X+c))**2).sum()/((Y-Y.mean())**2).sum()
    print("%-16s n=%d  slope %+.4f per 1000um  R2 %.3f" % (lab, len(X), m*1000, r2))
    Xt = X[X >= 600]; Yt = Y[X >= 600]
    if len(Xt) >= 2:
        mt, ct = np.polyfit(Xt, Yt, 1)
        print("%-16s thick end only (d*>=600): slope %+.4f  n=%d"
              % ("", mt*1000, len(Xt)))

# ------------------------------------------------------------------ figure
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})
fig, ax = plt.subplots(figsize=(7.6, 4.9))

ax.axhspan(EXP_LO, EXP_HI, color="#f2f2f2", zorder=0)
for y in (EXP_LO, EXP_HI):
    ax.axhline(y, color="0.78", lw=0.6, ls="-", zorder=1)
ax.axhline(EXP, color="0.45", lw=1.0, zorder=1)
ax.text(1288, EXP+0.0018, "experiment  %.3f" % EXP, ha="right", fontsize=7,
        color="0.3")
ax.text(1288, EXP_HI-0.0038, "experimental sampling interval", ha="right",
        fontsize=6.4, color="0.55")

xx = np.linspace(-40, 1300, 40)
ax.plot(xx, C3+M3*xx, color="k", lw=1.3, ls="--", zorder=2)
ax.text(1010, C3+M3*1010+0.0055,
        "correlation fitted to the 3D cases\n"
        "$x_{fc}$ = 0.9443 $-$ 0.0766 $\\delta^*_c$/1000$\\,\\mu$m",
        fontsize=6.8, rotation=-8.5, color="0.15", ha="center")

for n, lab in ((1, "1 lip cell"), (2, "2 lip cells")):
    ax.axvline(n*DX_L4, color="0.75", lw=0.7, ls=":", zorder=1)
    ax.text(n*DX_L4+7, 0.9485, lab, rotation=90, fontsize=6.2, color="0.5",
            va="top")

for grp, col, lab, mk in ((wall, CW, "wall tanh   $H$ = 2.26", "s"),
                          (shift, CS, "shifted tanh   $H$ = 6.0 - 14.0", "o"),
                          (trunc, CT, "truncation series  $s$ = 0, 1, 2", "^")):
    X = [q["ds"] for q in grp]; Y = [q["v"] for q in grp]
    ax.plot(X, Y, color=col, lw=1.0, alpha=0.55, zorder=3)
    for q in grp:
        ax.errorbar(q["ds"], q["v"], yerr=q["e"], fmt=mk, color=col,
                    ms=6.2 if q["new"] else 5.0,
                    mfc=col if q["new"] else "w", mec=col, mew=1.1,
                    capsize=2.4, lw=1.0, zorder=5)
if th:
    ax.errorbar(0.0, th, yerr=eh, fmt="P", color="#333333", ms=7.5, mec="w",
                mew=0.7, capsize=2.4, lw=1.0, zorder=6)
    ax.annotate("top hat (2D)", (0.0, th), textcoords="offset points",
                xytext=(9, -4), fontsize=6.8, color="#333333")

for nm, a, s, f, v in THREE_D:
    x = scalars(a*1e-6, f, s)["dstar_c"]*1e6
    ax.plot(x, v, "*", color="k", ms=13, zorder=7)
    ax.annotate(nm, (x, v), textcoords="offset points", xytext=(7, 5),
                fontsize=6.6)

for q in wall+shift+trunc:
    if q["nm"] in ("Wall100", "Wall340", "B1_wall680", "A2_a170", "Wall1150",
                   "Shift_s0"):
        off = {"Wall100": (8, -3), "Wall340": (-44, 1), "B1_wall680": (-62, -4),
               "A2_a170": (8, 4), "Wall1150": (6, -11),
               "Shift_s0": (-32, -11)}[q["nm"]]
        ax.annotate(q["nm"], (q["ds"], q["v"]), textcoords="offset points",
                    xytext=off, fontsize=6.4, color="0.25")

ax.set_xlabel(r"prescribed displacement thickness  $\delta^*_c$  [$\mu$m]")
ax.set_ylabel("$x_{fc}/D_e$")
ax.set_xlim(-55, 1300)
ax.set_ylim(0.848, 0.950)
ax.grid(True, ls="--", lw=0.5, color="0.93")
ax.set_axisbelow(True)
ax.legend(handles=[
    Line2D([], [], ls="", marker="*", color="k", ms=12,
           label="3D reference, lip 99.2 $\\mu$m"),
    Line2D([], [], color=CW, marker="s", ms=5.5,
           label="2D wall tanh,  $H$ = 2.26"),
    Line2D([], [], color=CS, marker="o", ms=5,
           label="2D shifted tanh,  $H$ = 6.0 - 14.0"),
    Line2D([], [], color=CT, marker="^", ms=5.5,
           label="2D truncation series,  $s$ = 0, 1, 2"),
    Line2D([], [], ls="", marker="s", color=CW, mfc="w", ms=5.5,
           label="open symbol: earlier case"),
    Line2D([], [], color="k", lw=1.3, ls="--", label="fit through the 3D cases")],
    fontsize=6.9, frameon=False, loc="lower left", handlelength=1.7,
    bbox_to_anchor=(0.005, 0.005),
    labelspacing=0.35)
fig.tight_layout()
fig.savefig(HERE/"two_ladders_L4.png", dpi=270)
fig.savefig(HERE/"two_ladders_L4.pdf")
print("\nwrote two_ladders_L4.png/pdf")
