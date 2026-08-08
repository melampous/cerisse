#!/usr/bin/env python3
"""x_fc against the prescribed displacement thickness, one colour per grid.

Three grids, identical in everything except the refinement of the nozzle lip:

  L3              198.4 um throughout - the resolution the 3D cases actually
                  carry where x_fc is measured, their level 5 being an annulus
                  that stops at x = 0.75 D_e
  L4              lip 99.2 um, the lip resolution of the 3D reference cases
  L5              lip 49.6 um, refinement forced on the lip band and the axis
  L5 lip-only     lip 49.6 um, the same lip band, axis refinement removed

Every profile appears once per grid it has been run on, so a vertical
separation between two colours at the same abscissa is a pure grid effect and
a horizontal trend within one colour is the physical response to the inlet.

The top hat is the discontinuous limit, delta* = 0: it prescribes no shear
layer at all, so it is the case for which the grid, not the boundary
condition, sets the lip gradient. It is drawn as a separate marker.

x_fc uses the thesis extraction: maximum raw density gradient between the
first expansion minimum and the first downstream principal maximum, refined by
a local quadratic fit, on the time-averaged centreline density of the finest
available level. Error bars are the standard error over seven disjoint
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
OUT = HERE/"sweep_out"
RHO_J, T0 = 1.6413, 4.0e-4
C3, M3 = 0.9443, -0.0766/1000.0
EXP, EXP_LO, EXP_HI = 0.888, 0.848, 0.928

# name, a[um], s, form
PROF = [("A0_tophat",     0.0, 3.0, "shifted"), ("Shift_s0",   170.0, 0.0, "shifted"),
        ("Wall100",     100.0, 3.0, "wall"),    ("C1_wall200", 200.0, 3.0, "wall"),
        ("Shift_s1",    170.0, 1.0, "shifted"), ("Wall340",    340.0, 3.0, "wall"),
        ("A1_a100",     100.0, 3.0, "shifted"), ("Shift_s2",   170.0, 2.0, "shifted"),
        ("Wall470",     470.0, 3.0, "wall"),    ("B1_wall680", 679.5, 3.0, "wall"),
        ("A2_a170",     170.0, 3.0, "shifted"), ("S4_a170_s4", 170.0, 4.0, "shifted"),
        ("Wall900",     900.0, 3.0, "wall"),    ("A3_a254",    254.0, 3.0, "shifted"),
        ("S5_a170_s5",  170.0, 5.0, "shifted"), ("Wall1150",  1150.0, 3.0, "wall"),
        ("S6_a170_s6",  170.0, 6.0, "shifted"), ("S7_a170_s7", 170.0, 7.0, "shifted")]

# label, colour, marker, filename patterns tried in order
GRIDS = [("L4   lip 99.2 $\\mu$m",              "#0072B2", "o",
          ("XFC_%s__L4.npz", "XFC_%s_L4.npz")),
         ("L5   lip 49.6 $\\mu$m",              "#D55E00", "s",
          ("XFC_%s__L5.npz", "XFC_%s.npz")),
         ("L5   lip 49.6 $\\mu$m, no axis refinement", "#009E73", "^",
          ("XFC_%s_LO__L5.npz", "XFC_%s_LO.npz")),
         ("L3   198.4 $\\mu$m everywhere", "#CC79A7", "v",
          ("XFC_%s_LO__L3.npz", "XFC_%s__L3.npz"))]

THREE_D = [("top hat", 0.0, 3.0, "shifted", 0.945),
           ("wall $a$=200", 200.0, 3.0, "wall", 0.930),
           ("shifted $a$=100", 100.0, 3.0, "shifted", 0.922),
           ("shifted $a$=254", 254.0, 3.0, "shifted", 0.883)]


def refine(x, y, i):
    if i <= 0 or i >= len(y)-1:
        return x[i]
    den = y[i-1]-2*y[i]+y[i+1]
    return x[i] if den == 0 else x[i]+0.5*(y[i-1]-y[i+1])/den*(x[1]-x[0])


def xfc(x, rho):
    e = np.where((x >= 0.35) & (x <= 1.0))[0]
    i0 = e[np.nanargmin(rho[e])]
    c = np.where((np.arange(x.size) > i0) & (x <= 1.7))[0]
    i1 = c[np.nanargmax(rho[c])]
    g = np.gradient(rho, x)
    s = np.arange(i0, i1+1)
    return refine(x, g, s[np.nanargmax(g[s])])


def measure(pats, nm):
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
    t = np.array([float(d[k+"_t"][0]) for k in ks])
    M = np.array([d[k+"_rho"] for k in ks])/RHO_J
    x = (np.arange(M.shape[1])+0.5)*dx/D
    W = t-T0
    sub = [(W[i+1]*M[i+1]-W[i]*M[i])/(t[i+1]-t[i])
           for i in range(len(t)-1) if t[i+1]-t[i] > 1.5e-4]
    v = np.array([xfc(x, q) for q in sub])
    return xfc(x, M[-1]), (np.nanstd(v, ddof=1)/np.sqrt(len(v)) if len(v) > 1
                           else np.nan)


DS = {nm: scalars(a*1e-6, f, s)["dstar_c"]*1e6 for nm, a, s, f in PROF}
data = {}
for lab, col, mk, pats in GRIDS:
    pts = []
    for nm, a, s, f in PROF:
        v, e = measure(pats, nm)
        if v is not None:
            pts.append((DS[nm], nm, v, e))
    data[lab] = sorted(pts)

print("%-12s %8s | %s" % ("case", "d*[um]",
                          " | ".join("%-17s" % l.split()[0] for l, *_ in GRIDS)))
for nm, a, s, f in sorted(PROF, key=lambda q: DS[q[0]]):
    cells = []
    for lab, *_ in GRIDS:
        m = [p for p in data[lab] if p[1] == nm]
        cells.append("%.4f +- %.4f" % (m[0][2], m[0][3]) if m else "%-17s" % "-")
    print("%-12s %8.1f | %s" % (nm, DS[nm], " | ".join(cells)))
for lab, *_ in GRIDS:
    P = data[lab]
    if len(P) < 3:
        print("%-40s n=%d" % (lab, len(P)))
        continue
    X = np.array([p[0] for p in P]); Y = np.array([p[2] for p in P])
    m, c = np.polyfit(X, Y, 1)
    r = Y-(C3+M3*X)
    print("%-40s n=%2d  slope %+.4f/1000um  rms from 3D line %.4f"
          % (lab, len(P), m*1000, np.sqrt((r**2).mean())))

# ------------------------------------------------------------------ figure
plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})
fig, ax = plt.subplots(figsize=(7.4, 4.7))

ax.axhline(EXP, color="0.55", lw=1.0, zorder=1)
for y in (EXP_LO, EXP_HI):
    ax.axhline(y, color="0.82", lw=0.6, ls=":", zorder=1)
ax.text(1292, EXP+0.0020, "experiment  %.3f" % EXP, ha="right", fontsize=7.2,
        color="0.35")

xx = np.linspace(-40, 1300, 40)
ax.plot(xx, C3+M3*xx, color="0.25", lw=1.2, ls="--", zorder=2)
ax.plot([q[1] for q in
         [(0, scalars(a*1e-6, f, s)["dstar_c"]*1e6) for _, a, s, f, _ in THREE_D]],
        [v for *_, v in THREE_D], ls="", marker="*", color="k", ms=13, zorder=7)
ax.text(1120, C3+M3*1120+0.0078, "fit through the 3D cases", fontsize=7,
        rotation=-8.5, color="0.25", ha="center")

for lab, col, mk, _ in GRIDS:
    P = [p for p in data[lab] if p[1] != "A0_tophat"]
    T = [p for p in data[lab] if p[1] == "A0_tophat"]
    if P:
        ax.plot([p[0] for p in P], [p[2] for p in P], color=col, lw=1.1,
                alpha=0.5, zorder=3)
        ax.errorbar([p[0] for p in P], [p[2] for p in P],
                    yerr=[p[3] for p in P], fmt=mk, color=col, ms=5.6,
                    mec="w", mew=0.7, capsize=2.4, lw=1.0, ls="", zorder=5)
    for p in T:                                   # top hat: open, ringed, larger
        ax.errorbar(p[0], p[2], yerr=p[3], fmt=mk, color=col, ms=10.5,
                    mfc="w", mec=col, mew=2.0, capsize=3.0, lw=1.2, zorder=6)

th = sorted(p[2] for lab in data for p in data[lab] if p[1] == "A0_tophat")
if len(th) >= 2:
    ax.plot([0.0, 0.0], [th[0], th[-1]], color="0.45", lw=0.9, ls="-",
            zorder=4)
    ax.annotate("", xy=(0, th[0]), xytext=(0, th[-1]),
                arrowprops=dict(arrowstyle="<->", lw=0.9, color="0.45",
                                shrinkA=9, shrinkB=9), zorder=4)
    ax.annotate("top hat,  $\\delta^*$ = 0 - the discontinuous limit.\n"
                "No shear layer is prescribed, so the lip gradient\n"
                "is set by the mesh: halving the cell moves $x_{fc}$\n"
                "by %.3f $D_e$, ten times the sampling error."
                % (th[-1]-th[0]),
                xy=(210, 0.8605), xytext=(210, 0.8605),
                textcoords="data", fontsize=7.4, color="0.2", va="bottom")

for nm, a, s, f, v in THREE_D:
    ax.annotate(nm, (scalars(a*1e-6, f, s)["dstar_c"]*1e6, v),
                textcoords="offset points", xytext=(8, 4), fontsize=7)

ax.set_xlabel(r"prescribed displacement thickness  $\delta^*_c$  [$\mu$m]")
ax.set_ylabel("$x_{fc}/D_e$")
ax.set_xlim(-55, 1300)
ax.set_ylim(0.845, 0.952)
ax.grid(True, ls="--", lw=0.5, color="0.93")
ax.set_axisbelow(True)

H = [Line2D([], [], ls="", marker="*", color="k", ms=12,
            label="3D reference,  lip 99.2 $\\mu$m")]
H += [Line2D([], [], color=c, marker=m, ms=5.6, mec="w", mew=0.7,
             label="2D  %s   (%d cases)" % (l, len(data[l])))
      for l, c, m, _ in GRIDS]
H += [Line2D([], [], ls="", marker="o", color="0.35", mfc="w", ms=9, mew=2.0,
             label="open ring: top hat, $\\delta^*$ = 0"),
      Line2D([], [], color="0.25", lw=1.2, ls="--",
             label="correlation fitted to the 3D cases")]
ax.legend(handles=H, fontsize=7.3, frameon=False, loc="upper right",
          handlelength=1.8, labelspacing=0.38, borderaxespad=0.8)
fig.tight_layout()
fig.savefig(HERE/"xfc_by_grid.png", dpi=270)
fig.savefig(HERE/"xfc_by_grid.pdf")
print("\nwrote xfc_by_grid.png/pdf")
