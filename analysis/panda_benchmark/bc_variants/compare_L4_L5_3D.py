#!/usr/bin/env python3
"""Do the 2D cases fall on the correlation fitted to the 3D data once they are
run at the 3D lip resolution?

The 3D reference cases resolve the lip with 99.2 um cells. The first 2D sweep
used 49.6 um, which left the thin profiles on a different point of the grid
convergence curve; the top hat alone moved 0.031 D_e between the two. This
compares, on the same axes:

  - the four 3D cases and the correlation fitted through them,
  - the 2D sweep at 49.6 um,
  - the 2D sweep repeated at 99.2 um.

The figure of merit is the root-mean-square distance of the 2D points from the
3D correlation, evaluated only on the four profiles that exist in both.
Everything uses the thesis extraction: raw-gradient maximum between the first
expansion minimum and the first downstream principal maximum, refined by a
local quadratic fit.
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
RHO_J = 1.6413
T0 = 4.0e-4
C3D_M, C3D_C = -1.9461, 0.9443            # 3D correlation, x_fc vs delta*/D_e

# name, a[um], s, form, 3D counterpart x_fc (thesis table) or None
SPEC = [("A0_tophat",    0.0, 3.0, "shifted", 0.945),
        ("C1_wall200", 200.0, 3.0, "wall",    0.930),
        ("A1_a100",    100.0, 3.0, "shifted", 0.922),
        ("B1_wall680", 679.5, 3.0, "wall",    None),
        ("A2_a170",    170.0, 3.0, "shifted", None),
        ("S4_a170_s4", 170.0, 4.0, "shifted", None),
        ("A3_a254",    254.0, 3.0, "shifted", 0.883),
        ("S5_a170_s5", 170.0, 5.0, "shifted", None),
        ("S6_a170_s6", 170.0, 6.0, "shifted", None),
        ("S7_a170_s7", 170.0, 7.0, "shifted", None)]


def refine(x, y, i):
    if i <= 0 or i >= len(y)-1:
        return x[i]
    den = y[i-1] - 2*y[i] + y[i+1]
    return x[i] if den == 0 else x[i] + 0.5*(y[i-1]-y[i+1])/den*(x[1]-x[0])


def xfc_thesis(x, rho):
    e = np.where((x >= 0.35) & (x <= 1.0))[0]
    i0 = e[np.nanargmin(rho[e])]
    c = np.where((np.arange(x.size) > i0) & (x <= 1.7))[0]
    i1 = c[np.nanargmax(rho[c])]
    g = np.gradient(rho, x)
    s = np.arange(i0, i1+1)
    return refine(x, g, s[np.nanargmax(g[s])])


def measure(nm):
    f = HERE/"sweep_out"/("XFC_%s.npz" % nm)
    if not f.exists():
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
    v = np.array([xfc_thesis(x, q) for q in sub])
    return xfc_thesis(x, M[-1]), np.nanstd(v, ddof=1)/np.sqrt(len(v))


rows = []
print("%-13s %8s | %8s %7s | %8s %7s | %8s | %8s"
      % ("case", "d*[um]", "2D L5", "+-", "2D L4", "+-", "3D", "L4-L5"))
print("-"*86)
for nm, a, s, form, x3 in SPEC:
    sc = scalars(a*1e-6, form, s)
    v5, s5 = measure(nm)
    v4, s4 = measure(nm+"_L4")
    rows.append(dict(nm=nm, a=a, s=s, form=form, ds=sc["dstar_c"]*1e6,
                     v5=v5, s5=s5, v4=v4, s4=s4, x3=x3))
    print("%-13s %8.1f | %8s %7s | %8s %7s | %8s | %8s"
          % (nm, sc["dstar_c"]*1e6,
             "%.4f" % v5 if v5 else "-", "%.4f" % s5 if s5 else "-",
             "%.4f" % v4 if v4 else "-", "%.4f" % s4 if s4 else "-",
             "%.3f" % x3 if x3 else "-",
             "%+.4f" % (v4-v5) if (v4 and v5) else "-"))

pair = [r for r in rows if r["x3"] and r["v4"] and r["v5"]]
if pair:
    d5 = np.array([r["v5"]-r["x3"] for r in pair])
    d4 = np.array([r["v4"]-r["x3"] for r in pair])
    print("\ndistance from the matching 3D case, over %d shared profiles"
          % len(pair))
    print("   2D at 49.6 um : rms %.4f   max %.4f" % (np.sqrt((d5**2).mean()),
                                                      np.abs(d5).max()))
    print("   2D at 99.2 um : rms %.4f   max %.4f" % (np.sqrt((d4**2).mean()),
                                                      np.abs(d4).max()))

for tag, key in (("49.6 um", "v5"), ("99.2 um", "v4")):
    lad = sorted([r for r in rows if r["a"] == 170.0 and r[key]],
                 key=lambda r: r["s"])
    if len(lad) >= 3:
        X = np.array([r["ds"] for r in lad]); Y = np.array([r[key] for r in lad])
        m, c = np.polyfit(X, Y, 1)
        r2 = 1-((Y-(m*X+c))**2).sum()/((Y-Y.mean())**2).sum()
        print("controlled series at %s : slope %+.4f per 1000 um  R2 %.3f"
              % (tag, m*1000, r2))
print("3D correlation                : slope %+.4f per 1000 um"
      % (C3D_M/(D*1e6)*1000))    # per D_e -> per 1000 um

# ------------------------------------------------------------------ figure
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})
fig, ax = plt.subplots(figsize=(7.0, 4.4))
xx = np.linspace(0, 1300, 50)
ax.plot(xx, C3D_C + C3D_M*xx/(D*1e6), color="k", lw=1.3, ls="--", zorder=2,
        label="correlation fitted to the 3D cases")
for r in rows:
    if r["x3"]:
        ax.plot(r["ds"], r["x3"], "*", color="k", ms=13, zorder=6)
    if r["v5"]:
        ax.errorbar(r["ds"], r["v5"], yerr=r["s5"], fmt="o", color="#1f4e9c",
                    ms=5, mec="w", mew=0.6, capsize=2.5, lw=1.0, zorder=4)
    if r["v4"]:
        ax.errorbar(r["ds"], r["v4"], yerr=r["s4"], fmt="s", color="#b3251e",
                    ms=5.5, mec="w", mew=0.6, capsize=2.5, lw=1.0, zorder=5)
    if r["v5"] and r["v4"]:
        ax.annotate("", xy=(r["ds"], r["v4"]), xytext=(r["ds"], r["v5"]),
                    arrowprops=dict(arrowstyle="->", lw=0.8, color="0.55"),
                    zorder=3)
    if r["x3"]:
        ax.annotate(r["nm"].split("_")[0], (r["ds"], r["x3"]),
                    textcoords="offset points", xytext=(7, 4), fontsize=6.6)
ax.set_xlabel(r"displacement thickness  $\delta^*_c$  [$\mu$m]")
ax.set_ylabel("$x_{fc}/D_e$")
ax.set_xlim(-60, 1320)
ax.grid(True, ls="--", lw=0.5, color="0.9"); ax.set_axisbelow(True)
ax.legend(handles=[Line2D([], [], ls="", marker="*", color="k", ms=12,
                          label="3D reference (99.2 $\\mu$m lip)"),
                   Line2D([], [], ls="", marker="s", color="#b3251e", ms=6,
                          label="2D at 99.2 $\\mu$m lip"),
                   Line2D([], [], ls="", marker="o", color="#1f4e9c", ms=5,
                          label="2D at 49.6 $\\mu$m lip"),
                   Line2D([], [], color="k", lw=1.3, ls="--",
                          label="3D correlation")],
          fontsize=7.2, frameon=False, loc="upper right", handlelength=1.6)
fig.tight_layout()
fig.savefig(HERE/"compare_L4_L5_3D.png", dpi=260)
fig.savefig(HERE/"compare_L4_L5_3D.pdf")
print("\nwrote compare_L4_L5_3D.png/pdf")
