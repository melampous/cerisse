#!/usr/bin/env python3
"""Every x_fc measurement on one axis, against the correlation fitted to the
3D cases and the experimental interval.

Background:
  - the least-squares line through the four 3D cases of the thesis table,
    refitted here rather than quoted, so its coefficients and scatter are
    visible;
  - the experimental value and its sampling interval, drawn as a horizontal
    band because the experiment has no prescribed displacement thickness -
    where that band crosses the 3D line is the delta* the measurement implies.

Points: the 3D cases, the 2D sweep at 49.6 um, the same sweep repeated at
99.2 um (the 3D lip resolution), and the earlier 2D production run at 24.8 um.
All x_fc use the thesis extraction: raw-gradient maximum between the first
expansion minimum and the first downstream principal maximum, quadratic
refinement.
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
T0_NEW, T0_L6 = 4.0e-4, 0.01574728874

# 3D reference: thesis table values, with the prescribed profile of each
THREE_D = [("top hat",     0.0, 3.0, "shifted", 0.945),
           ("wall a=200",200.0, 3.0, "wall",    0.930),
           ("a=100",     100.0, 3.0, "shifted", 0.922),
           ("a=254",     254.0, 3.0, "shifted", 0.883)]
EXP, EXP_LO, EXP_HI = 0.888, 0.848, 0.928

SPEC = [("A0_tophat",    0.0, 3.0, "shifted"), ("C1_wall200", 200.0, 3.0, "wall"),
        ("A1_a100",    100.0, 3.0, "shifted"), ("B1_wall680", 679.5, 3.0, "wall"),
        ("A2_a170",    170.0, 3.0, "shifted"), ("S4_a170_s4", 170.0, 4.0, "shifted"),
        ("A3_a254",    254.0, 3.0, "shifted"), ("S5_a170_s5", 170.0, 5.0, "shifted"),
        ("S6_a170_s6", 170.0, 6.0, "shifted"), ("S7_a170_s7", 170.0, 7.0, "shifted")]


def refine(x, y, i):
    if i <= 0 or i >= len(y)-1:
        return x[i]
    d = y[i-1]-2*y[i]+y[i+1]
    return x[i] if d == 0 else x[i]+0.5*(y[i-1]-y[i+1])/d*(x[1]-x[0])


def xfc(x, r):
    e = np.where((x >= 0.35) & (x <= 1.0))[0]; i0 = e[np.nanargmin(r[e])]
    c = np.where((np.arange(x.size) > i0) & (x <= 1.7))[0]
    i1 = c[np.nanargmax(r[c])]
    g = np.gradient(r, x); s = np.arange(i0, i1+1)
    return refine(x, g, s[np.nanargmax(g[s])])


def measure(fn, t0):
    f = HERE/"sweep_out"/fn if not str(fn).startswith("/") else Path(fn)
    if not f.exists():
        f = HERE/fn
    if not f.exists():
        return None, None
    d = np.load(f)
    ks = sorted([k[:-2] for k in d.files if k.endswith("_t")],
                key=lambda q: int("".join(c for c in q if c.isdigit())))
    dx = float(d[ks[0]+"_t"][1]); t = np.array([float(d[k+"_t"][0]) for k in ks])
    M = np.array([d[k+"_rho"] for k in ks])/RHO_J
    x = (np.arange(M.shape[1])+0.5)*dx/D; W = t-t0
    sub = [(W[i+1]*M[i+1]-W[i]*M[i])/(t[i+1]-t[i])
           for i in range(len(t)-1) if t[i+1]-t[i] > 1.5e-4]
    v = np.array([xfc(x, q) for q in sub])
    return xfc(x, M[-1]), (np.nanstd(v, ddof=1)/np.sqrt(len(v)) if len(v) > 1 else np.nan)


# ---- refit the 3D correlation from its own four points -------------------
X3 = np.array([scalars(a*1e-6, f, s)["dstar_c"]*1e6 for _, a, s, f, _ in THREE_D])
Y3 = np.array([v for *_, v in THREE_D])
m3, c3 = np.polyfit(X3, Y3, 1)
r2 = 1-((Y3-(m3*X3+c3))**2).sum()/((Y3-Y3.mean())**2).sum()
print("3D correlation refitted: x_fc = %.4f %+.6f * d*[um]   (%+.4f per 1000 um)"
      "   R2=%.4f" % (c3, m3, m3*1000, r2))
print("experiment %.3f  ->  implied d* = %.0f um  (interval %.0f - %.0f)"
      % (EXP, (EXP-c3)/m3, (EXP_HI-c3)/m3, (EXP_LO-c3)/m3))

rows = []
for nm, a, s, form in SPEC:
    ds = scalars(a*1e-6, form, s)["dstar_c"]*1e6
    v5, e5 = measure("XFC_%s.npz" % nm, T0_NEW)
    v4, e4 = measure("XFC_%s_L4.npz" % nm, T0_NEW)
    rows.append((nm, ds, v5, e5, v4, e4))
v6, e6 = measure("XFC_window_L6.npz", T0_L6)
ds6 = scalars(254e-6, "shifted", 3.0)["dstar_c"]*1e6

print("\n%-13s %8s | %8s %8s %8s" % ("case", "d*[um]", "2D 49.6", "2D 99.2", "resid vs 3D line"))
for nm, ds, v5, e5, v4, e4 in rows:
    pred = c3+m3*ds
    print("%-13s %8.1f | %8s %8s   %8s / %8s"
          % (nm, ds, "%.4f" % v5 if v5 else "-", "%.4f" % v4 if v4 else "-",
             "%+.4f" % (v5-pred) if v5 else "-", "%+.4f" % (v4-pred) if v4 else "-"))
if v6:
    print("%-13s %8.1f | 24.8um %.4f              %+.4f" % ("m142_2d_L6stat", ds6, v6, v6-(c3+m3*ds6)))

for lbl, idx in (("49.6 um", 2), ("99.2 um", 4)):
    d = [r[idx]-(c3+m3*r[1]) for r in rows if r[idx]]
    if d:
        print("rms residual from the 3D line, 2D at %s : %.4f  (n=%d)"
              % (lbl, np.sqrt(np.mean(np.square(d))), len(d)))

# ------------------------------------------------------------------ figure
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})
fig, ax = plt.subplots(figsize=(7.4, 4.8))
xx = np.linspace(-40, 1320, 60)
ax.axhspan(EXP_LO, EXP_HI, color="0.86", zorder=0)
ax.axhline(EXP, color="0.35", lw=1.1, zorder=1)
ax.text(1300, EXP+0.0015, "experiment %.3f" % EXP, ha="right", fontsize=7,
        color="0.25")
ax.text(1300, EXP_HI-0.004, "experimental sampling interval", ha="right",
        fontsize=6.6, color="0.45")
ax.plot(xx, c3+m3*xx, color="k", lw=1.4, ls="--", zorder=2)
ax.text(60, c3+m3*60+0.006,
        "3D correlation  $x_{fc}$ = %.4f $-$ %.4f$\\,\\delta^*$/1000$\\mu$m   ($R^2$=%.3f)"
        % (c3, -m3*1000, r2), fontsize=7, rotation=-11)

for nm, ds, v5, e5, v4, e4 in rows:
    if v5 and v4:
        ax.annotate("", xy=(ds, v4), xytext=(ds, v5),
                    arrowprops=dict(arrowstyle="->", lw=0.8, color="0.6"), zorder=3)
    if v5:
        ax.errorbar(ds, v5, yerr=e5, fmt="o", color="#1f4e9c", ms=5, mec="w",
                    mew=0.6, capsize=2.5, lw=1.0, zorder=4)
    if v4:
        ax.errorbar(ds, v4, yerr=e4, fmt="s", color="#b3251e", ms=5.5, mec="w",
                    mew=0.6, capsize=2.5, lw=1.0, zorder=5)
if v6:
    ax.errorbar(ds6, v6, yerr=e6, fmt="D", color="#2ca02c", ms=6, mec="w",
                mew=0.6, capsize=2.5, lw=1.0, zorder=5)
for nm, a, s, f, v in THREE_D:
    ax.plot(scalars(a*1e-6, f, s)["dstar_c"]*1e6, v, "*", color="k", ms=14,
            zorder=6)
    ax.annotate(nm, (scalars(a*1e-6, f, s)["dstar_c"]*1e6, v),
                textcoords="offset points", xytext=(8, 5), fontsize=6.8)

ax.set_xlabel(r"prescribed displacement thickness  $\delta^*_c$  [$\mu$m]")
ax.set_ylabel("$x_{fc}/D_e$")
ax.set_xlim(-50, 1320)
ax.grid(True, ls="--", lw=0.5, color="0.92"); ax.set_axisbelow(True)
ax.legend(handles=[
    Line2D([], [], ls="", marker="*", color="k", ms=13, label="3D  lip 99.2 $\\mu$m"),
    Line2D([], [], ls="", marker="s", color="#b3251e", ms=6, label="2D  lip 99.2 $\\mu$m"),
    Line2D([], [], ls="", marker="o", color="#1f4e9c", ms=5, label="2D  lip 49.6 $\\mu$m"),
    Line2D([], [], ls="", marker="D", color="#2ca02c", ms=6, label="2D  lip 24.8 $\\mu$m (earlier run)"),
    Line2D([], [], color="k", lw=1.4, ls="--", label="fit through the 3D cases")],
    fontsize=7, frameon=False, loc="lower left", handlelength=1.6)
fig.tight_layout()
fig.savefig(HERE/"master_xfc.png", dpi=260)
fig.savefig(HERE/"master_xfc.pdf")
print("\nwrote master_xfc.png/pdf")
