#!/usr/bin/env python3
"""The payoff: does x_fc follow the effective aperture, or the shear-layer
shape?

Reads one npz per case (centreline DensityMEAN from every statistics dump,
produced by collect_sweep.sh) and answers three questions:

  1. Along the s ladder - a = 170 um fixed, so delta_omega = 340 um and theta
     = 75 um are frozen while delta* grows as s*a - does x_fc slope or stay
     flat? The aperture law predicts a 0.052 D_e drop from s=3 to s=7; the
     delta_omega and theta laws both predict zero. Slope versus flat line.
  2. Across the a ladder and the wall-tanh cases, which inlet scalar best
     collapses x_fc.
  3. Per case, is the averaging window long enough? Each case's own dumps give
     the cumulative-window drift, so convergence is measured rather than
     assumed.

x_fc uses E3, the midpoint crossing of the rising branch, which never
differentiates the level-composited centreline and so carries 2.2x less
sampling noise than the argmax-of-raw-gradient definition.
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

RHO_J = 1.6413

HERE = Path(__file__).resolve().parent
OUT = HERE / "sweep_out"

# name, a[um], s, form, arm
CASES = [
    ("A0_tophat",     0.0, 3.0, "shifted", "A"),
    ("A1_a100",     100.0, 3.0, "shifted", "A"),
    ("A2_a170",     170.0, 3.0, "shifted", "A/S"),
    ("A3_a254",     254.0, 3.0, "shifted", "A"),
    ("S4_a170_s4",  170.0, 4.0, "shifted", "S"),
    ("S5_a170_s5",  170.0, 5.0, "shifted", "S"),
    ("S6_a170_s6",  170.0, 6.0, "shifted", "S"),
    ("S7_a170_s7",  170.0, 7.0, "shifted", "S"),
    ("B1_wall680",  679.5, 3.0, "wall",    "B"),
    ("C1_wall200",  200.0, 3.0, "wall",    "C"),
]


def bracket(x, r):
    e = np.where((x >= 0.35) & (x <= 1.0))[0]
    i0 = e[np.nanargmin(r[e])]
    c = np.where((np.arange(x.size) > i0) & (x <= 1.7))[0]
    return i0, c[np.nanargmax(r[c])]


def xfc_E3(x, r):
    """Midpoint crossing of the rising branch: no differentiation at all."""
    i0, i1 = bracket(x, r)
    half = 0.5 * (r[i0] + r[i1])
    k = np.where(r[i0:i1 + 1] >= half)[0]
    if not len(k):
        return np.nan
    k = k[0] + i0
    if k == i0:
        return x[i0]
    return x[k-1] + (half - r[k-1]) * (x[k] - x[k-1]) / (r[k] - r[k-1])


def load(nm):
    f = OUT / ("XFC_%s.npz" % nm)
    if not f.exists():
        return None
    d = np.load(f)
    keys = sorted([k[:-2] for k in d.files if k.endswith("_t")],
                  key=lambda s: int("".join(c for c in s if c.isdigit())))
    if not keys:
        return None
    dx = float(d[keys[0] + "_t"][1])
    t = np.array([float(d[k + "_t"][0]) for k in keys])
    M = np.array([d[k + "_rho"] for k in keys]) / RHO_J
    x = (np.arange(M.shape[1]) + 0.5) * dx / D
    return x, t, M


rows = []
print("%-14s %6s %4s %-8s %9s %9s %9s   %s"
      % ("case", "a[um]", "s", "form", "d*_c[um]", "x_fc", "drift", "dumps"))
for nm, a, s, form, arm in CASES:
    got = load(nm)
    sc = scalars(a * 1e-6, form, s)
    if got is None:
        print("%-14s %6.1f %4.1f %-8s %9.1f %9s %9s   %s"
              % (nm, a, s, form, sc["dstar_c"] * 1e6, "-", "-", "MISSING"))
        continue
    x, t, M = got
    v = np.array([xfc_E3(x, m) for m in M])          # cumulative windows
    final = v[-1]
    drift = np.ptp(v[-3:]) if len(v) >= 3 else np.nan  # last three windows
    rows.append(dict(nm=nm, a=a, s=s, form=form, arm=arm, xfc=final,
                     drift=drift, n=len(v), sc=sc, vser=v, tser=t))
    print("%-14s %6.1f %4.1f %-8s %9.1f %9.4f %9.4f   %d"
          % (nm, a, s, form, sc["dstar_c"] * 1e6, final, drift, len(v)))

if not rows:
    sys.exit("\nno case data yet - run collect_sweep.sh on the head node first")

# ------------------------------------------------- the s-ladder verdict
lad = sorted([r for r in rows if r["a"] == 170.0], key=lambda r: r["s"])
print("\nTHE DISCRIMINATOR - s ladder at a = 170 um")
print("  delta_omega and theta are frozen; only the aperture moves.")
print("  %4s %10s %10s %10s" % ("s", "d*_c[um]", "x_fc", "x_fc-x_fc(s=3)"))
if lad:
    base = lad[0]["xfc"]
    for r in lad:
        print("  %4.1f %10.1f %10.4f %+10.4f"
              % (r["s"], r["sc"]["dstar_c"] * 1e6, r["xfc"], r["xfc"] - base))
    span = abs(lad[-1]["xfc"] - lad[0]["xfc"])
    print("\n  measured span over the ladder : %.4f D_e" % span)
    print("  aperture law predicts         : 0.052 D_e")
    print("  delta_omega / theta laws      : 0.000 D_e")
    noise = np.nanmean([r["drift"] for r in lad])
    print("  per-case window drift         : %.4f D_e" % noise)
    if span > 5 * max(noise, 1e-4):
        print("  --> the ladder SLOPES: x_fc is set by the effective aperture.")
    elif span < 2 * max(noise, 1e-4):
        print("  --> the ladder is FLAT: the aperture does NOT set x_fc.")
    else:
        print("  --> inconclusive at this window length; extend T_END.")

# ------------------------------------------------- predictor collapse
print("\nPREDICTOR COLLAPSE over every case")
y = np.array([r["xfc"] for r in rows])
print("  %-22s %10s %10s %8s" % ("predictor", "slope", "intercept", "R2"))
best = None
for key, lab in [("dstar_c", "d*_c (aperture)"), ("Cd", "C_d"),
                 ("CJ", "C_J"), ("dw", "delta_omega"), ("theta_c", "theta_c")]:
    xv = np.array([r["sc"][key] for r in rows])
    xn = xv / D if key not in ("Cd", "CJ") else xv
    if np.ptp(xn) == 0:
        continue
    m, c = np.polyfit(xn, y, 1)
    r2 = 1 - np.sum((y - (m * xn + c))**2) / np.sum((y - y.mean())**2)
    print("  %-22s %+10.4f %10.4f %8.4f" % (lab, m, c, r2))
    if best is None or r2 > best[1]:
        best = (lab, r2)
print("  best: %s (R2=%.4f)" % best)

# ------------------------------------------------------------- figure
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.4, 3.3))
if lad:
    a1.plot([r["s"] for r in lad], [r["xfc"] for r in lad], "o-",
            color="#b3251e", ms=6, lw=1.4, mec="w", mew=0.7, label="measured")
    ap = np.array([0.9443 - 1.9461 * r["sc"]["dstar_c"] / D for r in lad])
    a1.plot([r["s"] for r in lad], ap, "s--", color="#7a7f85", ms=5, lw=1.1,
            label="aperture law (predicted)")
    a1.axhline(lad[0]["xfc"], color="#1f4e9c", lw=1.1, ls=":",
               label=r"$\delta_\omega$ / $\theta$ laws (flat)")
a1.set_xlabel("shift constant $s$   ($a$=170 $\\mu$m fixed)")
a1.set_ylabel("$x_{fc}/D_e$")
a1.grid(True, ls="--", lw=0.5, color="0.88"); a1.set_axisbelow(True)
a1.legend(fontsize=6.8, frameon=False, loc="best", handlelength=1.6)
a1.text(0.03, 0.05, "(a) the discriminator", transform=a1.transAxes,
        fontsize=7.5)

mk = {"A": "o", "A/S": "o", "S": "^", "B": "s", "C": "D"}
cl = {"A": "#b3251e", "A/S": "#b3251e", "S": "#ff6a00", "B": "#1f4e9c",
      "C": "#2ca02c"}
for r in rows:
    a2.plot(r["sc"]["dstar_c"] / D, r["xfc"], mk.get(r["arm"], "o"),
            color=cl.get(r["arm"], "k"), ms=6, mec="w", mew=0.7)
xs = np.array([r["sc"]["dstar_c"] / D for r in rows])
if np.ptp(xs) > 0:
    m, c = np.polyfit(xs, y, 1)
    xx = np.linspace(xs.min(), xs.max(), 30)
    a2.plot(xx, m * xx + c, color="0.45", lw=1.0)
a2.set_xlabel(r"$\delta^*_c/D_e$   (compressible displacement thickness)")
a2.set_ylabel("$x_{fc}/D_e$")
a2.grid(True, ls="--", lw=0.5, color="0.88"); a2.set_axisbelow(True)
a2.legend(handles=[Line2D([], [], ls="", marker=mk[k], color=cl[k], ms=6,
                          label={"A": "A  a ladder", "A/S": "A2 anchor",
                                 "S": "S  s ladder", "B": "B  wall, $\\delta^*$-matched",
                                 "C": "C  wall a=200"}[k])
                   for k in ["A", "S", "B", "C"]],
           fontsize=6.6, frameon=False, loc="best", handlelength=1.2)
a2.text(0.03, 0.05, "(b) collapse", transform=a2.transAxes, fontsize=7.5)
fig.tight_layout(w_pad=1.3)
fig.savefig(HERE / "sweep_result.png", dpi=260)
fig.savefig(HERE / "sweep_result.pdf")
print("\nwrote sweep_result.png/pdf")
