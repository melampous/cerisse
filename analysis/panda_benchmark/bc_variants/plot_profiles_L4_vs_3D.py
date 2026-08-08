#!/usr/bin/env python3
"""The inlet-profile sensitivity of x_fc: the four 3D reference cases and the
complete 2D sweep on the L4 grid, every point labelled with its profile.

Only two data sets appear here. The 2D points are all on one grid (lip cells
99.2 um), so differences between them are differences between inlet profiles
and nothing else. The 3D points are shown for reference together with the line
fitted through them; note that their x_fc is measured on 198.4 um cells,
because their finest level is a lip annulus that stops upstream of x_fc.

Inlet forms, all truncated at r = R with an adiabatic slip flange beyond:
  top hat        f = 1 for r <= R
  wall tanh      f = tanh(max(R - r, 0)/a)                      H = 2.26
  shifted tanh   f = 1/2 [1 - tanh((r - r0)/a)],  r0 = R - s a   H = 4.1 - 14.0

x_fc is the thesis extraction: maximum raw density gradient between the first
expansion minimum and the first downstream principal maximum, refined by a
local quadratic fit, on the time-averaged centreline density of the finest
available level. Error bars are the standard error of the disjoint sub-windows
of the 2.0 ms averaging interval.
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
EXP, EXP_LO, EXP_HI = 0.888, 0.848, 0.928

CW, CS, CT, CH = "#b3251e", "#1f4e9c", "#6a3d9a", "#111111"

# case, a[um], s, form, label, family, label offset (pts)
PROF = [
    ("A0_tophat",     0.0, 3.0, "shifted", "top hat",            "hat",   (11,  -4)),
    ("Shift_s0",    170.0, 0.0, "shifted", "shifted 170, s=0",   "trunc", (10,   9)),
    ("Wall100",     100.0, 3.0, "wall",    "wall 100",           "wall",  (8,  -12)),
    ("C1_wall200",  200.0, 3.0, "wall",    "wall 200",           "wall",  (-20, -14)),
    ("Shift_s1",    170.0, 1.0, "shifted", "shifted 170, s=1",   "trunc", (-96, -6)),
    ("Wall340",     340.0, 3.0, "wall",    "wall 340",           "wall",  (-18, -14)),
    ("A1_a100",     100.0, 3.0, "shifted", "shifted 100, s=3",   "shift", (-112, -3)),
    ("Shift_s2",    170.0, 2.0, "shifted", "shifted 170, s=2",   "trunc", (7,    7)),
    ("Wall470",     470.0, 3.0, "wall",    "wall 470",           "wall",  (8,  -11)),
    ("B1_wall680",  679.5, 3.0, "wall",    "wall 680",           "wall",  (-70,-12)),
    ("A2_a170",     170.0, 3.0, "shifted", "shifted 170, s=3",   "shift", (8,    5)),
    ("S4_a170_s4",  170.0, 4.0, "shifted", "shifted 170, s=4",   "shift", (-100, 5)),
    ("Wall900",     900.0, 3.0, "wall",    "wall 900",           "wall",  (7,  -12)),
    ("A3_a254",     254.0, 3.0, "shifted", "shifted 254, s=3",   "shift", (11, -13)),
    ("S5_a170_s5",  170.0, 5.0, "shifted", "shifted 170, s=5",   "shift", (-102, 3)),
    ("Wall1150",   1150.0, 3.0, "wall",    "wall 1150",          "wall",  (7,  -12)),
    ("S6_a170_s6",  170.0, 6.0, "shifted", "shifted 170, s=6",   "shift", (8,    4)),
    ("S7_a170_s7",  170.0, 7.0, "shifted", "shifted 170, s=7",   "shift", (-100, 4)),
]
COL = {"wall": CW, "shift": CS, "trunc": CT, "hat": CH}
MRK = {"wall": "s",  "shift": "o", "trunc": "^", "hat": "P"}

THREE_D = [("A0_tophat",  0.0, 3.0, "shifted", 0.945, "top hat",          (10, 6)),
           ("C1_wall200", 200.0, 3.0, "wall",  0.930, "wall 200",         (10, 5)),
           ("A1_a100",  100.0, 3.0, "shifted", 0.922, "shifted 100, s=3", (10, 4)),
           ("A3_a254",  254.0, 3.0, "shifted", 0.883, "shifted 254, s=3", (10, 4))]


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


def measure(nm, grid="L4"):
    pats = (("XFC_%s__L4.npz", "XFC_%s_L4.npz") if grid == "L4"
            else ("XFC_%s__L5.npz", "XFC_%s.npz"))
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
    # The reported value is x_fc of the full-window cumulative mean. Its
    # uncertainty is taken as how far that same quantity still moves as the
    # window lengthens - the largest departure of the running value from its
    # final value over the second half of the window - floored at half a cell,
    # which is the resolution of the extraction itself. This avoids the two
    # assumptions behind std/sqrt(n) over disjoint sub-windows: that the
    # sub-windows are independent, and that x_fc commutes with time averaging
    # (measured departure from commutation is 0.0047 D_e on average, as large
    # as the sub-window standard error itself).
    cum = np.array([xfc(x, m) for m in M])
    full = cum[-1]
    drift = np.abs(cum[len(cum)//2:] - full).max() if len(cum) > 2 else np.nan
    half_cell = 0.5*dx/D
    return full, max(drift, half_cell)


rows = []
for nm, a, s, form, lab, fam, off in PROF:
    sc = scalars(a*1e-6, form, s)
    v, e = measure(nm, "L4")
    v5, e5 = measure(nm, "L5")
    if v is None:
        continue
    rows.append(dict(nm=nm, lab=lab, fam=fam, off=off, a=a, s=s, form=form,
                     v5=v5, e5=e5,
                     ds=sc["dstar_c"]*1e6, dw=sc["dw"]*1e6, th=sc["theta_c"]*1e6,
                     H=sc["H"], Cd=sc["Cd"], v=v, e=e))
rows.sort(key=lambda q: q["ds"])

X3 = np.array([scalars(a*1e-6, f, s)["dstar_c"]*1e6 for _, a, s, f, *_ in THREE_D])
Y3 = np.array([q[4] for q in THREE_D])
m3, c3 = np.polyfit(X3, Y3, 1)
r2 = 1-((Y3-(m3*X3+c3))**2).sum()/((Y3-Y3.mean())**2).sum()

print("2D on the L4 grid, lip cells 99.2 um")
print("%-20s %8s %8s %8s %6s %7s | %8s %8s"
      % ("profile", "a[um]", "d*_c", "theta_c", "H", "C_d", "x_fc", "sigma"))
for r in rows:
    print("%-20s %8.0f %8.1f %8.1f %6s %7.4f | %8.4f %8.4f"
          % (r["lab"], r["a"], r["ds"], r["th"],
             "-" if not np.isfinite(r["H"]) else "%.2f" % r["H"],
             r["Cd"], r["v"], r["e"]))
print("\n3D reference, x_fc measured on 198.4 um cells")
for nm, a, s, f, v, lab, _ in THREE_D:
    print("%-20s %8.0f %8.1f %31.3f" % (lab, a, scalars(a*1e-6, f, s)["dstar_c"]*1e6, v))
print("\nline through the 3D points: x_fc = %.4f %+.4f d*/1000um   R2 = %.4f"
      % (c3, m3*1000, r2))
res = np.array([r["v"]-(c3+m3*r["ds"]) for r in rows])
print("2D L4 about that line: rms %.4f  max %.4f" % (np.sqrt((res**2).mean()),
                                                     np.abs(res).max()))

# ---------------------------------------------------------------- figure
plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})
fig = plt.figure(figsize=(12.8, 6.4))
gs = fig.add_gridspec(1, 2, width_ratios=[1.72, 1.0], wspace=0.05,
                      left=0.055, right=0.995, top=0.965, bottom=0.093)
ax = fig.add_subplot(gs[0, 0])
axt = fig.add_subplot(gs[0, 1]); axt.axis("off")

ax.axhline(EXP, color="0.55", lw=1.0, zorder=1)
for y in (EXP_LO, EXP_HI):
    ax.axhline(y, color="0.84", lw=0.6, ls=":", zorder=1)
ax.text(1295, EXP+0.0020, "experiment  %.3f" % EXP, ha="right", fontsize=7.2,
        color="0.35")

xx = np.linspace(-40, 1300, 40)
ax.plot(xx, c3+m3*xx, color="0.3", lw=1.2, ls="--", zorder=2)
ax.text(1090, c3+m3*1090+0.0080,
        "line through the 3D points\n$x_{fc}$ = %.4f $-$ %.4f $\\delta^*_c$/1000$\\,\\mu$m"
        % (c3, -m3*1000), fontsize=7, rotation=-8.5, color="0.3", ha="center")

for r in rows:                                   # L5 twin, and the grid shift
    if r["v5"] is not None:
        ax.plot([r["ds"], r["ds"]], [r["v"], r["v5"]], color="0.62", lw=0.9,
                zorder=3)
        ax.errorbar(r["ds"], r["v5"], yerr=r["e5"], fmt=MRK[r["fam"]],
                    color=COL[r["fam"]], ms=6.6 if r["fam"] == "hat" else 5.4,
                    mfc="w", mec=COL[r["fam"]], mew=1.4, capsize=2.2, lw=0.9,
                    ls="", zorder=4)
for r in rows:
    ax.errorbar(r["ds"], r["v"], yerr=r["e"], fmt=MRK[r["fam"]],
                color=COL[r["fam"]], ms=8.0 if r["fam"] == "hat" else 6.2,
                mec="w", mew=0.8, capsize=2.6, lw=1.0, ls="", zorder=5)
    ax.annotate(r["lab"], (r["ds"], r["v"]), textcoords="offset points",
                xytext=r["off"], fontsize=6.9, color=COL[r["fam"]])

for nm, a, s, f, v, lab, off in THREE_D:
    x = scalars(a*1e-6, f, s)["dstar_c"]*1e6
    ax.plot(x, v, "*", color="k", ms=15, zorder=8)
    ax.annotate(lab, (x, v), textcoords="offset points", xytext=off,
                fontsize=7.4, color="k")

ax.set_xlabel(r"prescribed displacement thickness  $\delta^*_c$  [$\mu$m]")
ax.set_ylabel("$x_{fc}/D_e$")
ax.set_xlim(-60, 1305)
ax.set_ylim(0.846, 0.952)
ax.grid(True, ls="--", lw=0.5, color="0.93")
ax.set_axisbelow(True)
ax.legend(handles=[
    Line2D([], [], ls="", marker="*", color="k", ms=13,
           label="3D reference (4 cases, 198.4 $\\mu$m at $x_{fc}$)"),
    Line2D([], [], ls="", marker="P", color=CH, ms=8, mec="w",
           label="2D L4  top hat,  $\\delta^*$ = 0"),
    Line2D([], [], ls="", marker="s", color=CW, ms=6.2, mec="w",
           label="2D L4  wall tanh,  $H$ = 2.26  (7)"),
    Line2D([], [], ls="", marker="o", color=CS, ms=6.2, mec="w",
           label="2D L4  shifted tanh,  $H$ = 6.0 - 14.0  (7)"),
    Line2D([], [], ls="", marker="^", color=CT, ms=6.2, mec="w",
           label="2D L4  truncation series  $s$ = 0, 1, 2  (3)"),
    Line2D([], [], ls="", marker="o", color="0.35", mfc="w", ms=5.6, mew=1.4,
           label="open symbol: same profile on 2D L5 (lip 49.6 $\\mu$m, 10)"),
    Line2D([], [], color="0.3", lw=1.2, ls="--", label="fit through the 3D points")],
    fontsize=7.2, frameon=False, loc="lower left", handlelength=1.7,
    labelspacing=0.38, borderaxespad=0.7)

# ------------------------------------------------------------ side table
L = ["%-20s %6s %5s %6s %8s %8s" % ("profile", "d*_c", "H", "C_d", "x_fc L4", "x_fc L5"),
     "%-20s %6s %5s %6s %8s %8s" % ("", "[um]", "", "", "99.2um", "49.6um")]
for r in rows:
    L.append("%-20s %6.1f %5s %6.4f %8.4f %8s" % (
        r["lab"], r["ds"],
        "-" if not np.isfinite(r["H"]) else "%.2f" % r["H"], r["Cd"], r["v"],
        "%.4f" % r["v5"] if r["v5"] is not None else "   -"))
L.append("")
L.append("3D reference (x_fc on 198.4 um cells)")
for nm, a, s, f, v, lab, _ in THREE_D:
    L.append("%-20s %6.1f %5s %6.4f %8.3f %8s" % (
        lab, scalars(a*1e-6, f, s)["dstar_c"]*1e6,
        "-" if a == 0 else "%.2f" % scalars(a*1e-6, f, s)["H"],
        scalars(a*1e-6, f, s)["Cd"], v, "   -"))
axt.text(0.0, 1.0, "\n".join(L), family="monospace", fontsize=7.0, va="top",
         ha="left", linespacing=1.60, transform=axt.transAxes)
axt.plot([0.0, 0.99], [0.958, 0.958], color="0.4", lw=0.8,
         transform=axt.transAxes, clip_on=False)

note = ("a is the tanh width in um; s the shift of the shifted form, r0 = R - s a.\n"
        "delta*_c is the compressible displacement thickness of the prescribed\n"
        "mass-flux profile; H = delta*/theta; C_d = A_eff/A.\n"
        "\n"
        "Error bar: how far x_fc of the running cumulative mean still moves over\n"
        "the second half of the averaging window, floored at half a cell (the\n"
        "resolution of the extraction: 0.0020 D_e on L4, 0.0010 on L5). This is\n"
        "defined on the reported quantity itself, unlike std/sqrt(n) over\n"
        "disjoint sub-windows, which assumes independence and that x_fc commutes\n"
        "with averaging - the measured departure from commutation is 0.0047 D_e.\n"
        "\n"
        "All 2D L4 points share one grid, one window (0.4 - 2.4 ms) and one\n"
        "extraction, so the spread between them is the inlet profile alone.")
axt.text(0.0, 0.020, note, family="monospace", fontsize=6.5, va="bottom",
         ha="left", color="0.32", linespacing=1.5, transform=axt.transAxes)

fig.savefig(HERE/"profiles_L4_vs_3D.png", dpi=250)
fig.savefig(HERE/"profiles_L4_vs_3D.pdf")
print("\nwrote profiles_L4_vs_3D.png/pdf")
