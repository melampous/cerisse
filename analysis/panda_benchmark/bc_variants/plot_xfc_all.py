#!/usr/bin/env python3
"""Every x_fc measurement of the campaign on one axis, with the case list.

Left  - x_fc against the prescribed compressible displacement thickness, one
        colour per grid configuration. Black stars are the three-dimensional
        reference cases; the black open diamonds are the two-dimensional runs
        on a grid that reproduces the 3D refinement hierarchy level by level.
Right - the full case list: every inlet profile, its integral scalars, and its
        x_fc on each grid it has been run on.

Grids (all 2D axisymmetric unless stated):

  L3          198.4 um throughout the first shock cell. This is the resolution
              the 3D cases actually carry where x_fc is measured, their finest
              level being a lip annulus that stops at x = 0.75 D_e.
  L4          lip 99.2 um, whole core at 99.2 um.
  L5          lip 49.6 um, refinement forced on the lip band and on the axis.
  L5 lip-only lip 49.6 um, the same lip band, axis strip removed.
  mirror      domain, base mesh and all five level footprints copied from the
              3D case m142_L5_walltanh; statistics over 6.5 - 10.4 ms as in 3D.

x_fc is the thesis extraction throughout: the maximum raw density gradient
between the first expansion minimum and the first downstream principal
maximum, refined by a local quadratic fit, on the time-averaged centreline
density of the finest available level.
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
RHO_J = 1.6413
C3, M3 = 0.9443, -0.0766/1000.0
EXP, EXP_LO, EXP_HI = 0.888, 0.848, 0.928

# name, a[um], s, form, short label
PROF = [("A0_tophat",     0.0, 3.0, "shifted"), ("Shift_s0",   170.0, 0.0, "shifted"),
        ("Wall100",     100.0, 3.0, "wall"),    ("C1_wall200", 200.0, 3.0, "wall"),
        ("Shift_s1",    170.0, 1.0, "shifted"), ("Wall340",    340.0, 3.0, "wall"),
        ("A1_a100",     100.0, 3.0, "shifted"), ("Shift_s2",   170.0, 2.0, "shifted"),
        ("Wall470",     470.0, 3.0, "wall"),    ("B1_wall680", 679.5, 3.0, "wall"),
        ("A2_a170",     170.0, 3.0, "shifted"), ("S4_a170_s4", 170.0, 4.0, "shifted"),
        ("Wall900",     900.0, 3.0, "wall"),    ("A3_a254",    254.0, 3.0, "shifted"),
        ("S5_a170_s5",  170.0, 5.0, "shifted"), ("Wall1150",  1150.0, 3.0, "wall"),
        ("S6_a170_s6",  170.0, 6.0, "shifted"), ("S7_a170_s7", 170.0, 7.0, "shifted")]

# key, legend label, colour, marker, T0 [s], filename patterns
GRIDS = [("L3",   "L3   198.4 $\\mu$m everywhere",            "#CC79A7", "v", 4.0e-4,
          ("XFC_%s_LO__L3.npz", "XFC_%s__L3.npz")),
         ("L4",   "L4   lip 99.2 $\\mu$m",                    "#0072B2", "o", 4.0e-4,
          ("XFC_%s__L4.npz", "XFC_%s_L4.npz")),
         ("L5",   "L5   lip 49.6 $\\mu$m",                    "#D55E00", "s", 4.0e-4,
          ("XFC_%s__L5.npz", "XFC_%s.npz")),
         ("L5lo", "L5   lip 49.6 $\\mu$m, no axis strip",     "#009E73", "^", 4.0e-4,
          ("XFC_%s_LO__L5.npz",)),
         ("MIR",  "2D on the mirrored 3D grid",               "#111111", "D", 6.5e-3,
          ("XFC_mirror_%s__MIR.npz",))]
MIRMAP = {"A0_tophat": "tophat", "C1_wall200": "wall200"}

THREE_D = {"A0_tophat": 0.945, "C1_wall200": 0.930,
           "A1_a100": 0.922, "A3_a254": 0.883}
LBL3D = {"A0_tophat": "top hat", "C1_wall200": "wall $a$=200",
         "A1_a100": "shifted $a$=100", "A3_a254": "shifted $a$=254"}


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


def measure(pats, nm, T0):
    for pat in pats:
        f = OUT/(pat % (MIRMAP.get(nm, nm) if "mirror" in pat else nm))
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


rows = []
for nm, a, s, form in PROF:
    sc = scalars(a*1e-6, form, s)
    r = dict(nm=nm, ds=sc["dstar_c"]*1e6, dw=sc["dw"]*1e6, H=sc["H"],
             three=THREE_D.get(nm))
    for key, lab, col, mk, T0, pats in GRIDS:
        if key == "MIR" and nm not in MIRMAP:
            r[key] = (None, None)
            continue
        r[key] = measure(pats, nm, T0)
    rows.append(r)
rows.sort(key=lambda q: q["ds"])

# --------------------------------------------------------------- console
hdr = "%-12s %7s %6s | %6s" % ("case", "d*[um]", "H", "3D")
for key, *_ in GRIDS:
    hdr += " %8s" % key
print(hdr)
for r in rows:
    line = "%-12s %7.1f %6s | %6s" % (
        r["nm"], r["ds"], "-" if not np.isfinite(r["H"]) else "%.2f" % r["H"],
        "%.3f" % r["three"] if r["three"] else "-")
    for key, *_ in GRIDS:
        v = r[key][0]
        line += " %8s" % ("%.4f" % v if v else "-")
    print(line)
for key, lab, *_ in GRIDS:
    P = [(r["ds"], r[key][0]) for r in rows if r[key][0]]
    if len(P) < 3:
        continue
    X = np.array([p[0] for p in P]); Y = np.array([p[1] for p in P])
    m, c = np.polyfit(X, Y, 1)
    res = Y-(C3+M3*X)
    print("%-34s n=%2d  slope %+.4f/1000um  rms from the 3D line %.4f"
          % (key, len(P), m*1000, np.sqrt((res**2).mean())))

# ---------------------------------------------------------------- figure
plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})
fig = plt.figure(figsize=(13.6, 6.3))
gs = fig.add_gridspec(1, 2, width_ratios=[1.52, 1.0], wspace=0.06,
                      left=0.055, right=0.995, top=0.965, bottom=0.095)
ax = fig.add_subplot(gs[0, 0])
axt = fig.add_subplot(gs[0, 1]); axt.axis("off")

ax.axhline(EXP, color="0.55", lw=1.0, zorder=1)
for y in (EXP_LO, EXP_HI):
    ax.axhline(y, color="0.82", lw=0.6, ls=":", zorder=1)
ax.text(1292, EXP+0.0020, "experiment  %.3f" % EXP, ha="right", fontsize=7.2,
        color="0.35")

xx = np.linspace(-40, 1300, 40)
ax.plot(xx, C3+M3*xx, color="0.25", lw=1.2, ls="--", zorder=2)
ax.text(1120, C3+M3*1120+0.0078, "fit through the 3D cases", fontsize=7,
        rotation=-8.5, color="0.25", ha="center")

for key, lab, col, mk, T0, _ in GRIDS:
    P = [r for r in rows if r[key][0] and r["nm"] != "A0_tophat"]
    T = [r for r in rows if r[key][0] and r["nm"] == "A0_tophat"]
    if P and key != "MIR":
        ax.plot([r["ds"] for r in P], [r[key][0] for r in P], color=col,
                lw=1.0, alpha=0.45, zorder=3)
    for r in P:
        ax.errorbar(r["ds"], r[key][0], yerr=r[key][1], fmt=mk, color=col,
                    ms=7.0 if key == "MIR" else 5.6,
                    mfc="w" if key == "MIR" else col, mec=col,
                    mew=1.6 if key == "MIR" else 0.7,
                    capsize=2.4, lw=1.0, ls="", zorder=6 if key == "MIR" else 5)
    for r in T:                                    # top hat: bigger, open
        ax.errorbar(r["ds"], r[key][0], yerr=r[key][1], fmt=mk, color=col,
                    ms=11.0, mfc="w", mec=col, mew=2.0, capsize=3.0, lw=1.2,
                    zorder=7)

th = sorted(r[k][0] for r in rows if r["nm"] == "A0_tophat"
            for k, *_ in GRIDS if r[k][0])
if len(th) >= 2:
    ax.annotate("", xy=(0, th[0]), xytext=(0, th[-1]),
                arrowprops=dict(arrowstyle="<->", lw=0.9, color="0.45",
                                shrinkA=10, shrinkB=10), zorder=4)
    ax.annotate("top hat, $\\delta^*$ = 0: the discontinuous limit.\n"
                "No shear layer is prescribed, so the lip gradient is set\n"
                "by the mesh and the four grids span %.3f $D_e$."
                % (th[-1]-th[0]),
                xy=(215, 0.8555), xytext=(215, 0.8555), textcoords="data",
                fontsize=7.4, color="0.2", va="bottom")

for nm, v in THREE_D.items():
    r = [q for q in rows if q["nm"] == nm][0]
    ax.plot(r["ds"], v, "*", color="k", ms=14, zorder=8)
    ax.annotate(LBL3D[nm], (r["ds"], v), textcoords="offset points",
                xytext=(9, 4), fontsize=7)

ax.set_xlabel(r"prescribed displacement thickness  $\delta^*_c$  [$\mu$m]")
ax.set_ylabel("$x_{fc}/D_e$")
ax.set_xlim(-55, 1300)
ax.set_ylim(0.845, 0.952)
ax.grid(True, ls="--", lw=0.5, color="0.93")
ax.set_axisbelow(True)

H = [Line2D([], [], ls="", marker="*", color="k", ms=13,
            label="3D reference (4 cases)")]
for key, lab, col, mk, T0, _ in GRIDS:
    n = sum(1 for r in rows if r[key][0])
    H.append(Line2D([], [], color=col, marker=mk, ms=6.4 if key == "MIR" else 5.6,
                    mfc="w" if key == "MIR" else col, mec=col,
                    mew=1.5 if key == "MIR" else 0.7,
                    ls="" if key == "MIR" else "-",
                    label="2D  %s   (%d)" % (lab, n)))
H += [Line2D([], [], ls="", marker="o", color="0.35", mfc="w", ms=9, mew=2.0,
             label="large open symbol: top hat, $\\delta^*$ = 0"),
      Line2D([], [], color="0.25", lw=1.2, ls="--",
             label="correlation fitted to the 3D cases")]
ax.legend(handles=H, fontsize=7.1, frameon=False, loc="upper right",
          handlelength=1.9, labelspacing=0.36, borderaxespad=0.7)

# ------------------------------------------------------------ case table
def cell(v):
    return "%.4f" % v if v else "  -   "


lines = ["%-12s %6s %5s %6s %6s %6s %6s %6s %6s"
         % ("case", "d*_c", "H", "3D", "L3", "L4", "L5", "L5lo", "MIR"),
         "%-12s %6s %5s %6s %6s %6s %6s %6s %6s"
         % ("", "[um]", "", "", "198.4", "99.2", "49.6", "49.6", "3D")]
for r in rows:
    lines.append("%-12s %6.1f %5s %6s %6s %6s %6s %6s %6s"
                 % (r["nm"], r["ds"],
                    "-" if not np.isfinite(r["H"]) else "%.2f" % r["H"],
                    "%.3f " % r["three"] if r["three"] else "  -   ",
                    cell(r["L3"][0]), cell(r["L4"][0]), cell(r["L5"][0]),
                    cell(r["L5lo"][0]), cell(r["MIR"][0])))

axt.text(0.0, 1.0, "\n".join(lines), family="monospace", fontsize=7.0,
         va="top", ha="left", linespacing=1.62, transform=axt.transAxes)
axt.plot([0.0, 0.985], [0.955, 0.955], color="0.4", lw=0.8,
         transform=axt.transAxes, clip_on=False)
axt.plot([0.0, 0.985], [0.916, 0.916], color="0.75", lw=0.6,
         transform=axt.transAxes, clip_on=False)

note = (
    "Inlet forms.  shifted tanh  f = 1/2 [1 - tanh((r-r0)/a)],  r0 = R - s a\n"
    "              wall tanh     f = tanh(max(R-r,0)/a)\n"
    "              top hat       a = 0\n"
    "H = delta*/theta is the shape factor: 2.26 for every wall tanh,\n"
    "6.0 to 14.0 along the shifted family, 1.4 to 4.1 for s = 0, 1, 2.\n"
    "\n"
    "Sampling error 0.0012 - 0.0124 D_e (error bars on the left).\n"
    "Statistics over 0.4 - 2.4 ms, except the mirror runs, which use the\n"
    "3D window 6.5 - 10.4 ms.\n"
    "\n"
    "Caveat on the 3D column: the four reference cases are not on one grid.\n"
    "wall a=200 and shifted a=100 reach 99.2 um at the lip only, and x_fc\n"
    "sits downstream of that patch, on 198.4 um; shifted a=254 is 198.4 um\n"
    "throughout and averaged over a different interval (5.06 ms). Their\n"
    "x_fc is quantised at 198.4 um = 0.0078 D_e by the extraction."
)
axt.text(0.0, 0.017, note, family="monospace", fontsize=6.4, va="bottom",
         ha="left", color="0.3", linespacing=1.5, transform=axt.transAxes)

fig.savefig(HERE/"xfc_all.png", dpi=250)
fig.savefig(HERE/"xfc_all.pdf")
print("\nwrote xfc_all.png/pdf")
