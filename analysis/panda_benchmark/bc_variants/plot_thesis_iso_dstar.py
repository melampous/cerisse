#!/usr/bin/env python3
"""Two thesis-format figures for the iso-delta* family.

Figure A follows Figure 3.18 (profile_controlled_families): three panels
showing the full radial profile, the lip region and the radial shear, with the
new family drawn against the two families already in the chapter.

Figure B follows Figure 3.19 (profile_xfc_correlations): x_fc against each
candidate exit-profile parameter, black stars for the three-dimensional
Cartesian reference and coloured symbols for the axisymmetric sweep, with the
linear fit to the three-dimensional values.

Notation follows the chapter: delta_s is the tanh width parameter, s the
shift in units of delta_s, delta*_c and theta_c the compressible displacement
and momentum thicknesses.

Two departures from the chapter figures, both forced by the data:

  wall family  only two wall-attached profiles were run at Level 5, so that
          family is thinner here than in the chapter, where the Level 4 sweep
          carried six. The remaining wall cases exist at Level 4 only.

  level   the chapter plots the axisymmetric sweep at Level 4. The new family
          was run at Level 5, and the two levels do not share an offset: for
          the shifted profile at delta_s = 170 um the difference is +0.0202
          D_e, for the wall profile at 679.5 um it is +0.0021 D_e. Mixing them
          would put a grid effect of the same size as the whole spread of the
          family into the comparison, so every axisymmetric point here is
          Level 5, taken from the same _LO__L5 extraction.

  abscissa  x_fc is computed on the abscissa the extractor writes, dx stored
          in the _t record, which spans 8.2 D_e. A hard-coded 24 D_e would
          scale every value by 2.93.
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
OUT = HERE/"sweep_out"
RHO_J = 1.6413
CW, CS, CD, CT = "#b3251e", "#1f4e9c", "#0b7d63", "#e8a33d"

# ------------------------------------------------------------------ measured
def refine(x, y, i):
    if i <= 0 or i >= len(y)-1:
        return x[i]
    d = y[i-1]-2*y[i]+y[i+1]
    return x[i] if d == 0 else x[i]+0.5*(y[i-1]-y[i+1])/d*(x[1]-x[0])


def xfc_of(f):
    d = np.load(OUT/f)
    ks = sorted([k[:-2] for k in d.files if k.endswith("_t")],
                key=lambda q: int("".join(c for c in q if c.isdigit())))
    dx = float(d[ks[0]+"_t"][1])
    r = d[ks[-1]+"_rho"]/RHO_J
    x = (np.arange(r.size)+0.5)*dx/D
    m = (x >= 0.35) & (x <= 1.05) & np.isfinite(r)
    xx, rr = x[m], r[m]
    i0 = int(np.argmin(rr))
    i1 = i0 + int(np.argmax(rr[np.arange(i0, len(rr))]))
    g = np.gradient(rr, xx)
    s = np.arange(i0, i1+1)
    return refine(xx, g, s[int(np.nanargmax(g[s]))])


#            label            delta_s   s     form      family   file
SWEEP = [("top hat",              0.0, 3.0, "shifted", "hat",  "A0_tophat_LO__L5"),
         ("wall 200",           200.0, 3.0, "wall",    "wall", "C1_wall200_LO__L5"),

         ("shift 100",          100.0, 3.0, "shifted", "shift", "A1_a100_LO__L5"),
         ("shift 254",          254.0, 3.0, "shifted", "shift", "A3_a254_LO__L5"),
         ("shift 170 s4",       170.0, 4.0, "shifted", "shift", "S4_a170_s4_LO__L5"),
         ("shift 170 s5",       170.0, 5.0, "shifted", "shift", "S5_a170_s5_LO__L5"),
         ("shift 170 s6",       170.0, 6.0, "shifted", "shift", "S6_a170_s6_LO__L5"),
         ("shift 170 s7",       170.0, 7.0, "shifted", "shift", "S7_a170_s7_LO__L5"),
         ("iso 100",            100.0, 5.195590, "shifted", "iso", "D1_a100_s5p196_LO__L5"),
         ("iso 130",            130.0, 3.965528, "shifted", "iso", "D2_a130_s3p966_LO__L5"),
         ("iso 170",            170.0, 3.000000, "shifted", "iso", "A2_a170_LO__L5"),
         ("iso 220",            220.0, 2.284208, "shifted", "iso", "D4_a220_s2p284_LO__L5"),
         ("iso 679.5 wall",     679.5, 3.0, "wall",    "iso",  "B1_wall680_LO__L5")]

# the four Cartesian references, x_fc from the chapter table
THREE = [("top hat",   0.0, 3.0, "shifted", 0.945),
         ("wall 200",  200.0, 3.0, "wall",   0.930),
         ("shift 100", 100.0, 3.0, "shifted", 0.922),
         ("shift 254", 254.0, 3.0, "shifted", 0.883)]

ROWS = []
for lab, ds, s, form, fam, f in SWEEP:
    sc = scalars(ds*1e-6, form, s)
    ROWS.append(dict(lab=lab, ds=ds, s=s, form=form, fam=fam,
                     dsc=sc["dstar_c"]*1e6, thc=sc["theta_c"]*1e6,
                     deff=sc["Deff"]/D, xfc=xfc_of("XFC_%s.npz" % f)))

print("%-16s %6s %8s %10s %10s %8s %9s"
      % ("case", "fam", "delta_s", "d*_c[um]", "th_c[um]", "Deff/D", "x_fc"))
for r in ROWS:
    print("%-16s %6s %8.1f %10.1f %10.1f %8.4f %9.4f"
          % (r["lab"], r["fam"], r["ds"], r["dsc"], r["thc"], r["deff"], r["xfc"]))

iso = [r for r in ROWS if r["fam"] == "iso"]
v = np.array([r["xfc"] for r in iso])
print("\niso-delta* set: n=%d  d*_c=%.1f um  x_fc %.4f +/- %.4f  spread %.4f D_e"
      % (len(iso), iso[0]["dsc"], v.mean(), v.std(ddof=1), np.ptp(v)))
print("   delta_omega spans %.0f to %.0f um, H spans %.2f to %.2f"
      % (min(2*r["ds"] if r["form"] == "shifted" else r["ds"] for r in iso),
         max(2*r["ds"] if r["form"] == "shifted" else r["ds"] for r in iso),
         min(r["dsc"]/r["thc"] for r in iso), max(r["dsc"]/r["thc"] for r in iso)))

# 3D fit, as in the chapter
p3 = np.array([scalars(a*1e-6, f, s)["dstar_c"]*1e6 for _, a, s, f, _ in THREE])
y3 = np.array([q[4] for q in THREE])
m3, c3 = np.polyfit(p3/1000.0, y3, 1)
r2 = 1-((y3-(c3+m3*p3/1000.0))**2).sum()/((y3-y3.mean())**2).sum()
print("\n3D fit  x_fc/D_e = %.4f %+.4f (d*_c/1000um)   R2 = %.4f" % (c3, m3, r2))

# ------------------------------------------------------------------ figure A
plt.rcParams.update({"font.size": 8.5, "axes.linewidth": 0.6,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True,
                     "legend.frameon": False, "axes.titlesize": 9})
FAMA = [(100.0, 5.195590), (130.0, 3.965528), (170.0, 3.000000), (220.0, 2.284208)]
GR = plt.get_cmap("viridis")(np.linspace(0.12, 0.80, len(FAMA)))
r = np.linspace(0.0, R, 400001)
rmm, Rmm = r*1e3, R*1e3

fig, AX = plt.subplots(1, 3, figsize=(9.2, 3.2))
fig.subplots_adjust(left=0.070, right=0.988, top=0.905, bottom=0.145, wspace=0.30)
for (ds, s), col in zip(FAMA, GR):
    f = fvel(r, ds*1e-6, "shifted", s)
    g = np.abs(np.gradient(f, rmm))
    lab = r"$\delta_s$ = %d $\mu$m, $s$ = %.2f" % (ds, s)
    AX[0].plot(f, rmm, color=col, lw=1.4, zorder=5, label=lab)
    AX[1].plot(f, rmm, color=col, lw=1.4, zorder=5)
    AX[2].plot(rmm, g, color=col, lw=1.4, zorder=5)
fw = fvel(r, 679.5e-6, "wall", 3.0)
AX[0].plot(fw, rmm, color=CW, lw=1.4, ls="--", zorder=4,
           label=r"wall, $\delta_s$ = 679.5 $\mu$m")
AX[1].plot(fw, rmm, color=CW, lw=1.4, ls="--", zorder=4)
AX[2].plot(rmm, np.abs(np.gradient(fw, rmm)), color=CW, lw=1.4, ls="--", zorder=4)

for k, (ax, ttl) in enumerate(zip(AX, ["(a)  full radial profile", "(b)  lip region",
                                       "(c)  finite radial shear"])):
    ax.grid(True, ls="--", lw=0.4, color="0.92")
    ax.set_axisbelow(True)
    ax.set_title(ttl, pad=4)
    if k < 2:
        ax.axhspan(Rmm, 20.0, color="0.915", zorder=0)
        ax.axhline(Rmm, color="0.4", lw=0.7, zorder=2)
        ax.set_xlim(-0.02, 1.06)
        ax.set_xlabel("profile function  $f = u/U_e$")
        ax.set_ylabel("radius  $r$  [mm]")
        ax.set_ylim((0.0, 13.2) if k == 0 else (10.9, 13.0))
    else:
        ax.axvspan(Rmm, 20.0, color="0.915", zorder=0)
        ax.axvline(Rmm, color="0.4", lw=0.7, zorder=2)
        ax.set_xlim(10.9, 13.0)
        ax.set_ylim(0, 5.45)
        ax.set_xlabel("radius  $r$  [mm]")
        ax.set_ylabel(r"$|\mathrm{d}f/\mathrm{d}r|$   [1/mm]")
AX[0].legend(loc="center left", bbox_to_anchor=(0.05, 0.42), fontsize=7.6,
             handlelength=1.6, labelspacing=0.32,
             title=r"$\delta^*_c$ = 533 $\mu$m throughout")
AX[0].get_legend().get_title().set_fontsize(7.8)
fig.savefig(HERE/"thesis_iso_dstar_families.pdf")
fig.savefig(HERE/"thesis_iso_dstar_families.png", dpi=300)
plt.close(fig)
print("\nwrote thesis_iso_dstar_families.pdf/png")

# ------------------------------------------------------------------ figure B
# tag goes in whichever corner the data leaves free in that panel
PAR = [("dsc", r"$\delta^*_c$   [$\mu$m]", "(a)", 0.965, 0.93, "right"),
       ("thc", r"$\theta_c$   [$\mu$m]", "(b)", 0.965, 0.93, "right"),
       ("deff", r"$D_{\mathrm{eff},q}/D_e=\sqrt{C_{d,\mathrm{eff},q}}$", "(c)", 0.035, 0.93, "left")]
fig, AX = plt.subplots(1, 3, figsize=(9.2, 3.05))
fig.subplots_adjust(left=0.070, right=0.988, top=0.965, bottom=0.315, wspace=0.24)
STY = dict(hat=(CT, "P", 8.0), wall=(CW, "s", 5.2),
           shift=(CS, "o", 5.2), iso=(CD, "D", 5.6))
for ax, (key, xlab, tag, tx, ty, ha) in zip(AX, PAR):
    ax.axhline(0.888, color="0.55", lw=0.7, zorder=1)
    for r_ in ROWS:
        col, mk, ms = STY[r_["fam"]]
        ax.plot(r_[key], r_["xfc"], ls="", marker=mk, ms=ms, color=col,
                mec="w" if r_["fam"] == "iso" else "none", mew=0.6, zorder=5)
    for lab, a, s, form, y in THREE:
        sc = scalars(a*1e-6, form, s)
        xv = dict(dsc=sc["dstar_c"]*1e6, thc=sc["theta_c"]*1e6,
                  deff=sc["Deff"]/D)[key]
        ax.plot(xv, y, ls="", marker="P" if a == 0 else "*",
                ms=9 if a == 0 else 12, color="#111111", zorder=6)
    if key in ("dsc", "deff"):
        xs = np.linspace(*ax.get_xlim(), 50)
        yy = c3 + m3*(xs/1000.0) if key == "dsc" else None
        if key == "dsc":
            ax.plot(xs, yy, "k--", lw=1.0, zorder=4)
        else:
            d3 = np.array([scalars(a*1e-6, f, s)["Deff"]/D for _, a, s, f, _ in THREE])
            md, cd = np.polyfit(d3, y3, 1)
            ax.plot(xs, cd+md*xs, "k--", lw=1.0, zorder=4)
    ax.set_xlabel(xlab)
    ax.grid(True, ls="--", lw=0.4, color="0.92")
    ax.set_axisbelow(True)
    ax.text(tx, ty, tag, transform=ax.transAxes, fontsize=8.5, va="top", ha=ha)
AX[0].set_ylabel("$x_{fc}/D_e$")
AX[0].text(0.035, 0.055, "3D fit  $R^2$ = %.4f" % r2, transform=AX[0].transAxes,
           fontsize=8.0)

H = [Line2D([], [], ls="", marker="*", ms=11, color="#111111",
            label="3D Cartesian, finite-thickness profiles"),
     Line2D([], [], ls="", marker="P", ms=8, color="#111111",
            label="top hat: 3D Cartesian"),
     Line2D([], [], ls="", marker="P", ms=8, color=CT, label="top hat: 2D L5"),
     Line2D([], [], ls="", marker="s", ms=5, color=CW,
            label="2D L5, wall-attached family"),
     Line2D([], [], ls="", marker="o", ms=5, color=CS,
            label="2D L5, shifted-tanh family"),
     Line2D([], [], ls="", marker="D", ms=5.4, color=CD, mec="w", mew=0.6,
            label=r"2D L5, iso-$\delta^*_c$ set, 533 $\mu$m (both families)"),
     Line2D([], [], ls="--", color="k", lw=1.0, label="linear fit to 3D values (a,c)"),
     Line2D([], [], ls="-", color="0.55", lw=0.7, label="experiment, $x_{fc}$")]
fig.legend(handles=H, loc="lower center", bbox_to_anchor=(0.5, 0.005), ncol=3,
           fontsize=7.6, handlelength=1.8, columnspacing=1.6, labelspacing=0.42)
fig.savefig(HERE/"thesis_xfc_iso_dstar.pdf")
fig.savefig(HERE/"thesis_xfc_iso_dstar.png", dpi=300)
print("wrote thesis_xfc_iso_dstar.pdf/png")
