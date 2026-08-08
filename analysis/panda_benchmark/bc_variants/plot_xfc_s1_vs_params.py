#!/usr/bin/env python3
"""Which inlet integral parameter controls the first shock cell?

Left column  : x_fc, the first-cell closure.
Right column : S_1 = x_2 - x_1, the first shock-cell spacing.
Rows         : the four candidate inlet parameters.

Both quantities are plotted against each parameter in turn, for the four 3D
reference cases and for the full 2D sweep on the L4 grid. A parameter that
controls the first cell should give a tight straight line in both columns; one
that does not should scatter.

x_fc is the maximum raw density gradient between the first expansion minimum
and the first downstream principal maximum. x_n are the peaks of the same mean
centreline density after a 0.05 D_e box filter, prominence >= 0.02, minimum
separation 0.5 D_e; S_1 = x_2 - x_1. Both follow the thesis definitions.
"""
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.signal import find_peaks

sys.path.insert(0, str(Path(__file__).resolve().parent))
from design_2d_sweep import scalars, D

HERE = Path(__file__).resolve().parent
OUT = HERE/"sweep_out"
RHO_J, T0 = 1.6413, 4.0e-4
CW, CS, CT, CH = "#b3251e", "#1f4e9c", "#6a3d9a", "#111111"

PROF = [("A0_tophat",     0.0, 3.0, "shifted", "hat"),
        ("Shift_s0",    170.0, 0.0, "shifted", "trunc"),
        ("Wall100",     100.0, 3.0, "wall",    "wall"),
        ("C1_wall200",  200.0, 3.0, "wall",    "wall"),
        ("Shift_s1",    170.0, 1.0, "shifted", "trunc"),
        ("Wall340",     340.0, 3.0, "wall",    "wall"),
        ("A1_a100",     100.0, 3.0, "shifted", "shift"),
        ("Shift_s2",    170.0, 2.0, "shifted", "trunc"),
        ("Wall470",     470.0, 3.0, "wall",    "wall"),
        ("B1_wall680",  679.5, 3.0, "wall",    "wall"),
        ("A2_a170",     170.0, 3.0, "shifted", "shift"),
        ("S4_a170_s4",  170.0, 4.0, "shifted", "shift"),
        ("Wall900",     900.0, 3.0, "wall",    "wall"),
        ("A3_a254",     254.0, 3.0, "shifted", "shift"),
        ("S5_a170_s5",  170.0, 5.0, "shifted", "shift"),
        ("Wall1150",   1150.0, 3.0, "wall",    "wall"),
        ("S6_a170_s6",  170.0, 6.0, "shifted", "shift"),
        ("S7_a170_s7",  170.0, 7.0, "shifted", "shift")]
COL = {"wall": CW, "shift": CS, "trunc": CT, "hat": CH}
MRK = {"wall": "s",  "shift": "o", "trunc": "^", "hat": "P"}

# 3D: thesis table 10_results_profile.tex
THREE = [("top hat",       0.0, 3.0, "shifted", 0.945, 1.244),
         ("wall 200",    200.0, 3.0, "wall",    0.930, 1.415),
         ("shifted 100", 100.0, 3.0, "shifted", 0.922, 1.508),
         ("baseline 254",254.0, 3.0, "shifted", 0.883, 1.300)]
EXP_XFC, EXP_S1 = 0.888, 1.359

PARAMS = [("dstar_c", r"displacement thickness  $\delta^*_c$  [$\mu$m]", 1e6),
          ("theta_c", r"momentum thickness  $\theta_c$  [$\mu$m]",       1e6),
          ("dw",      r"vorticity thickness  $\delta_\omega$  [$\mu$m]", 1e6),
          ("Deff",    r"effective diameter  $D_{\rm eff}/D_e=\sqrt{C_d}$", 1.0)]


def par(a, s, form, key):
    sc = scalars(a*1e-6, form, s)
    return np.sqrt(sc["Cd"]) if key == "Deff" else sc[key]


def refine(x, y, i):
    if i <= 0 or i >= len(y)-1:
        return x[i]
    den = y[i-1]-2*y[i]+y[i+1]
    return x[i] if den == 0 else x[i]+0.5*(y[i-1]-y[i+1])/den*(x[1]-x[0])


def boxfilt(y, x):
    w = max(1, int(round(0.05/(x[1]-x[0]))))
    k = np.ones(w)/w
    return np.convolve(y, k, mode="same")


def features(x, rho):
    """thesis extraction: x_fc, then the x_n peak train, S_1 = x_2 - x_1"""
    e = np.where((x >= 0.35) & (x <= 1.0))[0]
    i0 = e[np.nanargmin(rho[e])]
    c = np.where((np.arange(x.size) > i0) & (x <= 1.7))[0]
    i1 = c[np.nanargmax(rho[c])]
    g = np.gradient(rho, x)
    s = np.arange(i0, i1+1)
    xfc = refine(x, g, s[np.nanargmax(g[s])])
    sm = boxfilt(rho, x)
    dist = max(1, int(round(0.5/(x[1]-x[0]))))
    pk, _ = find_peaks(np.where(x > xfc, sm, -np.inf), prominence=0.02,
                       distance=dist)
    xn = np.array([refine(x, sm, i) for i in pk[:5]])
    S1 = xn[1]-xn[0] if len(xn) >= 2 else np.nan
    return xfc, (xn[0] if len(xn) else np.nan), S1


def load(nm):
    for pat in ("XFC_%s__L4.npz", "XFC_%s_L4.npz"):
        f = OUT/(pat % nm)
        if f.exists():
            break
    else:
        return None
    d = np.load(f)
    ks = sorted([k[:-2] for k in d.files if k.endswith("_t")],
                key=lambda q: int("".join(c for c in q if c.isdigit())))
    dx = float(d[ks[0]+"_t"][1])
    M = np.array([d[k+"_rho"] for k in ks])/RHO_J
    x = (np.arange(M.shape[1])+0.5)*dx/D
    return features(x, M[-1])


rows = []
for nm, a, s, form, fam in PROF:
    r = load(nm)
    if r is None:
        continue
    xfc, x1, S1 = r
    rows.append(dict(nm=nm, fam=fam, a=a, s=s, form=form, xfc=xfc, x1=x1, S1=S1))
print("2D L4:  %-14s %8s %8s %8s" % ("case", "x_fc", "x_1", "S_1"))
for r in rows:
    print("        %-14s %8.4f %8.4f %8.4f" % (r["nm"], r["xfc"], r["x1"], r["S1"]))

good = [r for r in rows if np.isfinite(r["S1"])]
print("\n2D L4 cases with a usable S_1: %d of %d" % (len(good), len(rows)))


def fit(x, y):
    m, c = np.polyfit(x, y, 1)
    r = y-(c+m*x)
    r2 = 1-(r**2).sum()/((y-y.mean())**2).sum()
    return m, c, r2


print("\n%-10s | %-28s | %-28s" % ("parameter", "x_fc   (3D n=4 / 2D n=%d)" % len(good),
                                   "S_1    (3D n=4 / 2D n=%d)" % len(good)))
STAT = {}
for key, lab, sc in PARAMS:
    p3 = np.array([par(a, s, f, key)*sc for _, a, s, f, *_ in THREE])
    y3f = np.array([q[4] for q in THREE]); y3s = np.array([q[5] for q in THREE])
    p2 = np.array([par(r["a"], r["s"], r["form"], key)*sc for r in good])
    y2f = np.array([r["xfc"] for r in good]); y2s = np.array([r["S1"] for r in good])
    STAT[key] = dict(p3=p3, y3f=y3f, y3s=y3s, p2=p2, y2f=y2f, y2s=y2s,
                     f3=fit(p3, y3f), s3=fit(p3, y3s),
                     f2=fit(p2, y2f), s2=fit(p2, y2s))
    print("%-10s | R2 3D %.4f   2D %.4f      | R2 3D %.4f   2D %.4f"
          % (key, STAT[key]["f3"][2], STAT[key]["f2"][2],
             STAT[key]["s3"][2], STAT[key]["s2"][2]))

# ---------------------------------------------------------------- figure
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})
fig, axs = plt.subplots(len(PARAMS), 2, figsize=(10.4, 11.4),
                        gridspec_kw={"hspace": 0.30, "wspace": 0.20})

for irow, (key, lab, sc) in enumerate(PARAMS):
    S = STAT[key]
    for icol, (yk3, yk2, fk3, fk2, ylab, expv) in enumerate(
            [("y3f", "y2f", "f3", "f2", "$x_{fc}/D_e$", EXP_XFC),
             ("y3s", "y2s", "s3", "s2", "$S_1/D_e$", EXP_S1)]):
        ax = axs[irow, icol]
        ax.axhline(expv, color="0.62", lw=0.9, zorder=1)
        if irow == 0:
            ax.text(0.985, expv, " experiment", transform=ax.get_yaxis_transform(),
                    ha="right", va="bottom", fontsize=6.8, color="0.4")
        xx = np.linspace(min(S["p2"].min(), S["p3"].min()),
                         max(S["p2"].max(), S["p3"].max()), 20)
        m2, c2, r22 = S[fk2]
        m3, c3, r23 = S[fk3]
        ax.plot(xx, c2+m2*xx, color="0.45", lw=1.1, ls="-", zorder=2)
        ax.plot(xx, c3+m3*xx, color="k", lw=1.2, ls="--", zorder=2)
        for r, p, y in zip(good, S["p2"], S[yk2]):
            ax.plot(p, y, MRK[r["fam"]], color=COL[r["fam"]],
                    ms=7.0 if r["fam"] == "hat" else 5.4, mec="w", mew=0.7,
                    zorder=5)
        ax.plot(S["p3"], S[yk3], "*", color="k", ms=13, zorder=7)
        ax.set_ylabel(ylab)
        if irow == len(PARAMS)-1:
            ax.set_xlabel(lab)
        else:
            ax.set_xlabel(lab, fontsize=8, labelpad=2)
        ax.grid(True, ls="--", lw=0.5, color="0.93"); ax.set_axisbelow(True)
        ax.text(0.035, 0.055,
                "3D  $R^2$ = %.4f\n2D L4  $R^2$ = %.4f" % (r23, r22),
                transform=ax.transAxes, fontsize=7.2, va="bottom",
                bbox=dict(fc="w", ec="0.85", lw=0.5, pad=2.2))

axs[0, 0].legend(handles=[
    Line2D([], [], ls="", marker="*", color="k", ms=12, label="3D reference (4)"),
    Line2D([], [], ls="", marker="P", color=CH, ms=7, mec="w", label="2D L4 top hat"),
    Line2D([], [], ls="", marker="s", color=CW, ms=5.4, mec="w", label="2D L4 wall tanh"),
    Line2D([], [], ls="", marker="o", color=CS, ms=5.4, mec="w", label="2D L4 shifted tanh"),
    Line2D([], [], ls="", marker="^", color=CT, ms=5.4, mec="w", label="2D L4 truncation"),
    Line2D([], [], color="k", lw=1.2, ls="--", label="fit, 3D"),
    Line2D([], [], color="0.45", lw=1.1, label="fit, 2D L4")],
    fontsize=6.8, frameon=False, loc="upper right", handlelength=1.6,
    labelspacing=0.30, borderaxespad=0.4)

fig.tight_layout()
fig.savefig(HERE/"xfc_s1_vs_params.png", dpi=250)
fig.savefig(HERE/"xfc_s1_vs_params.pdf")
print("\nwrote xfc_s1_vs_params.png/pdf")
