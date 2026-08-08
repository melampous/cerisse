#!/usr/bin/env python3
"""Two scalars per scheme that can be compared with the measurement directly.

The centreline density of an underexpanded jet is a decaying oscillation. Two
numbers describe it and both are measured quantities, not fitted shapes:

    L_s   the shock-cell spacing, from the mean separation of the density
          maxima;
    L_d   the decay length of the cell train, from a straight-line fit to
          log A(x), where A is the rolling envelope

              A(x) = 0.5 [ max(rho) - min(rho) ] over a window of one cell
                     length centred on x.

The envelope is taken with a rolling window rather than from peak-trough pairs
because peak detection on the AFD-HLLC cases picks up a secondary maximum near
the first recompression and reports a first-cell amplitude of 0.06 instead of
0.44, which then propagates into every derived number. The rolling form has no
detection step and is immune to that.

The comparison window is 1.4 < x/D_e < 6.4: it starts one cell length inside
the measurement so the window is always full, and ends one half-window before
the measurement stops at 6.92.

Reading the result
------------------
Agreement with the measurement at one grid does not establish that a scheme is
right. The cell train in the experiment decays because turbulent mixing
thickens the shear layer; in an under-resolved computation it also decays
because of numerical dissipation. A scheme whose dissipation happens to match
the missing mixing will reproduce L_d for a reason that has nothing to do with
the physics, and will stop doing so as soon as the grid changes. Separating the
two needs the same schemes on a finer grid: a scheme that is right keeps its
L_d under refinement, one that is compensating drifts.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

HERE = Path(__file__).resolve().parent
RHO_J = 1.6413
WIN = 1.15                      # rolling window, one cell length
QLO, QHI = 1.40, 6.40

CASES = [
    ("experiment", "Panda & Seasholtz (1999)", None, "k", "-", 0.0),
    ("llfweno", "LLF + WENO-Z5",
     "bc_variants/fields/SLM_llfweno_DensityMEAN.npz", "#0072B2", "-", 1.5),
    ("llfteno", "LLF + TENO5",
     "bc_variants/fields/SLM_llfteno_DensityMEAN.npz", "#0072B2", "--", 1.5),
    ("hllcweno", "AFD-HLLC + WENO-Z5",
     "bc_variants/fields/SLM_hllcweno_DensityMEAN.npz", "#D55E00", "-", 1.5),
    ("hllcteno", "AFD-HLLC + TENO5",
     "bc_variants/fields/SLM_hllcteno_DensityMEAN.npz", "#D55E00", "--", 1.5),
    ("skewjst", "skew-symmetric + JST",
     "bc_variants/sweep_out/SL_L5_s3d_skewjst_DensityMEAN.npz",
     "#009E73", "-", 2.1)]


def centreline(f):
    z = np.load(HERE/f)
    j = int(np.argmin(np.abs(z["y"])))
    return z["x"], z["field"][:, j]/RHO_J


def rolling_env(x, r, w=WIN):
    return np.array([0.5*np.ptp(r[(x >= q-w/2) & (x <= q+w/2)]) for q in x])


def spacing(x, r):
    m = (x >= 0.6) & (x <= 6.95)
    xx, rr = x[m], r[m]
    i, _ = find_peaks(rr, prominence=0.06,
                      distance=max(int(0.45/np.median(np.diff(xx))), 1))
    return (np.mean(np.diff(xx[i])) if len(i) > 1 else np.nan), xx[i]


e = np.load(HERE/"panda_m142den_axis_EPAPS.npz")
Q = np.arange(QLO, QHI+1e-9, 0.02)
R = {}
print("%-24s %9s %9s %9s %9s %9s"
      % ("source", "L_s", "L_d", "A(1.5)", "A(5.5)", "env rms"))
for key, name, f, col, ls, lw in CASES:
    if key == "experiment":
        x, r = e["x"], e["rho"]
    else:
        x, r = centreline(f)
    A = np.interp(Q, x, rolling_env(x, r))
    p = np.polyfit(Q, np.log(np.clip(A, 1e-6, None)), 1)
    Ls, _ = spacing(x, r)
    R[key] = dict(A=A, Ld=-1.0/p[0], Ls=Ls, name=name, col=col, ls=ls, lw=lw)
    d = A - R["experiment"]["A"] if key != "experiment" else np.zeros_like(A)
    R[key]["rms"] = float(np.sqrt(np.mean(d*d)))
    print("%-24s %9.3f %9.3f %9.3f %9.3f %9.4f"
          % (name, Ls, R[key]["Ld"], np.interp(1.5, Q, A),
             np.interp(5.5, Q, A), R[key]["rms"]))

E = R["experiment"]
print("\nrelative to the measurement")
print("%-24s %10s %10s" % ("scheme", "L_s", "L_d"))
for key, *_ in CASES[1:]:
    print("%-24s %+9.1f%% %+9.1f%%"
          % (R[key]["name"], 100*(R[key]["Ls"]/E["Ls"]-1),
             100*(R[key]["Ld"]/E["Ld"]-1)))

plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True,
                     "legend.frameon": False})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.6, 4.4),
                             gridspec_kw=dict(width_ratios=[1.6, 1.0]))
fig.subplots_adjust(left=0.078, right=0.985, top=0.925, bottom=0.135,
                    wspace=0.215)

for key, *_ in CASES:
    q = R[key]
    if key == "experiment":
        a1.plot(Q, q["A"], "o", ms=3.4, mfc="w", mec="k", mew=0.9, ls="",
                zorder=8, label=q["name"])
    else:
        a1.plot(Q, q["A"], color=q["col"], ls=q["ls"], lw=q["lw"], zorder=5,
                label="%s   $L_d$ = %.2f" % (q["name"], q["Ld"]))
a1.set_yscale("log")
a1.set_xlim(QLO, QHI)
a1.set_ylim(0.08, 0.62)
a1.set_xlabel("$x/D_e$")
a1.set_ylabel(r"cell amplitude $A$   [$\overline{\rho}/\rho_j$]")
a1.set_title("(a)  shock-cell envelope, one-cell rolling window",
             fontsize=10, pad=5)
a1.legend(loc="lower left", fontsize=8.4, handlelength=2.2)

for key, *_ in CASES[1:]:
    q = R[key]
    a2.plot(100*(q["Ls"]/E["Ls"]-1), 100*(q["Ld"]/E["Ld"]-1), "o", ms=9,
            mfc=q["col"], mec="k", mew=0.8, ls="", zorder=6,
            alpha=1.0 if q["ls"] == "-" else 0.55)
    off = (9, 9) if key == "skewjst" else (9, -3)
    a2.annotate(q["name"].replace(" + ", "+"),
                (100*(q["Ls"]/E["Ls"]-1), 100*(q["Ld"]/E["Ld"]-1)),
                textcoords="offset points", xytext=off, fontsize=8.2)
a2.plot([0], [0], "*", ms=16, mfc="w", mec="k", mew=1.0, zorder=7)
a2.annotate("measurement", (0, 0), textcoords="offset points",
            xytext=(-14, -19), fontsize=8.6, ha="center")
a2.axhline(0, color="0.6", lw=0.8, ls=":")
a2.axvline(0, color="0.6", lw=0.8, ls=":")
a2.set_xlim(-16, 6)
a2.set_ylim(-30, 360)
a2.set_xlabel(r"cell spacing $L_s$, error [%]")
a2.set_ylabel(r"decay length $L_d$, error [%]")
a2.set_title("(b)  both errors at once", fontsize=10, pad=5)
for ax in (a1, a2):
    ax.grid(True, ls="--", lw=0.5, color="0.93")
    ax.set_axisbelow(True)

fig.savefig(HERE/"scheme_metrics.png", dpi=250)
fig.savefig(HERE/"scheme_metrics.pdf")
print("\nwrote scheme_metrics.png/pdf")
