#!/usr/bin/env python3
"""Landmark peak-registration variant of Figure 3.15.

The four nozzle-profile centreline densities share the near-nozzle
structure but their shock-cell spacing differs slightly, so the extrema
drift apart with x (cumulative phase shift). A single constant shift
cannot fix this. Here the x-axis of every non-baseline curve is warped by
a MONOTONE map psi that sends its n-th extremum onto the baseline's n-th
extremum (alternating maxima and minima used as landmarks), with PCHIP
interpolation between landmarks and the endpoints (0,0), (x_last,x_last)
anchored. Amplitudes are left untouched: only the x-coordinate is
deformed, so the operation moves extrema into vertical alignment without
rescaling any value. This is a diagnostic registration, not a corrected
prediction.

  axis_4profile_aligned : (a) raw profiles, (b) same after warping every
  non-baseline curve onto the baseline extrema.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.signal import find_peaks
from scipy.interpolate import PchipInterpolator

HERE = Path(__file__).resolve().parent
D_E, RHO_J = 0.0254, 1.6413
XMAX = 8.0
CB = "#aab0b7"
VARIANTS = [
    ("tophat", "Top hat", "#0b6ef5", (0, (5, 2))),
    ("walltanh", "Wall tanh", "#2ca02c", (0, (3, 1.4, 1, 1.4))),
    ("shift100", "Thin shifted tanh", "#ff6a00", (0, (1, 1.4))),
]


def boxfilt(a, x, wD=0.05):
    w = max(3, int(round(wD / (x[1] - x[0]))))
    if w % 2 == 0:
        w += 1
    out = np.convolve(a, np.ones(w) / w, mode="same")
    out[:w // 2] = a[:w // 2]; out[-(w // 2):] = a[-(w // 2):]
    return out


def refine(x, y, i):
    if i <= 0 or i >= len(y) - 1:
        return x[i]
    d = y[i - 1] - 2 * y[i] + y[i + 1]
    return x[i] if d == 0 else x[i] + 0.5 * (y[i - 1] - y[i + 1]) / d * (x[1] - x[0])


def landmarks(x, y, prom=0.02, nmax=6):
    """Principal compression maxima x_1..x_n of the mean density (the
    thesis x_n metric): clean, well-separated, one per shock cell. These
    are the registration landmarks."""
    dist = max(1, int(round(0.5 / (x[1] - x[0]))))
    # first-cell closure: max positive gradient in [0.35, 1.65]
    exp = np.where((x >= 0.35) & (x <= 1.0))[0]
    imin = exp[np.nanargmin(y[exp])]
    comp = np.where((np.arange(x.size) > imin) & (x <= 1.65))[0]
    imax = comp[np.nanargmax(y[comp])]
    xfc = refine(x, np.gradient(y, x), np.arange(imin, imax + 1)[
        np.nanargmax(np.gradient(y, x)[imin:imax + 1])])
    pk, _ = find_peaks(np.where(x > xfc, y, -np.inf), prominence=prom,
                       distance=dist)
    return np.array([refine(x, y, i) for i in pk[:nmax]])


def load(v):
    d = np.load(HERE / f"AXIS8_{v}.npz")
    x = d["x"]
    return x, boxfilt(d["DensityMEAN"] / RHO_J, x)


def warp_to(xref_land, xsrc_land, x, y, xend):
    """f_tilde(x) = f(psi(x)) with psi monotone, psi(xref_n)=xsrc_n."""
    m = min(len(xref_land), len(xsrc_land))
    nodes_ref = np.concatenate(([0.0], xref_land[:m], [xend]))
    nodes_src = np.concatenate(([0.0], xsrc_land[:m], [xend]))
    order = np.argsort(nodes_ref)
    nodes_ref, nodes_src = nodes_ref[order], nodes_src[order]
    keep = np.concatenate(([True], np.diff(nodes_ref) > 1e-9))
    psi = PchipInterpolator(nodes_ref[keep], nodes_src[keep])
    return np.interp(psi(x), x, y), m


xb, rb = load("baseline")
ext_b = landmarks(xb, rb)
exp = np.load(HERE.parent / "panda_m142den_axis_EPAPS.npz")

plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})
fig, (axA, axB) = plt.subplots(
    2, 1, figsize=(7.6, 5.4), sharex=True,
    gridspec_kw={"height_ratios": [1.0, 1.0], "hspace": 0.08})
for ax in (axA, axB):
    ax.set_axisbelow(True); ax.grid(True, ls="--", lw=0.5, color="0.88")
    ax.tick_params(length=2.5)

mb = (xb > 0.05) & (xb <= XMAX)
# panel (a): raw
axA.plot(xb[mb], rb[mb], color=CB, lw=2.0, zorder=3)
# panel (b): aligned (baseline unchanged)
axB.plot(xb[mb], rb[mb], color=CB, lw=2.0, zorder=3)
print("baseline extrema:", np.round(ext_b, 3))
for v, lab, c, ls in VARIANTS:
    xv, rv = load(v)
    mv = (xv > 0.05) & (xv <= XMAX)
    axA.plot(xv[mv], rv[mv], color=c, lw=1.4, ls=ls, zorder=4)
    ext_v = landmarks(xv, rv)
    rw, m = warp_to(ext_b, ext_v, xv, rv, xv[-1])
    axB.plot(xv[mv], rw[mv], color=c, lw=1.4, ls=ls, zorder=4)
    print("%-9s extrema:" % v, np.round(ext_v, 3), " -> warped with %d landmarks" % m)
# experiment: raw markers in (a); warped markers in (b)
rexp = exp["rho"]; xexp = exp["x"]
ext_e = landmarks(xexp, boxfilt(rexp, xexp, 0.06), prom=0.03)
axA.plot(xexp, rexp, "o", mfc="w", mec="k", mew=0.9, ms=4.2, zorder=5)
# warp experiment x-positions onto baseline extrema
m = min(len(ext_b), len(ext_e))
if m >= 1:
    nodes_src = np.concatenate(([0.0], ext_e[:m], [xexp[-1]]))
    nodes_ref = np.concatenate(([0.0], ext_b[:m], [xexp[-1]]))
    o = np.argsort(nodes_src)
    inv = PchipInterpolator(nodes_src[o], nodes_ref[o])   # src -> ref
    xexp_w = inv(xexp)
else:
    xexp_w = xexp
axB.plot(xexp_w, rexp, "o", mfc="w", mec="k", mew=0.9, ms=4.2, zorder=5)

for ax, lab in ((axA, "(a) raw"), (axB, "(b) extrema registered to baseline")):
    ax.set_ylim(0.3, 1.55)
    ax.set_ylabel(r"$\langle\rho\rangle/\rho_j$")
    ax.text(0.14, 0.44, lab, fontsize=9, va="top")
for xv in ext_b:
    axB.axvline(xv, color="0.6", lw=0.5, ls=":", zorder=1)
axB.set_xlim(0, XMAX)
axB.set_xlabel(r"$x/D_e$")
handles = [
    Line2D([], [], marker="o", ls="", mfc="w", mec="k", ms=4.2,
           label=r"$M_j=1.42$ experiment"),
    Line2D([], [], color=CB, lw=2.0, label="Baseline (reference)"),
] + [Line2D([], [], color=c, lw=1.4, ls=ls, label=lab)
     for _, lab, c, ls in VARIANTS]
axA.legend(handles=handles, fontsize=8, frameon=False, ncol=3,
           loc="upper center", bbox_to_anchor=(0.5, 1.03),
           handlelength=2.0, columnspacing=1.3, handletextpad=0.5)
fig.subplots_adjust(left=0.10, right=0.975, top=0.93, bottom=0.095)
for ext in ("png", "pdf"):
    fig.savefig(HERE / f"axis_4profile_aligned.{ext}", dpi=250)
print("saved axis_4profile_aligned.png/pdf")
