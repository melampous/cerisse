#!/usr/bin/env python3
"""Compare the five formulations used in Section 3.5.3.

Four of them form a nominal formulation matrix and are encoded by it rather than by
four arbitrary hues:

    colour      flux path        LLF (blue)  vs  AFD + HLLC (vermillion)
    line style  reconstruction   WENO-Z5 (solid) vs TENO5 (dashed)

The fifth, the skew-symmetric split with JST scalar dissipation, sits outside
that design: it is a central formulation with an explicit artificial
dissipation term rather than a flux path combined with a shock-capturing
reconstruction, so neither factor applies to it. It carries its own colour
(bluish green) and is drawn with a heavier line so that it is not read as a
member of either pair.

Top     centreline mean-density profiles and successive maxima
Bottom  density at each maximum and its position relative to experiment

The feature extraction follows the common Chapter 3 definition. Numerical
profiles are smoothed with a 0.05 D_e moving box filter. Successive maxima
have a minimum prominence of 0.02 and a minimum separation of 0.5 D_e, and
their positions are refined by a three-point quadratic interpolation. The
experimental maxima are extracted from the unsmoothed sampled profile using
the same prominence, separation, and interpolation procedure.

The statistics windows are [6.5, 10.16] ms and [6.5, 9.30] ms for the LLF
WENO-Z5 and TENO5 cases, respectively, and [15.0, 17.6] ms for the two
AFD--HLLC cases and for the skew-symmetric case. This difference is considered
when interpreting the density values at the maxima.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.signal import find_peaks

HERE = Path(__file__).resolve().parent
PB = HERE.parent
OUT = HERE/"sweep_out"
RHO_J = 1.6413
D = 0.0254

C_LLF, C_HLLC, C_SKEW = "#0072B2", "#D55E00", "#009E73"


def native(x, r):
    """One sample per real cell.

    The axis line is stored on a uniform level-5 abscissa, but away from the
    finest patches the value comes from a coarser cell that has simply been
    replicated across 2, 4 or 8 output slots (level 4, 3, 2 - 198.4, 396.9 and
    793.7 um against the 99.2 um of the store). Drawing that as a curve on the
    fine abscissa renders a staircase and implies a resolution the data does
    not have. Collapsing each run of equal values to its midpoint gives one
    point per cell that actually exists, at the resolution it actually has.
    """
    m = np.isfinite(r)
    x, r = x[m], r[m]
    edge = np.r_[0, np.flatnonzero(np.diff(r) != 0) + 1, r.size]
    xc = np.array([0.5*(x[a] + x[b-1]) for a, b in zip(edge[:-1], edge[1:])])
    rc = r[edge[:-1]]
    return xc, rc


def from_stats3d(f):
    """Ring 0 hugs the axis but runs dry at 2 D_e; fill from rings 1 and 2."""
    z = np.load(PB/f)
    rings = z["ax_DensityMEAN"]/RHO_J
    r = rings[0].copy()
    for k in (1, 2):
        gap = ~np.isfinite(r)
        r[gap] = rings[k][gap]
    return z["xx"], r


def from_axpack(f):
    z = np.load(OUT/f)
    ts = sorted(float(k[1:].split("_")[0]) for k in z.files
                if k.endswith("_DensityMEAN"))
    r = z["t%.3f_DensityMEAN" % ts[-1]]/RHO_J
    return z["xx"], r


def from_axframe(f):
    z = np.load(OUT/f)
    return z["x"]/D, z["ax_DensityMEAN"]/RHO_J


CASES = [
    dict(lab="LLF–WENO-Z5",  c=C_LLF,  ls="-",  win="[6.5, 10.16] ms",
         d=lambda: from_stats3d("STATS3D_baseline.npz")),
    dict(lab="LLF–TENO5",    c=C_LLF,  ls="--", win="[6.5, 9.30] ms",
         d=lambda: from_stats3d("STATS3D_teno.npz")),
    dict(lab="AFD–HLLC–WENO-Z5", c=C_HLLC, ls="-",  win="[15.0, 17.6] ms",
         d=lambda: from_axpack("AX_L5_afdhllc_weno.npz")),
    dict(lab="AFD–HLLC–TENO5",   c=C_HLLC, ls="--", win="[15.0, 17.6] ms",
         d=lambda: from_axframe("AX_L5_s3d_afdhllc_plt15761.npz")),
    dict(lab="Skew–JST", c=C_SKEW, ls="-", lw=2.2,
         win="[15.0, 17.6] ms",
         d=lambda: from_axframe("AX_L5_s3d_skewjst_plt15594.npz")),
]

e = np.load(PB/"panda_m142den_axis_EPAPS.npz")
ex, er = e["x"], e["rho"]


def box_filter(values, x, width=0.05):
    """Moving box filter with the nominal width used in Section 3.4."""
    count = max(3, int(round(width / (x[1] - x[0]))))
    if count % 2 == 0:
        count += 1
    result = np.convolve(values, np.ones(count) / count, mode="same")
    result[:count // 2] = values[:count // 2]
    result[-(count // 2):] = values[-(count // 2):]
    return result


def refine(x, values, index):
    """Three-point quadratic refinement of a selected sample location."""
    if index <= 0 or index >= values.size - 1:
        return float(x[index])
    denominator = values[index - 1] - 2.0 * values[index] + values[index + 1]
    if denominator == 0.0:
        return float(x[index])
    offset = 0.5 * (values[index - 1] - values[index + 1]) / denominator
    return float(x[index] + offset * (x[1] - x[0]))


def closure_location(x, raw_density):
    expansion = np.where((x >= 0.35) & (x <= 1.0))[0]
    first_minimum = expansion[np.nanargmin(raw_density[expansion])]
    recompression = np.where(
        (np.arange(x.size) > first_minimum) & (x <= 1.7)
    )[0]
    first_maximum = recompression[np.nanargmax(raw_density[recompression])]
    gradient = np.gradient(raw_density, x)
    interval = np.arange(first_minimum, first_maximum + 1)
    selected = interval[np.nanargmax(gradient[interval])]
    return refine(x, gradient, selected)


def cells(x, raw_density, numerical, n=6):
    """Successive mean-density maxima using the Chapter 3 procedure."""
    values = box_filter(raw_density, x) if numerical else raw_density
    x_fc = closure_location(x, raw_density)
    distance = max(1, int(round(0.5 / (x[1] - x[0]))))
    peaks, _ = find_peaks(
        np.where(x > x_fc, values, -np.inf),
        prominence=0.02,
        distance=distance,
    )
    peaks = peaks[:n]
    positions = np.array([refine(x, values, i) for i in peaks])
    return positions, values[peaks]


epx, epr = cells(ex, er, numerical=False)
print("experiment cells   " + "  ".join("%.2f@%.2f" % (p, q)
                                        for q, p in zip(epx, epr)))
for c in CASES:
    c["fx"], c["fr"] = c["d"]()
    c["x"], c["r"] = native(c["fx"], c["fr"])
    c["px"], c["pr"] = cells(c["fx"], c["fr"], numerical=True)
    c["res"] = np.interp(ex, c["x"], c["r"], left=np.nan, right=np.nan) - er
    print("%-15s " % c["lab"] + "  ".join("%.2f@%.2f" % (p, q)
                                          for q, p in zip(c["px"], c["pr"]))
          + "   rms residual %.4f" % np.sqrt(np.nanmean(c["res"]**2)))

# ------------------------------------------------------------------ figure
plt.rcParams.update({"font.size": 12.2, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True,
                     "legend.frameon": False})
# Keep the source canvas close to the final thesis aspect and size.  The
# previous 12-inch canvas was reduced substantially in LaTeX, leaving text at
# only about 5 pt in the final document.
fig = plt.figure(figsize=(9.2, 6.15))
gs = fig.add_gridspec(2, 2, height_ratios=[1.30, 1.0], hspace=0.22,
                      wspace=0.22, left=0.078, right=0.985,
                      top=0.88, bottom=0.095)
a = fig.add_subplot(gs[0, :])
cc = fig.add_subplot(gs[1, 0])
dd = fig.add_subplot(gs[1, 1])

# centrelines overlaid, with the detected maxima marked on the curves so the
# indexing of the two lower panels can be checked against them by eye. The
# window starts below the first maximum (x ~ 1.15): cropping at 1.5 D_e hid
# n = 1 and made the visible peaks look one index out of step.
a.plot(ex, er, "o", color="k", ms=4.2, mfc="none", mew=1.0, zorder=6)
a.plot(epx, epr, "o", color="k", ms=9, mfc="none", mew=1.6, zorder=8)
for c in CASES:
    a.plot(c["x"], c["r"], color=c["c"], ls=c["ls"],
           lw=c.get("lw", 1.6), zorder=5)
    a.plot(c["px"], c["pr"], ls="", marker="v", color=c["c"], ms=7,
           mec="w", mew=0.8, zorder=7)
for i, px in enumerate(epx, start=1):
    a.text(px, 1.515, "$n$=%d" % i, ha="center", va="center", fontsize=11.2,
           color="0.30")
a.set_xlim(0.75, 7.2)
a.set_ylim(0.30, 1.56)
a.set_xlabel("$x/D_e$", labelpad=2)
a.set_ylabel(r"$\langle\rho\rangle/\rho_j$ on the axis")

# the decay ladder
n = np.arange(1, len(epr)+1)
cc.plot(n, epr, "o-", color="k", ms=5.4, mfc="none", mew=1.1, lw=1.2, zorder=6)
for c in CASES:
    k = np.arange(1, len(c["pr"])+1)
    cc.plot(k, c["pr"], color=c["c"], ls=c["ls"], lw=c.get("lw", 1.7),
            marker="s", ms=4.4, zorder=5)
cc.set_xlabel(r"index $n$ of successive mean-density maxima", labelpad=2)
cc.set_ylabel(r"$\langle\rho\rangle_{\mathrm{max},n}/\rho_j$")
cc.set_xticks(np.arange(1, 7))
# headroom at the top so the panel label does not sit on the n=1 markers
cc.set_ylim(0.83, 1.50)

# cell position error. The experiment is the datum, so its own points lie at
# zero by construction; they are drawn with the same open circles used in the
# two panels above rather than left implicit in the axis line, so the reference
# is a plotted series here as well.
dd.axhline(0.0, color="0.72", lw=0.9, zorder=2)
dd.plot(n, np.zeros_like(n, dtype=float), "o-", color="k", ms=5.4,
        mfc="none", mew=1.1, lw=1.2, zorder=6)
for c in CASES:
    k = min(len(c["px"]), len(epx))
    dd.plot(np.arange(1, k+1), c["px"][:k]-epx[:k], color=c["c"], ls=c["ls"],
            lw=c.get("lw", 1.7), marker="s", ms=4.4, zorder=5)
dd.set_xlabel(r"index $n$ of successive mean-density maxima", labelpad=2)
dd.set_ylabel(r"$\left(x_n-x_n^{\mathrm{exp}}\right)/D_e$")
dd.set_xticks(np.arange(1, 7))

for ax in (a, cc, dd):
    ax.grid(True, ls="--", lw=0.5, color="0.93")
    ax.set_axisbelow(True)

a.text(0.012, 0.965, "(a)", transform=a.transAxes, ha="left", va="top")
cc.text(0.018, 0.965, "(b)", transform=cc.transAxes, ha="left", va="top")
dd.text(0.018, 0.965, "(c)", transform=dd.transAxes, ha="left", va="top")

H = [Line2D([], [], color="k", ls="", marker="o", mfc="none", mew=1.0, ms=5.4,
            label=r"Experiment ($M_j=1.42$)")]
H += [Line2D([], [], color=c["c"], ls=c["ls"], lw=c.get("lw", 1.9),
             label=c["lab"]) for c in CASES]
# six entries: three per row keeps every label on one line at this width
fig.legend(handles=H, loc="upper center", bbox_to_anchor=(0.53, 1.001),
           ncol=3, fontsize=11.3, handlelength=2.1, columnspacing=1.30,
           labelspacing=0.35, frameon=False)

fig.savefig(HERE/"scheme_differences.png", dpi=250)
fig.savefig(HERE/"scheme_differences.pdf")
print("\nwrote scheme_differences.png/pdf")
