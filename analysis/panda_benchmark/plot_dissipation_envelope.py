#!/usr/bin/env python3
"""Numerical dissipation of the five Euler formulations, measured on the jet.

The centreline mean density is a decaying oscillation whose envelope is the
only quantity in this data set that integrates the dissipation of the
convective operator over the whole shock-cell train. Normalising the envelope
at one station removes the level offset left by the unequal statistics windows
and leaves the decay, which is what differs between the formulations.

    E(x)   = |hilbert(rho_bar - boxfilter(rho_bar, L_s))|, smoothed by the
             same box; L_s = 1.0 D_e is the shock-cell spacing, so the filter
             annihilates the fundamental exactly
    E_n(x) = E(x) / E(2 D_e)
    x_half = first x >= 2 D_e at which E_n falls to 0.5

x = 2 D_e is the anchor because all five cases still have 198 um cells there;
further downstream the AMR ladder coarsens identically for all of them, so
x_half is a comparative measure on a shared mesh and not a physical decay
length. It must not be quoted as one.

Panel (c) is the point of the figure. Four of the five form a two-by-two
design in flux path and reconstruction, so the effect of each factor can be
separated instead of asserted. The skew-symmetric case sits outside the design
and is drawn detached.

Caveats that belong in the caption
----------------------------------
The LLF pair was averaged over [6.5, 10.16] and [6.5, 9.30] ms, the AFD-HLLC
pair and the skew case over [15.0, 17.6] ms. Normalising at 2 D_e removes the
resulting level offset but not any difference in convergence; the bands on the
three cases that have a frame ensemble show what that scatter amounts to, and
the LLF pair has one frame each and therefore no band.

Agreement with the measurement is not a ranking of accuracy. The measured
envelope lies below four of the five: the least dissipative formulation is the
furthest from it. Whether the closest one is right for the right reason cannot
be settled at a single resolution.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.signal import hilbert

HERE = Path(__file__).resolve().parent
SW = HERE/"bc_variants"/"sweep_out"
RHO_J, D = 1.6413, 0.0254
LS = 1.00                      # shock-cell spacing used as the filter width
XA = 2.00                      # normalisation anchor
C_LLF, C_HLLC, C_SKEW = "#0072B2", "#D55E00", "#009E73"


def native(x, r):
    m = np.isfinite(r)
    x, r = x[m], r[m]
    e = np.r_[0, np.flatnonzero(np.diff(r) != 0)+1, r.size]
    return (np.array([0.5*(x[a]+x[b-1]) for a, b in zip(e[:-1], e[1:])]),
            r[e[:-1]])


def from_stats3d(f):
    z = np.load(HERE/f)
    rings = z["ax_DensityMEAN"]/RHO_J
    r = rings[0].copy()
    for k in (1, 2):
        g = ~np.isfinite(r)
        r[g] = rings[k][g]
    return z["xx"], r


def from_axpack(f, which=-1):
    z = np.load(SW/f)
    ts = sorted(float(k[1:].split("_")[0]) for k in z.files
                if k.endswith("_DensityMEAN"))
    return z["xx"], z["t%.3f_DensityMEAN" % ts[which]]/RHO_J


def from_axframe(f):
    z = np.load(SW/f)
    return z["x"]/D, z["ax_DensityMEAN"]/RHO_J


def box(v, n):
    n = n if n % 2 else n+1
    out = np.convolve(v, np.ones(n)/n, mode="same")
    out[:n//2] = v[:n//2]
    out[-(n//2):] = v[-(n//2):]
    return out


def envelope(x, r):
    xu = np.arange(0.60, min(x.max(), 9.90), 0.005)
    ru = np.interp(xu, *native(x, r))
    n = int(round(LS/0.005))
    E = box(np.abs(hilbert(ru - box(ru, n))), n)
    return xu, E/np.interp(XA, xu, E)


def x_half(xu, En):
    m = xu >= XA
    xx, ee = xu[m], En[m]
    k = np.flatnonzero(ee < 0.5)
    if not k.size:
        return np.nan
    i = k[0]
    if i == 0:
        return xx[0]
    return xx[i-1] + (0.5-ee[i-1])*(xx[i]-xx[i-1])/(ee[i]-ee[i-1])


CASES = [
    dict(k="llfweno",  lab="LLF + WENO-Z5",      c=C_LLF,  ls="-",  lw=1.7,
         flux=0, rec=0, load=lambda: from_stats3d("STATS3D_baseline.npz"),
         frames=None),
    dict(k="llfteno",  lab="LLF + TENO5",        c=C_LLF,  ls="--", lw=1.7,
         flux=0, rec=1, load=lambda: from_stats3d("STATS3D_teno.npz"),
         frames=None),
    dict(k="hllcweno", lab="AFD–HLLC + WENO-Z5", c=C_HLLC, ls="-",  lw=1.7,
         flux=1, rec=0, load=lambda: from_axpack("AX_L5_afdhllc_weno.npz"),
         frames=("pack", "AX_L5_afdhllc_weno.npz")),
    dict(k="hllcteno", lab="AFD–HLLC + TENO5",   c=C_HLLC, ls="--", lw=1.7,
         flux=1, rec=1,
         load=lambda: from_axframe("AX_L5_s3d_afdhllc_plt15761.npz"),
         frames=("glob", "AX_L5_s3d_afdhllc_plt*.npz")),
    dict(k="skewjst",  lab="Skew-symmetric + JST", c=C_SKEW, ls="-", lw=2.3,
         flux=None, rec=None,
         load=lambda: from_axframe("AX_L5_s3d_skewjst_plt15594.npz"),
         frames=("glob", "AX_L5_s3d_skewjst_plt*.npz")),
]

e = np.load(HERE/"panda_m142den_axis_EPAPS.npz")
EX, EE = envelope(e["x"], e["rho"])
print("%-22s %9s %9s %9s %9s" % ("source", "x_half", "E(4D)", "E(6D)", "band"))
print("%-22s %9.3f %9.3f %9.3f %9s"
      % ("experiment", x_half(EX, EE), np.interp(4.0, EX, EE),
         np.interp(6.0, EX, EE), "-"))

import glob
for c in CASES:
    c["x"], c["E"] = envelope(*c["load"]())
    c["xh"] = x_half(c["x"], c["E"])
    band = np.nan
    if c["frames"]:
        kind, pat = c["frames"]
        xs = []
        if kind == "glob":
            for f in sorted(glob.glob(str(SW/pat)))[-5:]:
                xs.append(x_half(*envelope(*from_axframe(Path(f).name))))
        else:
            z = np.load(SW/pat)
            nt = len([k for k in z.files if k.endswith("_DensityMEAN")])
            for i in range(max(0, nt-5), nt):
                xs.append(x_half(*envelope(*from_axpack(pat, i))))
        band = float(np.std(xs))
    c["band"] = band
    print("%-22s %9.3f %9.3f %9.3f %9s"
          % (c["lab"], c["xh"], np.interp(4.0, c["x"], c["E"]),
             np.interp(6.0, c["x"], c["E"]),
             "-" if band != band else "%.2f" % band))

# ---- two-by-two decomposition on the four that form the design
G = {(c["flux"], c["rec"]): c["xh"] for c in CASES if c["flux"] is not None}
fx = 0.5*((G[(1, 0)]+G[(1, 1)]) - (G[(0, 0)]+G[(0, 1)]))
rc = 0.5*((G[(0, 1)]+G[(1, 1)]) - (G[(0, 0)]+G[(1, 0)]))
ix = 0.5*((G[(1, 1)]-G[(1, 0)]) - (G[(0, 1)]-G[(0, 0)]))
print("\ntwo-by-two decomposition of x_half [D_e]")
print("  flux path   LLF -> AFD-HLLC    %+.3f" % fx)
print("  reconstruction WENO-Z5 -> TENO5 %+.3f" % rc)
print("  interaction                     %+.3f" % ix)
print("  -> the flux path moves x_half %.1f times as far as the reconstruction"
      % abs(fx/rc))

plt.rcParams.update({"font.size": 9.8, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True,
                     "legend.frameon": False})
fig = plt.figure(figsize=(12.2, 4.6))
gs = fig.add_gridspec(1, 2, width_ratios=[1.75, 1.0], wspace=0.20,
                      left=0.058, right=0.986, top=0.905, bottom=0.128)
a1 = fig.add_subplot(gs[0, 0])
a2 = fig.add_subplot(gs[0, 1])

a1.plot(EX, EE, "o", color="k", ms=3.0, mfc="none", mew=0.9, ls="", zorder=8,
        label="Panda & Seasholtz (1999)")
for c in CASES:
    a1.plot(c["x"], c["E"], color=c["c"], ls=c["ls"], lw=c["lw"], zorder=5,
            label="%s   $x_{1/2}$ = %.2f" % (c["lab"], c["xh"]))
    if c["band"] == c["band"]:
        a1.fill_betweenx([0.02, 1.4], c["xh"]-c["band"], c["xh"]+c["band"],
                         color=c["c"], alpha=0.13, lw=0, zorder=2)
a1.axhline(0.5, color="0.45", lw=0.9, ls=":", zorder=3)
a1.set_yscale("log")
a1.set_xlim(1.0, 9.0)
a1.set_ylim(0.03, 1.4)
a1.set_xlabel("$x/D_e$")
a1.set_ylabel(r"$E(x)\,/\,E(2D_e)$")
a1.set_title("(a)  shock-cell envelope of the centreline mean density,"
             " normalised at $2D_e$", fontsize=9.8, pad=5)
a1.legend(loc="lower left", fontsize=8.3, handlelength=2.2)

for c in CASES:
    if c["flux"] is None:
        continue
    a2.plot(c["flux"], c["xh"], "o", ms=9, mfc=c["c"], mec="k", mew=0.8,
            zorder=6, alpha=1.0 if c["rec"] == 0 else 0.5)
for rec, ls, lab in ((0, "-", "WENO-Z5"), (1, "--", "TENO5")):
    y = [G[(0, rec)], G[(1, rec)]]
    a2.plot([0, 1], y, ls=ls, color="0.35", lw=1.5, zorder=4, label=lab)
sk = [c for c in CASES if c["flux"] is None][0]
a2.plot([-0.34], [sk["xh"]], "s", ms=9, mfc=C_SKEW, mec="k", mew=0.8, zorder=6)
a2.annotate("skew+JST", (-0.34, sk["xh"]), textcoords="offset points",
            xytext=(0, -17), fontsize=8.4, ha="center")
a2.axhline(x_half(EX, EE), color="k", lw=1.0, ls=":", zorder=3)
a2.annotate("experiment", (-0.52, x_half(EX, EE)),
            textcoords="offset points", xytext=(2, 5), fontsize=8.4,
            ha="left", color="0.25")
a2.set_xlim(-0.55, 1.30)
a2.set_xticks([0, 1])
a2.set_xticklabels(["LLF", "AFD–HLLC"])
a2.set_xlabel("flux path")
a2.set_ylabel(r"$x_{1/2}/D_e$")
a2.set_title("(b)  the two factors separated", fontsize=9.8, pad=5)
a2.legend(loc="upper left", fontsize=8.6, handlelength=2.4, title="reconstruction",
          title_fontsize=8.4)
a2.text(0.97, 0.045,
        "flux path  %+.2f $D_e$\nreconstruction  %+.2f\ninteraction  %+.2f"
        % (fx, rc, ix), transform=a2.transAxes, ha="right", va="bottom",
        fontsize=8.4, color="0.25", linespacing=1.5)

for ax in (a1, a2):
    ax.grid(True, ls="--", lw=0.5, color="0.93")
    ax.set_axisbelow(True)

fig.savefig(HERE/"dissipation_envelope.png", dpi=250)
fig.savefig(HERE/"dissipation_envelope.pdf")
print("\nwrote dissipation_envelope.png/pdf")
