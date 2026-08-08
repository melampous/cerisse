#!/usr/bin/env python3
"""NPR 20 top hat at Level 4, flange against no flange.

One literal separates the two runs: r_flange goes from 0.1524 m (6 D_e, the
Panda plate) to r_jet, which empties the second radial zone of the z-low
boundary, so the adiabatic slip wall between the lip and 6 D_e is replaced by
the quiescent ambient Dirichlet and the jet entrains from a free boundary.
Grid, scheme, CFL, statistics window and refinement are byte-identical and the
binary was rebuilt rather than cloned (md5 0d18a41ecbd1 against f689f3485f92),
so the difference in the two fields is the boundary and nothing else.

Two figures.

n20_nofl_schlieren.png    numerical schlieren of the mean density, full view
                          and the Mach-disk region enlarged,
                              S = exp(-k |grad rho| / g_ref)
                          with g_ref the 99.5th percentile of |grad rho| over
                          x > 0.8 D_e, taken separately for each case. The
                          reference is a percentile and not the maximum because
                          the top-hat lip is a discontinuity in the imposed
                          profile and its gradient is about seventy times
                          anything downstream.

n20_nofl_centreline.png   centreline mean density with x_fc marked, and the
                          disk region enlarged.

x_fc is the steepest positive axial gradient of the centreline mean density
between the first expansion minimum and the following maximum. It is computed
here on the *native* grid: the extraction writes every level onto the finest
abscissa, so a level-2 cell appears as four identical output cells, and a
centred difference over that staircase produces a tied gradient pair whose
parabolic peak lands exactly half a cell downstream regardless of the data.
That artefact quantises x_fc to a cell edge and made the two cases return
byte-identical values. Collapsing each run of replicated values to its own
centre first removes it.

The refinement patches were sized for the NPR 3 condition, where the disk sits
near x = 1 D_e. At NPR 20 the disk stands at about 3 D_e, outside the finest
patch, and the local cell there is 397 um - level 2, four times the nominal
99.2 um. The figures state this rather than implying a Level-4 disk.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, ConnectionPatch

HERE = Path(__file__).resolve().parent
F, OUT = HERE/"fields2d", HERE/"sweep_out"
RHO_J, D = 1.6413, 0.0254
K = 5.0
ZW, ZH = 0.62, 0.78

CASES = [("with flange", "XFC_N20_tophat_LO__L4.npz",
          "FLD_N20_DensityMEAN.npz", "#0072B2"),
         ("no flange",   "XFC_N20_tophat_LO_nofl__L4.npz",
          "FLD_N20nofl_DensityMEAN.npz", "#D55E00")]


def centreline(f):
    d = np.load(OUT/f)
    ks = sorted([k[:-2] for k in d.files if k.endswith("_t")],
                key=lambda q: int("".join(c for c in q if c.isdigit())))
    dx = float(d[ks[0]+"_t"][1])
    r = d[ks[-1]+"_rho"]/RHO_J
    x = (np.arange(r.size)+0.5)*dx/D
    return x, r, float(np.atleast_1d(d[ks[-1]+"_t"])[0]), len(ks)


def native(x, r):
    """Collapse each run of replicated output cells onto its own centre, so the
    abscissa carries one point per cell of the level that actually owns it."""
    m = np.isfinite(r)
    x, r = x[m], r[m]
    e = np.r_[0, np.flatnonzero(np.diff(r) != 0)+1, r.size]
    xc = np.array([0.5*(x[a]+x[b-1]) for a, b in zip(e[:-1], e[1:])])
    return xc, r[e[:-1]]


def xfc_of(x, r, lo=1.70, hi=3.80):
    """Steepest positive gradient of the recompression, on the native grid.

    The difference is one-sided and evaluated at cell midpoints, which is the
    only consistent operator on a grid whose spacing jumps at every refinement
    boundary; a centred difference would mix two spacings into one estimate.
    """
    m = (x >= lo) & (x <= hi)
    xx, rr = x[m], r[m]
    i0 = int(np.argmin(rr))
    i1 = i0 + int(np.argmax(rr[i0:]))
    xm = 0.5*(xx[:-1]+xx[1:])
    g = np.diff(rr)/np.diff(xx)
    s = (xm > xx[i0]) & (xm < xx[i1])
    j = int(np.flatnonzero(s)[int(np.nanargmax(g[s]))])
    xr = xm[j]
    if 0 < j < len(g)-1:                       # parabolic peak, guarded
        d2 = g[j-1]-2*g[j]+g[j+1]
        if d2 < 0:
            h = 0.5*(g[j-1]-g[j+1])/d2
            if abs(h) <= 0.5:
                xr = xm[j] + h*0.5*(xm[j+1]-xm[j-1])
    return dict(xfc=xr, grad=g[j], hloc=xx[j+1]-xx[j],
                xmin=xx[i0], rmin=rr[i0], xmax=xx[i1], rmax=rr[i1])


def grad_rf(a, h, axis):
    """Derivative of a composited AMR field along one axis, differencing over
    the cell that actually owns the data.

    The composite writes every level onto the finest mesh, so a level-2 cell
    appears as four identical entries. Differencing that with the finest
    spacing gives zero inside each block and the whole jump across its edge,
    which is the comb-shaped texture that otherwise covers every coarse patch
    and swamps the shock it is meant to show. Here each line is run-length
    encoded, the derivative is taken between the centres of neighbouring runs,
    and the result is broadcast back over the run - so the operator uses the
    local cell size everywhere, whatever level owns it.

    Bit-identical neighbours are the signature of replication because the
    composite copies values; two genuinely equal cells on the finest level are
    merged as well, but only where the field is flat and the derivative is
    zero either way.
    """
    a = a if axis == 1 else a.T
    out = np.zeros_like(a)
    n = a.shape[1]
    idx = np.arange(n)
    for k in range(a.shape[0]):
        row = a[k]
        e = np.r_[0, np.flatnonzero(np.diff(row) != 0)+1, n]
        if e.size < 3:
            continue
        c = 0.5*(e[:-1]+e[1:]-1)*h              # run centres, physical
        v = row[e[:-1]]
        g = np.gradient(v, c) if v.size > 2 else np.diff(v)/np.diff(c)
        out[k] = np.repeat(g if v.size > 2 else np.r_[g, g[-1]], np.diff(e))
    return out if axis == 1 else out.T


def schlieren(f):
    z = np.load(F/f)
    rho = z["F"].astype(np.float64)
    dr = float(z["dr"])
    r, x = z["r"]/D, z["z"]/D
    gz = grad_rf(rho, dr, 1)
    gr = grad_rf(rho, dr, 0)
    g = np.sqrt(gr*gr + gz*gz)
    gref = np.nanpercentile(g[:, x > 0.8], 99.5)
    s = np.exp(-K*np.clip(g/gref, 0.0, 1.0))
    return x, np.r_[-r[::-1], r], np.vstack([s[::-1, :], s]), \
        float(z["time"])*1e3, gref


Q = {}
print("%-12s %8s %7s %9s %9s %9s %8s %8s"
      % ("case", "t[ms]", "frames", "x_fc/De", "h_loc/De", "grad", "x_min", "rho_min"))
for lab, cf, ff, col in CASES:
    x, r, t, n = centreline(cf)
    xn, rn = native(x, r)
    q = xfc_of(xn, rn)
    q.update(lab=lab, x=xn, r=rn, t=t, n=n, col=col, fld=ff)
    Q[lab] = q
    print("%-12s %8.3f %7d %9.4f %9.4f %9.2f %8.4f %8.4f"
          % (lab, t*1e3, n, q["xfc"], q["hloc"], q["grad"], q["xmin"], q["rmin"]))

A, B = Q["with flange"], Q["no flange"]
dv = B["xfc"]-A["xfc"]
CR, AS = 0.6455*np.sqrt(20.0), 0.67*np.sqrt(20.0)
print("\n  x_fc/D_e   flange %.4f   no flange %.4f   difference %+.4f D_e (%+.2f %%)"
      % (A["xfc"], B["xfc"], dv, 100*dv/A["xfc"]))
print("  local cell at the disk %.4f D_e (%.0f um, level 2), so the difference is"
      " %.2f cells and is not resolved" % (A["hloc"], A["hloc"]*D*1e6,
                                           abs(dv)/A["hloc"]))
print("  /Crist %.4f and %.4f   /Ashkenas-Sherman %.4f and %.4f"
      % (A["xfc"]/CR, B["xfc"]/CR, A["xfc"]/AS, B["xfc"]/AS))
print("  first expansion minimum  rho/rho_j %.4f at %.4f  ->  %.4f at %.4f"
      "   (%+.4f D_e = %.1f cells)"
      % (A["rmin"], A["xmin"], B["rmin"], B["xmin"],
         B["xmin"]-A["xmin"], (B["xmin"]-A["xmin"])/A["hloc"]))

plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True,
                     "legend.frameon": False})

# ---------------------------------------------------------------- schlieren
fig, AX = plt.subplots(2, 2, figsize=(11.6, 5.9),
                       gridspec_kw=dict(width_ratios=[2.45, 1.0]))
fig.subplots_adjust(left=0.055, right=0.988, top=0.955, bottom=0.088,
                    wspace=0.105, hspace=0.225)
print("\n%-12s %12s" % ("case", "g_ref(p99.5)"))
for row, lab in enumerate(("with flange", "no flange")):
    q = Q[lab]
    x, rr, ss, t, gref = schlieren(q["fld"])
    print("%-12s %12.1f" % (lab, gref))
    xf = q["xfc"]
    a, b = AX[row, 0], AX[row, 1]
    for ax in (a, b):
        ax.pcolormesh(x, rr, ss, cmap="gray", shading="nearest",
                      rasterized=True, vmin=0.0, vmax=1.0)
        ax.axvline(xf, color="#e8482c", lw=1.0, ls="--", zorder=5)
        ax.set_aspect("equal")
    a.set_xlim(0, 6.0)
    a.set_ylim(-1.6, 1.6)
    a.set_ylabel("$r/D_e$")
    a.text(0.010, 0.90, "NPR$_0$ = 20, %s,   $t$ = %.2f ms" % (lab, t),
           transform=a.transAxes, fontsize=10,
           bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="none",
                     alpha=0.86))
    a.add_patch(Rectangle((xf-ZW, -ZH), 2*ZW, 2*ZH, fill=False,
                          edgecolor="#e8482c", lw=1.1, zorder=6))
    b.set_xlim(xf-ZW, xf+ZW)
    b.set_ylim(-ZH, ZH)
    b.set_ylabel("$r/D_e$", labelpad=1)
    bb = dict(boxstyle="round,pad=0.22", fc="white", ec="none", alpha=0.86)
    b.text(0.035, 0.905, r"$x_{fc}$ = %.4f" % xf, transform=b.transAxes,
           fontsize=9.5, color="#e8482c", va="top", bbox=bb)
    b.text(0.035, 0.055, "local cell %.0f $\\mu$m (level 2)" % (q["hloc"]*D*1e6),
           transform=b.transAxes, fontsize=8.2, color="0.25", bbox=bb)
    for sp in b.spines.values():
        sp.set_color("#e8482c")
        sp.set_linewidth(1.1)
    for yy in (-ZH, ZH):
        fig.add_artist(ConnectionPatch(
            xyA=(xf+ZW, yy), coordsA=a.transData,
            xyB=(xf-ZW, yy), coordsB=b.transData,
            color="#e8482c", lw=0.7, ls=":", zorder=1))
AX[1, 0].set_xlabel("$x/D_e$")
AX[1, 1].set_xlabel("$x/D_e$")
fig.savefig(HERE/"n20_nofl_schlieren.png", dpi=300)
fig.savefig(HERE/"n20_nofl_schlieren.pdf")

# --------------------------------------------------------------- centreline
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.4, 4.2),
                             gridspec_kw=dict(width_ratios=[1.9, 1.0]))
fig.subplots_adjust(left=0.062, right=0.988, top=0.925, bottom=0.135,
                    wspace=0.20)
for lab in ("with flange", "no flange"):
    q = Q[lab]
    a1.plot(q["x"], q["r"], color=q["col"], lw=1.5, zorder=5,
            label="%s,  $x_{fc}$ = %.4f" % (lab, q["xfc"]))
    a2.plot(q["x"], q["r"], color=q["col"], lw=1.6, zorder=5)
    for ax in (a1, a2):
        ax.axvline(q["xfc"], color=q["col"], lw=1.0, ls="--", zorder=3)
    a2.plot([q["xmin"]], [q["rmin"]], "v", color=q["col"], ms=7, mec="w",
            mew=0.7, zorder=7)
a1.set_xlim(0, 7.0)
a1.set_ylim(0, 5.0)
a1.set_ylabel(r"$\overline{\rho}/\rho_j$   on the axis")
a1.set_title(r"(a)  NPR$_0$ = 20 top hat, Level 4, centreline mean density",
             fontsize=10, pad=5)
a1.legend(loc="upper right", fontsize=9.0, handlelength=1.9)
lo = min(Q[k]["xfc"] for k in Q) - 0.7
hi = max(Q[k]["xfc"] for k in Q) + 0.7
a2.set_xlim(lo, hi)
a2.set_ylim(0.15, 1.20)
a2.set_ylabel(r"$\overline{\rho}/\rho_j$")
a2.set_title("(b)  the Mach-disk region", fontsize=10, pad=5)
a2.text(0.035, 0.93,
        "$\\Delta x_{fc}$ = %+.4f $D_e$\nlocal cell %.4f $D_e$\n= %.2f cells,"
        " not resolved" % (dv, A["hloc"], abs(dv)/A["hloc"]),
        transform=a2.transAxes, va="top", fontsize=9, color="0.25",
        linespacing=1.5)
for ax in (a1, a2):
    ax.set_xlabel("$x/D_e$")
    ax.grid(True, ls="--", lw=0.5, color="0.93")
    ax.set_axisbelow(True)
fig.savefig(HERE/"n20_nofl_centreline.png", dpi=250)
fig.savefig(HERE/"n20_nofl_centreline.pdf")
print("\nwrote n20_nofl_schlieren.png/pdf and n20_nofl_centreline.png/pdf")
