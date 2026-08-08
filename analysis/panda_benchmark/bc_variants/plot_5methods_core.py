#!/usr/bin/env python3
"""Supersonic core and numerical schlieren across five scheme combinations.

Rows are LLF+WENO-Z5, LLF+TENO5, AFD-HLLC+WENO-Z5, AFD-HLLC+TENO5 and
skew-symmetric+JST. The first four come from the SLM_ extraction, which kept
nine variables including all three velocity components and three planes in k;
the skew-JST case was extracted earlier with a smaller variable set - mean
density, mean pressure and the axial mean velocity only.

Consequences, stated rather than worked around:

  * The Mach number of the mean field is formed from the AXIAL component alone
    for every row, M_bar = |u_bar| / sqrt(gamma p_bar / rho_bar), so all five
    rows are the same construction. Against the full three-component value,
    available for the first four, the axial-only form differs by at most 0.006
    on the axis and 0.051 anywhere inside |y| < 0.8 D_e; the 99th percentile of
    the difference is 0.023 to 0.027. On the axis this is negligible. Off axis
    the sonic line can shift slightly where the transverse mean velocity is not
    small, which is the shear layer.
  * Vorticity, Q and lambda_2 are not computed here. They need the transverse
    velocity components and the k planes, which the skew-JST extraction does
    not carry; those figures stay at four schemes.

The two extractions share the cell size and the y grid exactly. The skew-JST
slice ends at x = 6.998 D_e against 7.201, so every panel is cut at 6.99 and
all integral measures use the common window.

The schlieren gradient is taken with the same run-length operator on all five,
so no row benefits from the rfmap that only the SLM_ files carry: each line is
run-length encoded, the derivative is taken between the centres of neighbouring
runs and broadcast back. Differencing the composite with the finest spacing
instead would return zero inside every replicated coarse block and the whole
jump at its edge, which is the comb texture that otherwise covers each coarse
patch.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

HERE = Path(__file__).resolve().parent
FLD, SWP = HERE/"fields", HERE/"sweep_out"
GAM, D = 1.4, 0.0254
YMAX = 1.20
XLO, XHI = 0.05, 6.99

CASES = [("llfweno",  "LLF + WENO-Z5",       FLD, "SLM_llfweno_%s.npz",  "#0072B2", "-"),
         ("llfteno",  "LLF + TENO5",         FLD, "SLM_llfteno_%s.npz",  "#0072B2", "--"),
         ("hllcweno", "AFD-HLLC + WENO-Z5",  FLD, "SLM_hllcweno_%s.npz", "#D55E00", "-"),
         ("hllcteno", "AFD-HLLC + TENO5",    FLD, "SLM_hllcteno_%s.npz", "#D55E00", "--"),
         ("skewjst",  "skew-symmetric + JST", SWP, "SL_L5_s3d_skewjst_%s.npz",
          "#009E73", "-")]


def grad_rf(a, h, axis):
    """Derivative along one axis of a composited AMR field, differencing over
    the cell that actually owns the data. Bit-identical neighbours mark a
    replicated coarse cell because the composite copies values."""
    a = a if axis == 1 else a.T
    out = np.zeros_like(a)
    n = a.shape[1]
    for k in range(a.shape[0]):
        row = a[k]
        e = np.r_[0, np.flatnonzero(np.diff(row) != 0)+1, n]
        if e.size < 3:
            continue
        c = 0.5*(e[:-1]+e[1:]-1)*h
        v = row[e[:-1]]
        g = np.gradient(v, c) if v.size > 2 else np.diff(v)/np.diff(c)
        out[k] = np.repeat(g if v.size > 2 else np.r_[g, g[-1]], np.diff(e))
    return out if axis == 1 else out.T


def load(dirpath, pat):
    def one(v):
        z = np.load(dirpath/(pat % v))
        return z["field"].astype(np.float64), z
    rm, z = one("DensityMEAN")
    pm, _ = one("pressureMEAN")
    um, _ = one("x_velocityMEAN")
    d = float(z["dx_um"])*1e-6
    M = np.abs(um)/np.sqrt(GAM*pm/rm)
    g = np.hypot(grad_rf(rm, d, 0), grad_rf(rm, d, 1))
    return dict(M=M, grho=g, x=z["x"], y=z["y"], t=float(z["time"])*1e3)


Q = {}
for tag, name, dp, pat, col, ls in CASES:
    q = load(dp, pat)
    q.update(name=name, col=col, ls=ls)
    j0 = int(np.argmin(np.abs(q["y"])))
    w = (q["x"] >= XLO) & (q["x"] <= XHI)
    xx, aa = q["x"][w], q["M"][w, j0]
    s = aa < 1.0
    e = np.r_[0, np.flatnonzero(np.diff(s.astype(int)) != 0)+1, s.size]
    q["pkt"] = [(xx[b], xx[c-1], aa[b:c].min())
                for b, c in zip(e[:-1], e[1:]) if s[b]]
    q["xax"], q["ax"] = xx, aa
    Q[tag] = q

print("common window %.2f < x/D_e < %.2f\n" % (XLO, XHI))
print("%-22s %7s %8s %7s %8s %8s %9s"
      % ("scheme", "M_max", "M_ax(6.9)", "n_pkt", "L_sub", "M_min", "y_sonic"))
for tag, name, *_ in CASES:
    q = Q[tag]
    L = sum(b-a for a, b, _ in q["pkt"])
    mm = min((m for _, _, m in q["pkt"]), default=np.nan)
    ys = np.nanmax(np.abs(q["y"][np.any(q["M"] >= 1.0, axis=0)]))
    print("%-22s %7.3f %8.3f %7d %8.4f %8.4f %9.3f"
          % (name, q["M"].max(), q["ax"][np.argmin(np.abs(q["xax"]-6.9))],
             len(q["pkt"]), L, mm, ys))

# shared schlieren reference across all five rows
pool = np.concatenate([Q[t]["grho"][Q[t]["x"] > 0.8].ravel() for t, *_ in CASES])
GREF = np.nanpercentile(pool, 99.5)
for t, *_ in CASES:
    Q[t]["sch"] = np.exp(-5.0*np.clip(Q[t]["grho"]/GREF, 0.0, 1.0))
print("\nshared schlieren reference g_ref = %.1f" % GREF)

plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True,
                     "legend.frameon": False})


def stack5(key, cmap, lo, hi, cbar, out, sonic=False, light=True):
    fig, AX = plt.subplots(5, 1, figsize=(9.4, 10.2), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.078, right=0.872, top=0.990, bottom=0.048,
                        hspace=0.085)
    nrm = Normalize(lo, hi)
    for ax, (tag, name, *_) in zip(AX, CASES):
        q = Q[tag]
        m = np.abs(q["y"]) <= YMAX
        im = ax.pcolormesh(q["x"], q["y"][m], q[key].T[m], cmap=cmap, norm=nrm,
                           shading="nearest", rasterized=True)
        if sonic:
            ax.contour(q["x"], q["y"][m], q["M"].T[m], levels=[1.0],
                       colors="w", linewidths=1.0, zorder=4)
        ax.set_xlim(0, XHI)
        ax.set_ylim(-YMAX, YMAX)
        ax.set_aspect("equal")
        ax.set_ylabel("$y/D_e$", labelpad=1)
        ax.text(0.009, 0.86, name, transform=ax.transAxes,
                color="w" if light else "k", fontsize=10)
    AX[-1].set_xlabel("$x/D_e$", labelpad=2)
    if sonic:
        AX[0].plot([], [], "-", color="w", lw=1.0,
                   label=r"sonic line, $\overline{M}=1$")
        AX[0].legend(loc="lower left", fontsize=8.6, labelcolor="w",
                     bbox_to_anchor=(0.005, 0.02))
    cax = fig.add_axes([0.888, 0.048, 0.018, 0.942])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(cbar, fontsize=11.5, labelpad=4)
    cb.ax.tick_params(labelsize=9)
    fig.savefig(HERE/(out+".png"), dpi=250)
    fig.savefig(HERE/(out+".pdf"))
    plt.close(fig)
    print("wrote %s.png/pdf" % out)


stack5("M", "viridis", 0.0, 2.40, r"$\overline{M}$", "core5_mach", sonic=True)
stack5("sch", "gray", 0.0, 1.0, "numerical schlieren of "
       r"$\overline{\rho}$,   $\exp(-k|\nabla\overline{\rho}|/g_{ref})$",
       "core5_schlieren", light=False)

# ------------------------------------------------------------- centreline
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.4, 4.2),
                             gridspec_kw=dict(width_ratios=[1.95, 1.0]))
fig.subplots_adjust(left=0.060, right=0.988, top=0.925, bottom=0.135,
                    wspace=0.185)
for tag, name, *_ in CASES:
    q = Q[tag]
    L = sum(b-a for a, b, _ in q["pkt"])
    a1.plot(q["xax"], q["ax"], color=q["col"], ls=q["ls"], lw=1.4, zorder=5,
            label="%s   $L_{sub}$ = %.3f $D_e$" % (name, L))
    a2.plot(q["xax"], q["ax"], color=q["col"], ls=q["ls"], lw=1.7, zorder=5)
for ax in (a1, a2):
    ax.axhline(1.0, color="0.35", lw=1.0, ls=":", zorder=3)
    ax.set_xlabel("$x/D_e$")
    ax.grid(True, ls="--", lw=0.5, color="0.93")
    ax.set_axisbelow(True)
a1.set_xlim(0, XHI)
a1.set_ylim(0.7, 2.5)
a1.set_ylabel(r"$\overline{M}$   on the axis")
a1.set_title(r"(a)  Mach number of the mean field on the axis, five schemes",
             fontsize=10, pad=5)
a1.legend(loc="upper right", fontsize=8.4, handlelength=2.4)
a2.set_xlim(0.85, 1.75)
a2.set_ylim(0.82, 1.35)
a2.set_ylabel(r"$\overline{M}$")
a2.set_title("(b)  the subsonic pocket, first recompression", fontsize=10, pad=5)
a2.axhspan(0.82, 1.0, color="0.86", alpha=0.55, zorder=1, lw=0)
a2.text(0.030, 0.055, r"$\overline{M}<1$", transform=a2.transAxes,
        fontsize=9.5, color="0.30")
fig.savefig(HERE/"core5_centreline.png", dpi=250)
fig.savefig(HERE/"core5_centreline.pdf")
print("wrote core5_centreline.png/pdf")
