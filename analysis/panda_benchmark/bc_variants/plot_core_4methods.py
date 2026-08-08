#!/usr/bin/env python3
"""Supersonic core of the four scheme combinations.

Rows are ordered LLF+WENO-Z5, LLF+TENO5, AFD-HLLC+WENO-Z5, AFD-HLLC+TENO5, so
the two flux paths sit as adjacent pairs and the reconstruction alternates
inside each pair. One colour bar serves all four rows, so a difference between
rows is a difference in the data and never in the mapping.

The field plotted is the Mach number of the mean field,

    M_bar = |u_bar| / sqrt(gamma p_bar / rho_bar),

which is what the stored statistics support: the code accumulates means of the
primitive variables, not a mean of the pointwise Mach number, and the two are
not the same quantity in a compressible shear layer. The sonic line is the
M_bar = 1 contour of that field and the supersonic core is the region it
encloses.

On the axis the core is not simply connected. The mean Mach number falls below
one immediately behind the first recompression - the incipient Mach disk at
this pressure ratio, eta_0 = 3.27, which sits between the Muraoka-Hiejima
first-disk band and the Gibbings regular-to-Mach transition - and recovers
within a fraction of a diameter. Everywhere else along 0.05 < x/D_e < 7.2 the
axis is supersonic in all four schemes, so the core does not terminate inside
the slice and no core length can be quoted from these data. What the four
schemes do differ in is that subsonic pocket, which the second figure resolves.

    core_4methods_mach.png       M_bar with the sonic line, four schemes
    core_4methods_centreline.png centreline M_bar and the subsonic pocket
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
F = HERE/"fields"
GAM = 1.4
YMAX = 1.20
XLO, XHI = 0.05, 7.20

TAGS = ["llfweno", "llfteno", "hllcweno", "hllcteno"]
ROWS = ["LLF + WENO-Z5", "LLF + TENO5", "AFD-HLLC + WENO-Z5", "AFD-HLLC + TENO5"]
# flux path carries the colour, reconstruction the line style, so the 2x2
# design of the campaign is readable straight off the legend
STY = dict(llfweno=("#0072B2", "-"), llfteno=("#0072B2", "--"),
           hllcweno=("#D55E00", "-"), hllcteno=("#D55E00", "--"))


def mach(tag):
    def ld(v):
        z = np.load(F/("SLM_%s_%s.npz" % (tag, v)))
        return z["field"].astype(np.float64), z
    rm, z = ld("DensityMEAN")
    pm, _ = ld("pressureMEAN")
    um, _ = ld("x_velocityMEAN")
    vm, _ = ld("y_velocityMEAN")
    wm, _ = ld("z_velocityMEAN")
    M = np.sqrt(um*um + vm*vm + wm*wm)/np.sqrt(GAM*pm/rm)
    return dict(M=M, x=z["x"], y=z["y"], t=float(z["time"])*1e3)


def pockets(x, a):
    """Runs of the axis where the mean Mach number is below one."""
    w = (x >= XLO) & (x <= XHI)
    xx, aa = x[w], a[w]
    s = aa < 1.0
    e = np.r_[0, np.flatnonzero(np.diff(s.astype(int)) != 0)+1, s.size]
    return [(xx[b], xx[c-1], aa[b:c].min())
            for b, c in zip(e[:-1], e[1:]) if s[b]], xx, aa


Q = {}
print("%-18s %7s %7s %7s %7s %8s %8s %8s"
      % ("scheme", "M_max", "M_ax_7", "n_pkt", "L_pkt", "M_min", "x_start", "x_end"))
for t in TAGS:
    q = mach(t)
    j0 = int(np.argmin(np.abs(q["y"])))
    p, xx, aa = pockets(q["x"], q["M"][:, j0])
    q.update(ax=aa, xax=xx, pkt=p, j0=j0)
    Q[t] = q
    L = sum(b-a for a, b, _ in p)
    print("%-18s %7.3f %7.3f %7d %7.4f %8.4f %8.4f %8.4f"
          % (ROWS[TAGS.index(t)], q["M"].max(),
             aa[np.argmin(np.abs(xx-7.0))], len(p), L,
             min(m for _, _, m in p), p[0][0], p[-1][1]))

plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True,
                     "legend.frameon": False})

# --------------------------------------------- mean Mach field + sonic line
fig, AX = plt.subplots(4, 1, figsize=(9.4, 8.4), sharex=True, sharey=True)
fig.subplots_adjust(left=0.078, right=0.872, top=0.988, bottom=0.058,
                    hspace=0.085)
nrm = Normalize(0.0, 2.40)
for ax, tag, name in zip(AX, TAGS, ROWS):
    q = Q[tag]
    m = np.abs(q["y"]) <= YMAX
    im = ax.pcolormesh(q["x"], q["y"][m], q["M"].T[m], cmap="viridis",
                       norm=nrm, shading="nearest", rasterized=True)
    ax.contour(q["x"], q["y"][m], q["M"].T[m], levels=[1.0], colors="w",
               linewidths=1.0, zorder=4)
    ax.set_xlim(0, 7.2)
    ax.set_ylim(-YMAX, YMAX)
    ax.set_aspect("equal")
    ax.set_ylabel("$y/D_e$", labelpad=1)
    ax.text(0.009, 0.86, name, transform=ax.transAxes, color="w", fontsize=10)
AX[-1].set_xlabel("$x/D_e$", labelpad=2)
AX[0].plot([], [], "-", color="w", lw=1.0, label=r"sonic line, $\overline{M}=1$")
AX[0].legend(loc="lower left", fontsize=8.6, labelcolor="w",
             bbox_to_anchor=(0.005, 0.02))
cax = fig.add_axes([0.888, 0.058, 0.018, 0.930])
cb = fig.colorbar(im, cax=cax)
cb.set_label(r"$\overline{M}$", fontsize=11.5, labelpad=4)
cb.ax.tick_params(labelsize=9)
fig.savefig(HERE/"core_4methods_mach.png", dpi=250)
fig.savefig(HERE/"core_4methods_mach.pdf")
plt.close(fig)

# ------------------------------------------------ centreline and the pocket
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.4, 4.2),
                             gridspec_kw=dict(width_ratios=[1.95, 1.0]))
fig.subplots_adjust(left=0.060, right=0.988, top=0.925, bottom=0.135,
                    wspace=0.185)
for tag, name in zip(TAGS, ROWS):
    q = Q[tag]
    c, ls = STY[tag]
    L = sum(b-a for a, b, _ in q["pkt"])
    a1.plot(q["xax"], q["ax"], color=c, ls=ls, lw=1.4, zorder=5,
            label="%s   $L_{sub}$ = %.3f $D_e$" % (name, L))
    a2.plot(q["xax"], q["ax"], color=c, ls=ls, lw=1.7, zorder=5)
for ax in (a1, a2):
    ax.axhline(1.0, color="0.35", lw=1.0, ls=":", zorder=3)
    ax.set_xlabel("$x/D_e$")
    ax.grid(True, ls="--", lw=0.5, color="0.93")
    ax.set_axisbelow(True)
a1.set_xlim(0, 7.2)
a1.set_ylim(0.7, 2.5)
a1.set_ylabel(r"$\overline{M}$   on the axis")
a1.set_title(r"(a)  Mach number of the mean field on the axis",
             fontsize=10, pad=5)
a1.legend(loc="upper right", fontsize=8.8, handlelength=2.4)
a2.set_xlim(0.85, 1.75)
a2.set_ylim(0.82, 1.35)
a2.set_ylabel(r"$\overline{M}$")
a2.set_title("(b)  the subsonic pocket, first recompression",
             fontsize=10, pad=5)
a2.axhspan(0.82, 1.0, color="0.86", alpha=0.55, zorder=1, lw=0)
a2.text(0.030, 0.055, r"$\overline{M}<1$", transform=a2.transAxes,
        fontsize=9.5, color="0.30")
fig.savefig(HERE/"core_4methods_centreline.png", dpi=250)
fig.savefig(HERE/"core_4methods_centreline.pdf")
print("\nwrote core_4methods_mach.png/pdf and core_4methods_centreline.png/pdf")
