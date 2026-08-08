#!/usr/bin/env python3
"""The subsonic pocket behind the first recompression, five formulations.

Three columns of the same region, one row per formulation:

    mean Mach number      M_bar = |u_bar| / sqrt(gamma p_bar / rho_bar), with
                          the sonic line M_bar = 1 drawn in white and the
                          extent of the subsonic pocket on the axis marked
    vorticity magnitude   |omega| D_e / U_j, instantaneous
    numerical schlieren   exp(-k |grad rho| / g_ref), instantaneous

This is the only place in the flow where the five formulations disagree about
the topology rather than the amplitude: the axis is supersonic everywhere else
in the slice, and the pocket is where the incipient Mach disk at eta_0 = 3.27
sits. Its length spans a factor of 2.4 across the set.

All three quantities now use the three velocity components for every row: the
skew-symmetric case was re-extracted with the same nine variables as the other
four plus instantaneous density and pressure, so the axial-only Mach
approximation used before is no longer needed and no row is computed
differently from any other.

One colour bar per column serves all five rows, so a difference down a column
is a difference in the data. The schlieren reference is shared across the five
as well, taken as the 99.5th percentile of |grad rho| over x > 0.8 D_e pooled
over the set; a per-row reference would hide exactly the sharpness difference
the column is for.

Derivatives use the run-length-aware operator: the slice is a finest-wins AMR
composite, so a plain difference returns zero inside every replicated coarse
block and the whole jump at its edge.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

HERE = Path(__file__).resolve().parent
F = HERE/"fields"
RHO_J, U_J, D, GAM = 1.6413, 414.2, 0.0254, 1.4
ZX = (0.80, 1.80)
ZY = (-0.65, 0.65)

TAGS = ["llfweno", "llfteno", "hllcweno", "hllcteno", "skewjst"]
ROWS = ["LLF + WENO-Z5", "LLF + TENO5", "AFD-HLLC + WENO-Z5",
        "AFD-HLLC + TENO5", "Skew-symmetric + JST"]


def load(tag, var, stack=False):
    z = np.load(F/("SLM_%s_%s.npz" % (tag, var)))
    return (z["kstack"] if stack else z["field"]).astype(np.float64), z


def ddx(f, rf, d, axis):
    out = np.zeros_like(f)
    for r in np.unique(rf[rf > 0]):
        r = int(r)
        g = (np.roll(f, -r, axis=axis) - np.roll(f, r, axis=axis))/(2.0*r*d)
        sl = [slice(None)]*f.ndim
        sl[axis] = slice(0, r)
        g[tuple(sl)] = 0.0
        sl[axis] = slice(-r, None)
        g[tuple(sl)] = 0.0
        out = np.where(rf == r, g, out)
    return out


def fields(tag):
    rm, z = load(tag, "DensityMEAN")
    pm, _ = load(tag, "pressureMEAN")
    um, _ = load(tag, "x_velocityMEAN")
    vm, _ = load(tag, "y_velocityMEAN")
    wm, _ = load(tag, "z_velocityMEAN")
    ri, _ = load(tag, "Density")
    rf, d = z["rfmap"], float(z["dx_um"])*1e-6

    U, _ = load(tag, "x_velocity", True)
    V, _ = load(tag, "y_velocity", True)
    W, _ = load(tag, "z_velocity", True)
    dz = 2.0*rf*d
    wx = ddx(W[1], rf, d, 1) - (V[2]-V[0])/dz
    wy = (U[2]-U[0])/dz - ddx(W[1], rf, d, 0)
    wz = ddx(V[1], rf, d, 0) - ddx(U[1], rf, d, 1)

    M = np.sqrt(um*um + vm*vm + wm*wm)/np.sqrt(GAM*pm/rm)
    x, y = z["x"], z["y"]
    j0 = int(np.argmin(np.abs(y)))
    a = M[:, j0]
    w = (x >= 0.6) & (x <= 2.2)
    xx, aa = x[w], a[w]
    s = aa < 1.0
    e = np.r_[0, np.flatnonzero(np.diff(s.astype(int)) != 0)+1, s.size]
    pkt = [(xx[b], xx[c-1]) for b, c in zip(e[:-1], e[1:]) if s[b]]

    return dict(M=M, omg=np.sqrt(wx*wx + wy*wy + wz*wz)*D/U_J,
                grho=np.hypot(ddx(ri, rf, d, 0), ddx(ri, rf, d, 1)),
                x=x, y=y, pkt=pkt, t=float(z["time"])*1e3)


Q = {t: fields(t) for t in TAGS}
pool = np.concatenate([Q[t]["grho"][Q[t]["x"] > 0.8].ravel() for t in TAGS])
GREF = np.nanpercentile(pool, 99.5)
for t in TAGS:
    Q[t]["sch"] = np.exp(-5.0*np.clip(Q[t]["grho"]/GREF, 0.0, 1.0))

print("subsonic pocket on the axis, and the fields in the window")
print("%-22s %26s %8s %9s %9s"
      % ("scheme", "pocket(s) x/D_e", "length", "|w| p99", "M_min"))
for t, n in zip(TAGS, ROWS):
    q = Q[t]
    j0 = int(np.argmin(np.abs(q["y"])))
    m = np.ix_((q["x"] >= ZX[0]) & (q["x"] <= ZX[1]),
               (q["y"] >= ZY[0]) & (q["y"] <= ZY[1]))
    L = sum(b-a for a, b in q["pkt"])
    txt = "  ".join("%.3f-%.3f" % p for p in q["pkt"]) or "none"
    w = (q["x"] >= 0.6) & (q["x"] <= 2.2)
    print("%-22s %26s %8.3f %9.2f %9.4f"
          % (n, txt, L, np.percentile(q["omg"][m], 99),
             q["M"][w, j0].min()))
print("\nshared schlieren reference g_ref = %.1f" % GREF)

plt.rcParams.update({"font.size": 9.0, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})

COLS = [dict(key="M", cmap="viridis", lo=0.0, hi=2.4, sonic=True,
             cb=r"$\overline{M}$", ttl="mean Mach number", light=True),
        dict(key="omg", cmap="magma", lo=0.0, hi=25.0, sonic=False,
             cb=r"$|\omega|\,D_e/U_j$", ttl="vorticity magnitude", light=True),
        dict(key="sch", cmap="gray", lo=0.0, hi=1.0, sonic=False,
             cb="numerical schlieren", ttl="numerical schlieren", light=False)]

fig, AX = plt.subplots(5, 3, figsize=(10.6, 12.4), sharex=True, sharey=True)
fig.subplots_adjust(left=0.068, right=0.885, top=0.955, bottom=0.052,
                    wspace=0.085, hspace=0.075)

for c, spec in enumerate(COLS):
    nrm = Normalize(spec["lo"], spec["hi"])
    for r, (tag, name) in enumerate(zip(TAGS, ROWS)):
        q = Q[tag]
        ax = AX[r, c]
        my = (q["y"] >= ZY[0]-0.05) & (q["y"] <= ZY[1]+0.05)
        im = ax.pcolormesh(q["x"], q["y"][my], q[spec["key"]].T[my],
                           cmap=spec["cmap"], norm=nrm, shading="nearest",
                           rasterized=True)
        if spec["sonic"]:
            ax.contour(q["x"], q["y"][my], q["M"].T[my], levels=[1.0],
                       colors="w", linewidths=0.9, zorder=4)
        for a, b in q["pkt"]:                       # pocket extent on the axis
            ax.plot([a, b], [0, 0], "-", color="#ff3b30", lw=2.4, zorder=6,
                    solid_capstyle="butt")
        ax.set_xlim(*ZX)
        ax.set_ylim(*ZY)
        ax.set_aspect("equal")
        if c == 0:
            ax.set_ylabel("$y/D_e$", labelpad=1)
            ax.text(0.030, 0.055, name, transform=ax.transAxes, fontsize=8.8,
                    color="w", zorder=8)
        if r == 0:
            ax.set_title(spec["ttl"], fontsize=10.0, pad=5)
        if r == 4:
            ax.set_xlabel("$x/D_e$", labelpad=2)
    cax = fig.add_axes([0.895, 0.052 + (2-c)*0.302, 0.016, 0.268])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(spec["cb"], fontsize=9.6, labelpad=3)
    cb.ax.tick_params(labelsize=8.2)

AX[0, 0].plot([], [], "-", color="#ff3b30", lw=2.4,
              label=r"$\overline{M}<1$ on the axis")
AX[0, 0].plot([], [], "-", color="w", lw=0.9, label=r"sonic line")
AX[0, 0].legend(loc="upper left", fontsize=7.6, labelcolor="w",
                frameon=False, handlelength=1.6,
                bbox_to_anchor=(0.015, 0.995))

fig.savefig(HERE/"pocket_5x3.png", dpi=250)
fig.savefig(HERE/"pocket_5x3.pdf")
print("wrote pocket_5x3.png/pdf")
