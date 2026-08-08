#!/usr/bin/env python3
"""Mid-plane fields for the two AFD + HLLC cases: mean density, Mach number of
the mean field, density fluctuation, and instantaneous vorticity magnitude.

Each panel carries both cases: WENO-Z5 above the axis, TENO5 below it, sharing
one colour scale, so the difference between the two reconstructions is read
directly across the axis rather than by moving the eye between plots.

Derivatives on an AMR composite
-------------------------------
The slice is a finest-wins composite on a level-5 abscissa, so away from the
finest patches one coarse cell is replicated across rf = 2, 4, 8 ... output
cells. A plain gradient on that array returns zero inside every replicated
block and a spike at its edge - a grid of artificial filaments, not vorticity.
Every derivative here is therefore taken over +/- rf output cells and divided
by 2*rf*dx, with rf read per cell from the map the extractor recorded, which
is the derivative at the resolution the data actually has. The same applies to
d/dz, for which the extractor kept one plane either side offset by a cell of
the owning level, not of the output grid.

The vorticity panel is honest everywhere in the sense that it reports what
each level can resolve, but the finest structure only exists where the finest
level does; the shear layer and the first cells are level 5, the far field is
not.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

HERE = Path(__file__).resolve().parent
F = HERE/"fields"
RHO_J = 1.6413
U_J = 414.2
D = 0.0254
GAM = 1.4

TOP, BOT = "hllcweno", "hllcteno"
LAB = {TOP: "WENO-Z5", BOT: "TENO5"}


def load(tag, var, stack=False):
    z = np.load(F/("SLM_%s_%s.npz" % (tag, var)))
    a = (z["kstack"] if stack else z["field"]).astype(np.float64)
    return a, z


def ddx(f, rf, d, axis):
    """Centred difference over +/- rf cells, so the stencil spans one real
    cell of the level that owns each point."""
    out = np.zeros_like(f)
    for r in np.unique(rf[rf > 0]):
        r = int(r)
        fp = np.roll(f, -r, axis=axis)
        fm = np.roll(f, r, axis=axis)
        g = (fp - fm)/(2.0*r*d)
        sl = [slice(None)]*f.ndim
        sl[axis] = slice(0, r)
        g[tuple(sl)] = 0.0
        sl[axis] = slice(-r, None)
        g[tuple(sl)] = 0.0
        out = np.where(rf == r, g, out)
    return out


def quantities(tag):
    rm, z = load(tag, "DensityMEAN")
    rs, _ = load(tag, "DensitySQR")
    pm, _ = load(tag, "pressureMEAN")
    um, _ = load(tag, "x_velocityMEAN")
    vm, _ = load(tag, "y_velocityMEAN")
    wm, _ = load(tag, "z_velocityMEAN")
    rf = z["rfmap"]
    d = float(z["dx_um"])*1e-6

    rho = rm/RHO_J
    mach = np.sqrt(um**2 + vm**2 + wm**2)/np.sqrt(GAM*pm/rm)
    rrms = np.sqrt(np.clip(rs - rm**2, 0.0, None))/RHO_J

    U, _ = load(tag, "x_velocity", True)
    V, _ = load(tag, "y_velocity", True)
    W, _ = load(tag, "z_velocity", True)
    dz = (2.0*rf*d)
    wx = ddx(W[1], rf, d, 1) - (V[2]-V[0])/dz
    wy = (U[2]-U[0])/dz - ddx(W[1], rf, d, 0)
    wz = ddx(V[1], rf, d, 0) - ddx(U[1], rf, d, 1)
    omg = np.sqrt(wx**2 + wy**2 + wz**2)*D/U_J
    return dict(rho=rho, mach=mach, rrms=rrms, omg=omg,
                x=z["x"], y=z["y"], t=float(z["time"])*1e3, rf=rf)


Q = {t: quantities(t) for t in (TOP, BOT)}
for t in (TOP, BOT):
    q = Q[t]
    print("%-9s t=%.3f ms  rho %.2f-%.2f  Mach %.2f-%.2f  rho_rms %.3f-%.3f  "
          "|w| p50/p99/max %.1f/%.1f/%.1f"
          % (t, q["t"], q["rho"].min(), q["rho"].max(), q["mach"].min(),
             q["mach"].max(), q["rrms"].min(), q["rrms"].max(),
             np.percentile(q["omg"], 50), np.percentile(q["omg"], 99),
             q["omg"].max()))
print("level coverage of the window: " + ", ".join(
    "rf=%d %.1f%%" % (r, 100*np.mean(Q[TOP]["rf"] == r))
    for r in np.unique(Q[TOP]["rf"]) if r > 0))

PAN = [("rho",  r"$\overline{\rho}/\rho_j$",           "inferno", 0.35, 1.50),
       ("mach", r"$\overline{M}$",                     "viridis", 0.00, 2.40),
       ("rrms", r"$\rho_{rms}/\rho_j$",                "cividis", 0.00, 0.18),
       ("omg",  r"$|\omega|\,D_e/U_j$",                "magma",   0.00, 25.0)]

plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})
fig, AX = plt.subplots(2, 2, figsize=(12.4, 5.9))
fig.subplots_adjust(left=0.048, right=0.972, top=0.975, bottom=0.082,
                    wspace=0.255, hspace=0.215)

for ax, (key, lab, cmap, lo, hi) in zip(AX.ravel(), PAN):
    nrm = Normalize(lo, hi)
    qt, qb = Q[TOP], Q[BOT]
    up = qt["y"] >= 0
    dn = qb["y"] < 0
    im = ax.pcolormesh(qt["x"], qt["y"][up], qt[key].T[up], cmap=cmap,
                       norm=nrm, shading="nearest", rasterized=True)
    ax.pcolormesh(qb["x"], qb["y"][dn], qb[key].T[dn], cmap=cmap, norm=nrm,
                  shading="nearest", rasterized=True)
    ax.axhline(0.0, color="w", lw=0.8)
    ax.set_xlim(0, 7.2)
    ax.set_ylim(-1.6, 1.6)
    ax.set_aspect("equal")
    ax.set_xlabel("$x/D_e$", labelpad=1)
    ax.set_ylabel("$y/D_e$", labelpad=1)
    ax.text(0.012, 0.90, LAB[TOP], transform=ax.transAxes, color="w",
            fontsize=9.0)
    ax.text(0.012, 0.055, LAB[BOT], transform=ax.transAxes, color="w",
            fontsize=9.0)
    cb = fig.colorbar(im, ax=ax, pad=0.012, fraction=0.028, aspect=11)
    cb.set_label(lab, fontsize=10.5, labelpad=3)
    cb.ax.tick_params(labelsize=8)

fig.savefig(HERE/"scheme_fields.png", dpi=250)
fig.savefig(HERE/"scheme_fields.pdf")
print("\nwrote scheme_fields.png/pdf")
