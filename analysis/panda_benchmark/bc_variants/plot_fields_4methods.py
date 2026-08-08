#!/usr/bin/env python3
"""One flow-field map per quantity, all four schemes stacked, shared scale.

Rows are the four combinations of flux path and reconstruction, ordered
LLF+WENO-Z5, LLF+TENO5, AFD-HLLC+WENO-Z5, AFD-HLLC+TENO5, so the two flux
paths sit as adjacent pairs and the reconstruction alternates inside each
pair. One colour bar serves all four rows: a difference between rows is a
difference in the data, never in the mapping.

Four quantities: mean density, Mach number of the mean field, density
fluctuation, and instantaneous vorticity magnitude.

All four schemes now come from one extraction path. The AFD pair was sliced
from its own L5 statistics frames; the LLF pair had its statistics plotfiles
pruned and was recovered from the surviving checkpoints by a zero-advance
restart (see dump_llf2.sh), which is why its slices exist on the same grid
with the same nine variables.

Derivatives on an AMR composite
-------------------------------
The slice is a finest-wins composite on a level-5 abscissa, so away from the
finest patches one coarse cell is replicated across rf = 2, 4 or 8 output
cells. A plain gradient there returns zero inside every replicated block and a
spike at its edge - a lattice of artificial filaments, not vorticity. Every
derivative is therefore taken over +/- rf output cells and divided by 2*rf*dx,
rf read per cell from the map the extractor recorded. d/dz uses the planes the
extractor kept one cell of the owning level either side, not one cell of the
output grid.
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
YMAX = 1.20

TAGS = ["llfweno", "llfteno", "hllcweno", "hllcteno"]
ROWS = ["LLF + WENO-Z5", "LLF + TENO5", "AFD-HLLC + WENO-Z5", "AFD-HLLC + TENO5"]


def load(tag, var, stack=False):
    z = np.load(F/("SLM_%s_%s.npz" % (tag, var)))
    return (z["kstack"] if stack else z["field"]).astype(np.float64), z


def ddx(f, rf, d, axis):
    """Centred difference spanning one real cell of the owning level."""
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


def quantities(tag):
    rm, z = load(tag, "DensityMEAN")
    rs, _ = load(tag, "DensitySQR")
    pm, _ = load(tag, "pressureMEAN")
    um, _ = load(tag, "x_velocityMEAN")
    vm, _ = load(tag, "y_velocityMEAN")
    wm, _ = load(tag, "z_velocityMEAN")
    rf = z["rfmap"]
    d = float(z["dx_um"])*1e-6

    U, _ = load(tag, "x_velocity", True)
    V, _ = load(tag, "y_velocity", True)
    W, _ = load(tag, "z_velocity", True)
    dz = 2.0*rf*d
    wx = ddx(W[1], rf, d, 1) - (V[2]-V[0])/dz
    wy = (U[2]-U[0])/dz - ddx(W[1], rf, d, 0)
    wz = ddx(V[1], rf, d, 0) - ddx(U[1], rf, d, 1)

    return dict(rho=rm/RHO_J,
                mach=np.sqrt(um**2+vm**2+wm**2)/np.sqrt(GAM*pm/rm),
                rrms=np.sqrt(np.clip(rs - rm*rm, 0.0, None))/RHO_J,
                omg=np.sqrt(wx**2+wy**2+wz**2)*D/U_J,
                x=z["x"], y=z["y"], t=float(z["time"])*1e3)


Q = {t: quantities(t) for t in TAGS}
for t, name in zip(TAGS, ROWS):
    q = Q[t]
    print("%-20s t=%6.3f ms  rho %.2f-%.2f  M %.2f-%.2f  rms %.3f-%.3f  "
          "|w| p99=%.1f max=%.1f"
          % (name, q["t"], q["rho"].min(), q["rho"].max(), q["mach"].min(),
             q["mach"].max(), q["rrms"].min(), q["rrms"].max(),
             np.percentile(q["omg"], 99), q["omg"].max()))

FIGS = [dict(key="rho",  cbar=r"$\overline{\rho}/\rho_j$", cmap="inferno",
             lo=0.35, hi=1.50, out="fields4_rhoMEAN"),
        dict(key="mach", cbar=r"$\overline{M}$", cmap="viridis",
             lo=0.00, hi=2.40, out="fields4_mach"),
        dict(key="rrms", cbar=r"$\rho_{rms}/\rho_j$", cmap="cividis",
             lo=0.00, hi=0.19, out="fields4_rhoRMS"),
        dict(key="omg",  cbar=r"$|\omega|\,D_e/U_j$", cmap="magma",
             lo=0.00, hi=25.0, out="fields4_vorticity")]

plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})

for spec in FIGS:
    fig, AX = plt.subplots(4, 1, figsize=(9.4, 8.4), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.078, right=0.872, top=0.988, bottom=0.058,
                        hspace=0.085)
    nrm = Normalize(spec["lo"], spec["hi"])
    dark = spec["key"] == "omg"
    for ax, tag, name in zip(AX, TAGS, ROWS):
        q = Q[tag]
        m = np.abs(q["y"]) <= YMAX
        im = ax.pcolormesh(q["x"], q["y"][m], q[spec["key"]].T[m],
                           cmap=spec["cmap"], norm=nrm, shading="nearest",
                           rasterized=True)
        ax.set_xlim(0, 7.2)
        ax.set_ylim(-YMAX, YMAX)
        ax.set_aspect("equal")
        ax.set_ylabel("$y/D_e$", labelpad=1)
        ax.text(0.009, 0.86, name, transform=ax.transAxes,
                color="w" if dark or spec["key"] != "rrms" else "k",
                fontsize=10)
    AX[-1].set_xlabel("$x/D_e$", labelpad=2)
    cax = fig.add_axes([0.888, 0.058, 0.018, 0.930])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(spec["cbar"], fontsize=11.5, labelpad=4)
    cb.ax.tick_params(labelsize=9)
    fig.savefig(HERE/(spec["out"]+".png"), dpi=250)
    fig.savefig(HERE/(spec["out"]+".pdf"))
    plt.close(fig)
    print("wrote %s.png/pdf" % spec["out"])
