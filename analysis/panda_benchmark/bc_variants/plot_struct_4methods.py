#!/usr/bin/env python3
"""Structure loss across the four scheme combinations, on one slice.

Rows are LLF+WENO-Z5, LLF+TENO5, AFD-HLLC+WENO-Z5, AFD-HLLC+TENO5, so the two
flux paths sit as adjacent pairs and the reconstruction alternates inside each
pair. Every figure uses one colour bar for all four rows, so a difference
between rows is a difference in the data and never in the mapping. Same grid,
same slice plane, same normalisation.

    struct4_vorticity.png   |omega| D_e/U_j, full view and the shear layer
    struct4_schlieren.png   numerical schlieren of the MEAN density
    struct4_lambda2.png     lambda_2 of the mean-square strain-rotation tensor
    struct4_q.png           Q criterion

What the slice can and cannot support
-------------------------------------
The extractor kept three planes of each velocity component, one cell of the
owning level either side of the slice, so the full velocity gradient tensor
du_i/dx_j is available on the plane and Q and lambda_2 are computed exactly,
not from an in-plane approximation. They are computed here as a *map on the
slice*. An isosurface is a different object and needs the volume; these files
do not carry it.

The stored instantaneous variables are the three velocity components only.
There is no instantaneous density, so the schlieren here is formed from the
mean density and reports the sharpness of the shock cells and of the mean
shear layer. It cannot show instantaneous contact surfaces or the braid
structure between rollers - that needs an instantaneous density field, which
is not in these slices.

Derivatives on an AMR composite
-------------------------------
The slice is a finest-wins composite on a level-5 abscissa, so away from the
finest patches one coarse cell is replicated across rf = 2, 4 or 8 output
cells. A plain gradient there returns zero inside every replicated block and a
spike at its edge, which is a lattice of artificial filaments rather than
vorticity. Every derivative is taken over +/- rf output cells and divided by
2 rf dx, with rf read per cell from the map the extractor recorded; d/dz uses
the planes the extractor kept one cell of the owning level either side.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent
F = HERE/"fields"
RHO_J, U_J, D = 1.6413, 414.2, 0.0254
YMAX = 1.20
XHI = 7.20
ZX = (1.45, 3.35)          # shear-layer enlargement
ZY = (0.20, 0.88)

TAGS = ["llfweno", "llfteno", "hllcweno", "hllcteno"]
ROWS = ["LLF + WENO-Z5", "LLF + TENO5", "AFD-HLLC + WENO-Z5", "AFD-HLLC + TENO5"]


def load(tag, var, stack=False):
    z = np.load(F/("SLM_%s_%s.npz" % (tag, var)))
    return (z["kstack"] if stack else z["field"]).astype(np.float64), z


def ddx(f, rf, d, axis):
    """Centred difference spanning one real cell of the level that owns it."""
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
    rf, d = z["rfmap"], float(z["dx_um"])*1e-6
    U, _ = load(tag, "x_velocity", True)
    V, _ = load(tag, "y_velocity", True)
    W, _ = load(tag, "z_velocity", True)
    dz = 2.0*rf*d

    # full velocity gradient tensor on the plane, J[i][j] = du_i/dx_j
    J = [[ddx(U[1], rf, d, 0), ddx(U[1], rf, d, 1), (U[2]-U[0])/dz],
         [ddx(V[1], rf, d, 0), ddx(V[1], rf, d, 1), (V[2]-V[0])/dz],
         [ddx(W[1], rf, d, 0), ddx(W[1], rf, d, 1), (W[2]-W[0])/dz]]

    wx, wy, wz = J[2][1]-J[1][2], J[0][2]-J[2][0], J[1][0]-J[0][1]
    omg = np.sqrt(wx*wx + wy*wy + wz*wz)

    S = [[0.5*(J[i][j]+J[j][i]) for j in range(3)] for i in range(3)]
    O = [[0.5*(J[i][j]-J[j][i]) for j in range(3)] for i in range(3)]
    # Q = 0.5(|Omega|^2 - |S|^2), the standard second invariant
    Q = 0.5*(sum(O[i][j]**2 for i in range(3) for j in range(3))
             - sum(S[i][j]**2 for i in range(3) for j in range(3)))

    # lambda_2: the median eigenvalue of S^2 + Omega^2, negative in a core
    M = np.empty(J[0][0].shape + (3, 3))
    for i in range(3):
        for j in range(3):
            M[..., i, j] = sum(S[i][k]*S[k][j] + O[i][k]*O[k][j]
                               for k in range(3))
    lam2 = np.linalg.eigvalsh(0.5*(M + np.swapaxes(M, -1, -2)))[..., 1]

    grho = np.hypot(ddx(rm, rf, d, 0), ddx(rm, rf, d, 1))
    sc = U_J/D
    return dict(omg=omg/sc, Q=Q/sc**2, lam2=lam2/sc**2, grho=grho,
                x=z["x"], y=z["y"], t=float(z["time"])*1e3)


Q = {t: quantities(t) for t in TAGS}
print("%-19s %8s %8s %9s %10s %10s"
      % ("scheme", "|w|_p99", "|w|_max", "Q_p99", "lam2_p01", "|grad rho|_p99.5"))
for t, n in zip(TAGS, ROWS):
    q = Q[t]
    print("%-19s %8.2f %8.2f %9.1f %10.1f %10.1f"
          % (n, np.percentile(q["omg"], 99), q["omg"].max(),
             np.percentile(q["Q"], 99), np.percentile(q["lam2"], 1),
             np.percentile(q["grho"], 99.5)))

# schlieren: one reference for all four rows so the rows stay comparable, taken
# as the 99.5th percentile over x > 0.8 D_e pooled across the four schemes. The
# lip is a discontinuity in the imposed profile and its gradient is far above
# anything downstream, so a maximum-based reference would leave every shock in
# the first per cent of the range.
pool = np.concatenate([Q[t]["grho"][Q[t]["x"] > 0.8].ravel() for t in TAGS])
GREF = np.nanpercentile(pool, 99.5)
for t in TAGS:
    Q[t]["sch"] = np.exp(-5.0*np.clip(Q[t]["grho"]/GREF, 0.0, 1.0))
print("\nshared schlieren reference g_ref = %.1f (99.5th percentile, x > 0.8)"
      % GREF)

plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True,
                     "legend.frameon": False})


def stack4(key, cmap, lo, hi, cbar, out, light_text=True, zoom=False):
    n = 8 if zoom else 4
    fig, AX = plt.subplots(n, 1, figsize=(9.4, 8.4 if not zoom else 12.6),
                           sharex=False)
    fig.subplots_adjust(left=0.078, right=0.872, top=0.988, bottom=0.048,
                        hspace=0.10 if not zoom else 0.16)
    nrm = Normalize(lo, hi)
    for k, (tag, name) in enumerate(zip(TAGS, ROWS)):
        q = Q[tag]
        m = np.abs(q["y"]) <= YMAX
        ax = AX[k if not zoom else 2*k]
        im = ax.pcolormesh(q["x"], q["y"][m], q[key].T[m], cmap=cmap, norm=nrm,
                           shading="nearest", rasterized=True)
        ax.set_xlim(0, XHI)
        ax.set_ylim(-YMAX, YMAX)
        ax.set_aspect("equal")
        ax.set_ylabel("$y/D_e$", labelpad=1)
        ax.text(0.009, 0.86, name, transform=ax.transAxes,
                color="w" if light_text else "k", fontsize=10)
        if zoom:
            ax.add_patch(Rectangle((ZX[0], ZY[0]), ZX[1]-ZX[0], ZY[1]-ZY[0],
                                   fill=False, edgecolor="#e8482c", lw=1.0,
                                   zorder=6))
            az = AX[2*k+1]
            mz = (q["y"] >= ZY[0]) & (q["y"] <= ZY[1])
            az.pcolormesh(q["x"], q["y"][mz], q[key].T[mz], cmap=cmap,
                          norm=nrm, shading="nearest", rasterized=True)
            az.set_xlim(*ZX)
            az.set_ylim(*ZY)
            az.set_aspect("equal")
            az.set_ylabel("$y/D_e$", labelpad=1)
            for sp in az.spines.values():
                sp.set_color("#e8482c")
                sp.set_linewidth(1.0)
        if k < 3:
            ax.set_xticklabels([])
    AX[-1].set_xlabel("$x/D_e$", labelpad=2)
    cax = fig.add_axes([0.888, 0.048, 0.018, 0.940])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(cbar, fontsize=11.5, labelpad=4)
    cb.ax.tick_params(labelsize=9)
    fig.savefig(HERE/(out+".png"), dpi=250)
    fig.savefig(HERE/(out+".pdf"))
    plt.close(fig)
    print("wrote %s.png/pdf" % out)


stack4("omg", "magma", 0.0, 25.0, r"$|\omega|\,D_e/U_j$",
       "struct4_vorticity", zoom=True)
stack4("sch", "gray", 0.0, 1.0, "numerical schlieren of "
       r"$\overline{\rho}$,   $\exp(-k|\nabla\overline{\rho}|/g_{ref})$",
       "struct4_schlieren", light_text=False)
# Limits from the pooled distribution of all four schemes, not from one of
# them: p0.1 and p99.9 of Q are -138 and +191, of lambda_2 -187 and +89, so a
# symmetric +/-80 keeps the vortex cores on scale and saturates only the
# shock-crossing tail. A wider range would flatten every roller to the
# midpoint and hide exactly the difference the figure is for.
stack4("lam2", "RdBu_r", -80.0, 80.0,
       r"$\lambda_2\,(D_e/U_j)^2$", "struct4_lambda2", light_text=False)
stack4("Q", "RdBu_r", -80.0, 80.0, r"$Q\,(D_e/U_j)^2$", "struct4_q",
       light_text=False)
