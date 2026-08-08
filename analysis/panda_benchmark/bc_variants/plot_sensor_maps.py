#!/usr/bin/env python3
"""Where the two dissipative branches actually fire, evaluated offline.

Both switches that decide a face's dissipation are algebraic functions of the
instantaneous primitive state in a one-dimensional stencil, so they can be
reproduced exactly from the slice now that instantaneous density and pressure
are on disk for all five formulations. This turns the mechanism argument from
a reading of the source into a measurement.

LLF first-order fallback            archived implementation in src/rhs/Weno_old.h
    for the two cells c adjacent to the face, along the face normal,
        pocket |= rho[c] < 0.1 * max(rho[c-1], rho[c+1])
        pocket |= p[c]   < 0.1 * max(p[c-1],   p[c+1])
    A flagged face in that historical implementation drops the whole
    WENO/TENO reconstruction and takes the
    first-order Rusanov flux. The comment in the source calls it a
    rarefaction/near-vacuum detector, not a shock detector: a shock is a
    monotone jump, not a local minimum.

AFD smoothness sensor               src/rhs/AfdIBM.h:1113-1149, threshold 0.08
    sound_scale = max |c| over the six cells iv-3 .. iv+2
    for each of rho, p, u, v, w:
        field_scale = max |q| over the same six
        floor = sound_scale for the velocities, field_scale otherwise
        for the four triplets (a,b,c) at offsets (-3,-2,-1) .. (0,1,2):
            |a - 2b + c| / (|a| + 2|b| + |c| + floor) > 0.08  ->  NOT smooth
    A non-smooth face is sent to characteristic_llf_high_order, which is
    algebraically the LLF operator. With cns.afd_shock_llf = 1 and
    cns.afd_smoothness_threshold = 0.08, both defaults and neither overridden
    in the campaign inputs, this fraction is the share of faces on which
    AFD-HLLC is running as LLF.

Two caveats, both printed on the figure
---------------------------------------
The stencils are evaluated with offsets scaled by the local refinement ratio
read from rfmap, so each test spans the cells of the level that owns the data,
as the solver does. Without that scaling every replicated coarse block would
return a second difference of exactly zero and be reported as smooth, which
would understate the non-smooth fraction wherever the mesh is coarse - that is
most of the domain beyond 3 D_e.

The solver evaluates these per level at every Runge-Kutta stage. What is
measured here is one snapshot on a composite, so it is an instantaneous
activation fraction, not a time-integrated count.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

HERE = Path(__file__).resolve().parent
F = HERE/"fields"
GAM, YMAX, XHI = 1.4, 1.20, 7.20
RAREFY, THRESH = 0.10, 0.08

TAGS = ["llfweno", "llfteno", "hllcweno", "hllcteno", "skewjst"]
ROWS = ["LLF + WENO-Z5", "LLF + TENO5", "AFD-HLLC + WENO-Z5",
        "AFD-HLLC + TENO5", "Skew-symmetric + JST"]


def load(tag, var):
    z = np.load(F/("SLM_%s_%s.npz" % (tag, var)))
    return z["field"].astype(np.float64), z


def shift(a, k, rf, axis):
    """a evaluated k cells of the OWNING level away, along one axis."""
    out = np.empty_like(a)
    for r in np.unique(rf[rf > 0]):
        r = int(r)
        s = np.roll(a, -k*r, axis=axis)
        if k*r > 0:
            sl = [slice(None)]*a.ndim
            sl[axis] = slice(-k*r, None)
            s[tuple(sl)] = np.nan
        elif k*r < 0:
            sl = [slice(None)]*a.ndim
            sl[axis] = slice(0, -k*r)
            s[tuple(sl)] = np.nan
        out = np.where(rf == r, s, out)
    return out


def sensors(tag):
    rho, z = load(tag, "Density")
    p, _ = load(tag, "pressure")
    u, _ = load(tag, "x_velocity")
    v, _ = load(tag, "y_velocity")
    w, _ = load(tag, "z_velocity")
    rf = z["rfmap"]
    c = np.sqrt(GAM*p/rho)

    pocket = np.zeros(rho.shape, bool)
    nonsmooth = np.zeros(rho.shape, bool)
    for ax in (0, 1):
        # --- LLF rarefaction pocket, the two cells adjacent to the face
        for side in (-1, 0):
            for q in (rho, p):
                qc = shift(q, side, rf, ax)
                qm = shift(q, side-1, rf, ax)
                qp = shift(q, side+1, rf, ax)
                with np.errstate(invalid="ignore"):
                    pocket |= np.nan_to_num(
                        qc < RAREFY*np.maximum(qm, qp), nan=False)
        # --- AFD smoothness, six-cell stencil iv-3 .. iv+2
        cs = np.nanmax(np.stack([np.abs(shift(c, m-3, rf, ax))
                                 for m in range(6)]), axis=0)
        for fld, isvel in ((rho, False), (p, False),
                           (u, True), (v, True), (w, True)):
            st = [shift(fld, m-3, rf, ax) for m in range(6)]
            fs = np.nanmax(np.abs(np.stack(st)), axis=0)
            floor = cs if isvel else fs
            for m in range(4):
                a, b, cc = st[m], st[m+1], st[m+2]
                den = np.abs(a) + 2*np.abs(b) + np.abs(cc) + floor
                with np.errstate(invalid="ignore", divide="ignore"):
                    bad = np.abs(a - 2*b + cc)/den > THRESH
                nonsmooth |= np.nan_to_num(bad, nan=False)
    return dict(pocket=pocket, nonsmooth=nonsmooth, x=z["x"], y=z["y"],
                t=float(z["time"])*1e3, rf=rf)


Q = {t: sensors(t) for t in TAGS}
inb = None
print("activation fractions inside |y| < 1.2 D_e, x < 7.2 D_e")
print("%-22s %14s %16s" % ("scheme", "LLF fallback", "AFD non-smooth"))
for t, n in zip(TAGS, ROWS):
    q = Q[t]
    if inb is None:
        inb = np.ix_(q["x"] <= XHI, np.abs(q["y"]) <= YMAX)
    print("%-22s %13.4f%% %15.2f%%"
          % (n, 100*q["pocket"][inb].mean(), 100*q["nonsmooth"][inb].mean()))

print("\nAFD non-smooth fraction by station (the AFD pair is the one it governs)")
print("%-22s %8s %8s %8s %8s %8s" % ("scheme", "0-1D", "1-2D", "2-3D", "3-5D", "5-7D"))
for t, n in zip(TAGS, ROWS):
    q = Q[t]
    row = []
    for a, b in ((0, 1), (1, 2), (2, 3), (3, 5), (5, 7.2)):
        m = np.ix_((q["x"] >= a) & (q["x"] < b), np.abs(q["y"]) <= YMAX)
        row.append(100*q["nonsmooth"][m].mean())
    print("%-22s %7.1f%% %7.1f%% %7.1f%% %7.1f%% %7.1f%%" % (n, *row))

plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})

SPEC = [dict(key="nonsmooth", out="sensor5_afd_nonsmooth", col="#D55E00",
             ttl=r"AFD smoothness sensor fires  (face runs as LLF), "
                 r"threshold 0.08",
             note="fraction of the panel flagged"),
        dict(key="pocket", out="sensor5_llf_fallback", col="#8B1A1A",
             ttl=r"LLF first-order fallback fires  "
                 r"($\rho$ or $p$ below $0.1\times$ its neighbours)",
             note="fraction of the panel flagged")]

for s in SPEC:
    cmap = ListedColormap(["#f7f7f7", s["col"]])
    nrm = BoundaryNorm([0, 0.5, 1], 2)
    fig, AX = plt.subplots(5, 1, figsize=(9.4, 10.4), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.078, right=0.985, top=0.952, bottom=0.052,
                        hspace=0.085)
    for ax, t, n in zip(AX, TAGS, ROWS):
        q = Q[t]
        m = np.abs(q["y"]) <= YMAX
        ax.pcolormesh(q["x"], q["y"][m], q[s["key"]].T[m].astype(float),
                      cmap=cmap, norm=nrm, shading="nearest", rasterized=True)
        ax.set_xlim(0, XHI)
        ax.set_ylim(-YMAX, YMAX)
        ax.set_aspect("equal")
        ax.set_ylabel("$y/D_e$", labelpad=1)
        ax.text(0.009, 0.87, "%s      %.2f %% of the panel"
                % (n, 100*q[s["key"]][inb].mean()), transform=ax.transAxes,
                fontsize=9.6)
    AX[-1].set_xlabel("$x/D_e$", labelpad=2)
    fig.suptitle(s["ttl"], fontsize=10.5, y=0.985)
    fig.savefig(HERE/(s["out"]+".png"), dpi=250)
    fig.savefig(HERE/(s["out"]+".pdf"))
    plt.close(fig)
    print("wrote %s.png/pdf" % s["out"])
