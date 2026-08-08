#!/usr/bin/env python3
"""yz cross-sections of the 3D L5 solution at the eight radial-comparison
stations (x/D_e = 0.6, 0.75, 0.9, 1.05, 1.25, 1.55, 2.0, 3.05), 2 x 4
panels per figure, same station order as paperfig_cmp_radial_rho:
  paperfig_yz8_rho       : instantaneous rho/rho_j            (viridis)
  paperfig_yz8_schlieren : exp(-8 |grad rho| / |grad rho|_99.9), full
                           gradient, global reference           (gray)
  paperfig_yz8_vort      : instantaneous |omega| D_e/U_j       (viridis)
  paperfig_yz8_rhomean   : <rho>/rho_j                        (viridis)
  paperfig_yz8_rhostd    : rho_rms/rho_j                       (magma)
Data: SIM_L5_yz8_multi.npz (Density/DensityMEAN/DensitySQR planes) and
SIM_L5_yz8_struct.npz (grho), plt09170, t = 10.16 ms."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams.update({
    "font.size": 9, "axes.linewidth": 0.6,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
})
HERE = Path(__file__).parent
RHO_J = 1.6413
D_E = 0.0254
U_J = 414.2
R = 0.85
STATIONS = [0.6, 0.75, 0.9, 1.05, 1.25, 1.55, 2.0, 3.05]


def key(st):
    return str(float(st)).replace(".", "p")


dm = np.load(HERE / "SIM_L5_yz8_multi.npz", allow_pickle=True)
WANT = list(dm["want"])


def gv(st, name):
    return dm["x" + key(st)][WANT.index(name)].astype(np.float64)


ds = np.load(HERE / "SIM_L5_yz8_struct.npz")
gkeys = {k for k in ds.files if k.startswith("grho_")}


def grho(st):
    for cand in ("grho_" + key(st), "grho_x" + key(st)):
        if cand in gkeys:
            return ds[cand].astype(np.float64)
    raise KeyError("no grho key for station %s in %s" % (st, sorted(gkeys)))


def struct(st, name):
    """Return one of the structural fields stored at a requested station."""
    for cand in (name + "_" + key(st), name + "_x" + key(st)):
        if cand in ds.files:
            return ds[cand].astype(np.float64)
    raise KeyError("no %s key for station %s" % (name, st))


def draw(fields, outname, cblabel, vmin, vmax, cmap):
    fig, axs = plt.subplots(2, 4, figsize=(7.2, 3.65), sharex=True,
                            sharey=True)
    for k, (stn, ax) in enumerate(zip(STATIONS, axs.flat)):
        F = fields[k]
        im = ax.imshow(F, origin="lower", extent=[-R, R, -R, R],
                       vmin=vmin, vmax=vmax, cmap=cmap, aspect="equal",
                       interpolation="nearest", rasterized=True)
        th = np.linspace(0, 2 * np.pi, 200)
        lc = "k" if cmap in ("gray", "gray_r") else "w"
        ax.plot(0.5 * np.cos(th), 0.5 * np.sin(th), color=lc, lw=0.5,
                ls=":")
        # Keep the station annotation outside the data region so it never
        # obscures the shock or shear-layer structure.
        ax.set_title("(%s) $x/D_e=%g$" % (chr(97 + k), stn),
                     loc="left", fontsize=8.2, pad=2.0)
        ax.set_xlim(-R, R); ax.set_ylim(-R, R)
        ax.set_xticks([-0.5, 0.0, 0.5])
        ax.set_yticks([-0.5, 0.0, 0.5])
        ax.tick_params(length=2.5, labelsize=8)
    for ax in axs[1]:
        ax.set_xlabel("$z/D_e$", labelpad=2)
    for ax in axs[:, 0]:
        ax.set_ylabel("$y/D_e$", labelpad=1)
    # Leave a full outer margin for the vertical colour-bar tick labels and
    # its rotated label.  This avoids clipping without changing the saved
    # canvas size or relying on bbox_inches="tight".
    fig.subplots_adjust(left=0.075, right=0.855, top=0.93, bottom=0.135,
                        wspace=0.08, hspace=0.28)
    cax = fig.add_axes([0.880, 0.155, 0.018, 0.70])
    cb = fig.colorbar(im, cax=cax, orientation="vertical")
    cb.set_label(cblabel, fontsize=8.5, labelpad=3)
    cb.ax.tick_params(length=2.5, labelsize=8)
    fig.savefig(HERE / (outname + ".png"), dpi=300)
    fig.savefig(HERE / (outname + ".pdf"), dpi=300)
    print("wrote", outname)


# instantaneous and mean density, common scale
rho_i = [gv(st, "Density") / RHO_J for st in STATIONS]
rho_m = [gv(st, "DensityMEAN") / RHO_J for st in STATIONS]
draw(rho_i, "paperfig_yz8_rho", "$\\rho/\\rho_j$", 0.35, 1.45, "viridis")
draw(rho_m, "paperfig_yz8_rhomean", "$\\langle\\rho\\rangle/\\rho_j$",
     0.35, 1.45, "viridis")

# numerical schlieren, one global reference over all eight stations
gs = [grho(st) for st in STATIONS]
gref = np.nanpercentile(np.concatenate([g.ravel() for g in gs]), 99.9)
draw([np.exp(-8.0 * g / gref) for g in gs], "paperfig_yz8_schlieren",
     "$\\exp(-8\\,|\\nabla\\rho|/|\\nabla\\rho|_{99.9})$", 0.0, 1.0,
     "gray")

# instantaneous vorticity magnitude, using the scale adopted in Chapter 3
vort = [struct(st, "om") * D_E / U_J for st in STATIONS]
draw(vort, "paperfig_yz8_vort", "$|\\omega|D_e/U_j$",
     0.0, 12.0, "viridis")

# density rms, common ceiling at the global 99.5th percentile
std = [np.sqrt(np.maximum(gv(st, "DensitySQR")
                          - gv(st, "DensityMEAN")**2, 0)) / RHO_J
       for st in STATIONS]
VS = float(np.ceil(np.nanpercentile(
    np.concatenate([s.ravel() for s in std]), 99.5) * 20) / 20)
draw(std, "paperfig_yz8_rhostd", "$\\rho_{\\mathrm{rms}}/\\rho_j$",
     0.0, VS, "magma")
