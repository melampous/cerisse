#!/usr/bin/env python3
"""2D-vs-3D vorticity-dynamics figures (matched finest spacing 99.2 um:
axisymmetric L4 at t=8.0 ms vs 3D L5 z=0 plane at t=10.16 ms; single
instantaneous snapshots, level-by-level gradients, finest-wins).

  paperfig_dim_vortcomp   : schlieren | cylindrical omega_theta (signed);
                            row 1 = axisymmetric, row 2 = 3D.
  paperfig_dim_production : P*_s, P*_d, P*_b maps, common symlog scale.
  paperfig_dim_integrals  : (a) net Pi_k(x), (b) activity Pi_k(x),
                            (c) F_3D(x). Axisymmetric weight 2*pi*r dr;
                            3D from the +-y-averaged z=0 plane (two
                            azimuthal samples, plane-based estimate)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from pathlib import Path

plt.rcParams.update({
    "font.size": 9, "axes.linewidth": 0.6,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
})
HERE = Path(__file__).parent
D, UJ = 0.0254, 414.2
S1 = D / UJ
SP = D**3 / UJ**3
XW, RW = 4.2, 1.25

d3 = np.load(HERE / "VBUD3D_L5.npz")
d2 = np.load(HERE / "VBUD2D_L4.npz")
x3, y3 = d3["x"], d3["y"]
r2, a2 = d2["r"], d2["xax"]

def mirror2d(F, flip=False):
    """(r, x) one-sided -> (x, +-r) mirrored, matching the 3D plane."""
    G = F.T                                    # (x, r)
    lower = -G[:, ::-1] if flip else G[:, ::-1]
    return np.concatenate([lower, G], axis=1)  # (x, 2r)

y2m = np.concatenate([-r2[::-1], r2])

# ================= Fig 1: vorticity components ======================
sch2 = mirror2d(d2["grho"])
sch3 = d3["grho"]
gref = np.nanpercentile(np.concatenate([sch2.ravel(), sch3.ravel()]), 99.9)
# Plot the scalar cylindrical azimuthal component on both sides of the axis.
# On z=0, omega_theta = sgn(y) omega_z.  The axisymmetric scalar is even
# under reflection.  This differs from plotting the Cartesian vector
# component omega_z, which changes sign across the axis.
omth2 = mirror2d(d2["omth"], flip=False) * S1
sgn_y3 = np.where(y3 >= 0.0, 1.0, -1.0)
omth3 = d3["omz"] * sgn_y3[None, :] * S1

fig, axs = plt.subplots(2, 2, figsize=(6.5, 4.4), sharex=True,
                        sharey=True)
VLIM = 12.0
ROWS = [("Axisymmetric cylindrical\nLevel 4, $t=8.00$ ms",
         (np.exp(-8 * sch2 / gref), 0, 1, "gray"),
         (omth2, -VLIM, VLIM, "RdBu_r")),
        ("Three-dimensional Cartesian\nLevel 5, $z=0$, $t=10.16$ ms",
         (np.exp(-8 * sch3 / gref), 0, 1, "gray"),
         (omth3, -VLIM, VLIM, "RdBu_r"))]
HEADS = ["Numerical schlieren", "$\\omega_\\theta D_e/U_j$"]
ims = {}
for rI, (rowlab, *panels) in enumerate(ROWS):
    xx, yy = (a2, y2m) if rI == 0 else (x3, y3)
    for cI, (F, v0, v1, cm) in enumerate(panels):
        ax = axs[rI, cI]
        ims[cI] = ax.imshow(F.T, origin="lower",
                            extent=[xx[0], xx[-1], yy[0], yy[-1]],
                            vmin=v0, vmax=v1, cmap=cm, aspect="equal",
                            interpolation="nearest", rasterized=True)
        if rI == 0:
            ax.set_title(HEADS[cI], fontsize=8.5, pad=3)
        tcol = "k" if cm == "gray" or (cm == "RdBu_r") else "w"
        ax.text(0.02, 0.95, "(%s)" % chr(97 + 2 * rI + cI),
                transform=ax.transAxes, ha="left", va="top", fontsize=8,
                color=tcol,
                bbox=dict(facecolor="white" if tcol == "k" else "black",
                          edgecolor="none", alpha=0.5, pad=1.0))
        ax.set_xlim(0, XW); ax.set_ylim(-RW, RW)
        ax.tick_params(length=2.2, labelsize=7.5)
for rI in range(2):
    axs[rI, 0].set_ylabel("$y/D_e$", labelpad=1)
for ax in axs[1]:
    ax.set_xlabel("$x/D_e$", labelpad=2)
fig.subplots_adjust(left=0.095, right=0.950, top=0.93, bottom=0.20,
                    wspace=0.08, hspace=0.06)
for x0, key, lab in ((0.13, 0, "Numerical schlieren"),
                     (0.575, 1, "$\\omega_\\theta D_e/U_j$")):
    cax = fig.add_axes([x0, 0.085, 0.34, 0.02])
    cb = fig.colorbar(
        ims[key],
        cax=cax,
        orientation="horizontal",
        extend="both" if key == 1 else "neither",
    )
    cb.set_label(lab, fontsize=8, labelpad=1)
    cb.ax.tick_params(length=2, labelsize=7.5)
fig.savefig(HERE / "paperfig_dim_vortcomp.png", dpi=280)
fig.savefig(HERE / "paperfig_dim_vortcomp.pdf", dpi=280)
print("wrote paperfig_dim_vortcomp")

# ================= Fig 2: production maps ===========================
norm = SymLogNorm(linthresh=10.0, vmin=-1e3, vmax=1e3, base=10)
P2 = [mirror2d(d2[k]) * SP for k in ("Ps", "Pd", "Pb")]
P3 = [d3[k] * SP for k in ("Ps", "Pd", "Pb")]
HEADS = ["$P_s^{*}$ (stretching)", "$P_d^{*}$ (dilatation)",
         "$P_b^{*}$ (baroclinic)"]
fig, axs = plt.subplots(2, 3, figsize=(8.6, 4.4), sharex=True,
                        sharey=True)
im = None
for rI, (Ps, xx, yy, rowlab) in enumerate(
        [(P2, a2, y2m, "2D axisym. L4"), (P3, x3, y3, "3D L5 ($z=0$)")]):
    for cI, F in enumerate(Ps):
        ax = axs[rI, cI]
        im = ax.imshow(F.T, origin="lower",
                       extent=[xx[0], xx[-1], yy[0], yy[-1]],
                       norm=norm, cmap="RdBu_r", aspect="equal",
                       interpolation="nearest", rasterized=True)
        if rI == 0:
            ax.set_title(HEADS[cI], fontsize=8.5, pad=3)
        ax.text(0.02, 0.95, "(%s)" % chr(97 + 3 * rI + cI),
                transform=ax.transAxes, ha="left", va="top", fontsize=8,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.6,
                          pad=1.0))
        ax.set_xlim(0, XW); ax.set_ylim(-RW, RW)
        ax.tick_params(length=2.2, labelsize=7.5)
    axs[rI, 0].set_ylabel("$r/D_e$", labelpad=1)
    axs[rI, 0].text(-0.30, 0.5, rowlab, transform=axs[rI, 0].transAxes,
                    rotation=90, ha="center", va="center", fontsize=7.5)
for ax in axs[1]:
    ax.set_xlabel("$x/D_e$", labelpad=2)
fig.subplots_adjust(left=0.09, right=0.985, top=0.93, bottom=0.20,
                    wspace=0.08, hspace=0.06)
cax = fig.add_axes([0.30, 0.085, 0.40, 0.02])
cb = fig.colorbar(im, cax=cax, orientation="horizontal")
cb.set_label("$P_k^{*} = P_k\\,D_e^3/U_j^3$  (symlog)", fontsize=8,
             labelpad=1)
cb.ax.tick_params(length=2, labelsize=7)
fig.savefig(HERE / "paperfig_dim_production.png", dpi=280)
fig.savefig(HERE / "paperfig_dim_production.pdf", dpi=280)
print("wrote paperfig_dim_production")

# ================= Fig 3: axial integrals ===========================
def boxfilt(a, x, wD):
    w = max(3, int(round(wD / (x[1] - x[0]))))
    if w % 2 == 0:
        w += 1
    out = np.convolve(a, np.ones(w) / w, mode="same")
    out[:w//2] = a[:w//2]; out[-(w//2):] = a[-(w//2):]
    return out

# axisymmetric: integral 2*pi*r dr over the one-sided (r, x) arrays
rr_m = r2 * D
w2 = 2 * np.pi * rr_m * (r2[1] - r2[0]) * D          # ring areas
def int2d(F):
    return np.nansum(F * w2[:, None], axis=0)        # (x,)

# 3D: fold +-y, same axisymmetric weight (plane-based estimate)
jc = len(y3) // 2
rfold3 = y3[jc:]
wr3 = 2 * np.pi * (rfold3 * D) * (y3[1] - y3[0]) * D
def int3d(F):
    Ff = 0.5 * (F[:, jc:] + F[:, :jc][:, ::-1])
    return np.nansum(Ff * wr3[None, :], axis=1)      # (x,)

CN = D / UJ**3
COLS = {"Ps": "#b3251e", "Pd": "#1f4e9c", "Pb": "#2a7f62"}
LABK = {"Ps": "$\\Pi_s$", "Pd": "$\\Pi_d$", "Pb": "$\\Pi_b$"}
fig, (aN, aA, aF) = plt.subplots(1, 3, figsize=(8.6, 2.8))
for ax in (aN, aA, aF):
    ax.set_axisbelow(True)
    ax.grid(True, ls="--", lw=0.5, color="0.87")
    ax.tick_params(length=2.5, labelsize=8)
for k in ("Ps", "Pd", "Pb"):
    aN.plot(x3, boxfilt(CN * int3d(d3[k]), x3, 0.1), color=COLS[k],
            lw=1.2, label=LABK[k] + " 3D")
    aN.plot(a2, boxfilt(CN * int2d(d2[k]), a2, 0.1), color=COLS[k],
            lw=1.0, ls="--", label=LABK[k] + " 2D")
    aA.plot(x3, boxfilt(CN * int3d(np.abs(d3[k])), x3, 0.1),
            color=COLS[k], lw=1.2)
    aA.plot(a2, boxfilt(CN * int2d(np.abs(d2[k])), a2, 0.1),
            color=COLS[k], lw=1.0, ls="--")
aN.axhline(0, color="0.3", lw=0.6)
aN.set_xlabel("$x/D_e$", labelpad=2)
aN.set_ylabel("$\\Pi_k^{\\mathrm{net}}$", labelpad=2)
aN.set_xlim(0, XW)
aN.legend(fontsize=6.0, frameon=False, ncol=2, loc="upper left",
          handlelength=1.4, columnspacing=0.8)
aN.text(0.03, 0.04, "(a)", transform=aN.transAxes, fontsize=9)
aA.set_xlabel("$x/D_e$", labelpad=2)
aA.set_ylabel("$\\Pi_k^{\\mathrm{act}}$", labelpad=2)
aA.set_xlim(0, XW)
aA.set_yscale("log")
aA.text(0.03, 0.04, "(b)", transform=aA.transAxes, fontsize=9)
# F_3D
naz2 = d3["omx"]**2 + d3["omy"]**2
tot = naz2 + d3["omz"]**2
F3D = int3d(naz2) / np.maximum(int3d(tot), 1e-30)
aF.plot(x3, boxfilt(F3D, x3, 0.1), color="k", lw=1.3,
        label="3D L5")
aF.axhline(2.0 / 3.0, color="0.4", lw=0.7, ls=":")
aF.text(3.9, 0.675, "isotropic $2/3$", ha="right", va="bottom",
        fontsize=7, color="0.35")
aF.axhline(0, color="#1f4e9c", lw=1.0, ls="--")
aF.text(3.9, 0.02, "2D axisym. $\\equiv 0$", ha="right", va="bottom",
        fontsize=7, color="#1f4e9c")
aF.set_xlabel("$x/D_e$", labelpad=2)
aF.set_ylabel("$F_{\\mathrm{3D}}$", labelpad=2)
aF.set_xlim(0, XW); aF.set_ylim(0, 0.85)
aF.text(0.03, 0.90, "(c)", transform=aF.transAxes, fontsize=9)
fig.tight_layout(w_pad=1.0)
fig.savefig(HERE / "paperfig_dim_integrals.png", dpi=300)
fig.savefig(HERE / "paperfig_dim_integrals.pdf")
print("wrote paperfig_dim_integrals")
# numbers for the section text
i1 = np.argmin(np.abs(x3 - 1.0)); i2 = np.argmin(np.abs(x3 - 2.0))
i4 = np.argmin(np.abs(x3 - 4.0))
F3s = boxfilt(F3D, x3, 0.1)
print("F_3D at x=1,2,4:", F3s[i1].round(3), F3s[i2].round(3),
      F3s[i4].round(3))
half = x3[np.argmax(F3s > 0.5)] if (F3s > 0.5).any() else np.nan
print("F_3D crosses 0.5 at x =", round(float(half), 2))
PsA3 = CN * int3d(np.abs(d3["Ps"])); PsA2 = CN * int2d(np.abs(d2["Ps"]))
print("activity ratio Pi_s 3D/2D at x=2: %.1f"
      % (boxfilt(PsA3, x3, 0.1)[i2]
         / boxfilt(PsA2, a2, 0.1)[np.argmin(np.abs(a2 - 2.0))]))
