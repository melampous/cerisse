#!/usr/bin/env python3
"""Warp-consistent radial comparison (zero per-station free parameters).

The 5-node PCHIP axial alignment P (sim -> exp coordinate) built on the
CENTRELINE (plot_axis_rho_pchip_warp.py, identical construction) is used
to select the axial sampling position of the SIMULATED radial profiles:
for an experimental station x_s the simulation is sampled at

    x_sample = P^{-1}(x_s)   (monotone map, numerically inverted),

i.e. the first shock cell is almost unmoved and the sampling point moves
progressively upstream downstream, following the non-linear stretching.
No radial information enters the map: if the warp-sampled profiles match
the measured ones, the phase error is a coherent axial stretching of the
whole field.  Diagnostic only; density values are not rescaled.

Per station the r.m.s. against the M_j=1.42 survey points (azimuthal-mean
profile, both experiment sides) is reported for (i) nominal sampling,
(ii) warp-predicted sampling, (iii) the per-station free best fit within
+-0.3 D_e (lower bound of what any axial shift can achieve).
Figure: same composition as paperfig_cmp_radial_rho (M_j=1.42 survey
filled, M_j=1.44 traverse mirrored light, y-cut red dashed, z-cut blue
dashed, +-0.05 D_e envelope) with the cuts sampled at P^{-1}(x_s)."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.interpolate import PchipInterpolator
from scipy.signal import find_peaks

HERE = Path(__file__).resolve().parent
D_E = 0.0254
RHO_J = 1.6413
N_PEAKS = 5
C3D, C2D = "#b3251e", "#1f4e9c"


def moving_average(values, coordinates, width):
    count = max(3, int(round(width / (coordinates[1] - coordinates[0]))))
    if count % 2 == 0:
        count += 1
    result = np.convolve(values, np.ones(count) / count, mode="same")
    result[:count // 2] = values[:count // 2]
    result[-(count // 2):] = values[-(count // 2):]
    return result


def refined_peaks(x, y, n, prominence=0.02, xmin=0.5):
    idx, _ = find_peaks(y, prominence=prominence)
    idx = idx[x[idx] > xmin][:n]
    out = []
    for i in idx:
        d = y[i - 1] - 2.0 * y[i] + y[i + 1]
        s = 0.5 * (y[i - 1] - y[i + 1]) / d if abs(d) > 0 else 0.0
        out.append(x[i] + np.clip(s, -1, 1) * 0.5 * (x[i + 1] - x[i - 1]))
    return np.array(out)


# ---------------- axis warp, identical construction --------------------
sim = np.load(HERE / "AX3T_L5.npz", allow_pickle=True)
exp_ax = np.load(HERE.parent / "panda_m142den_axis_EPAPS.npz")
x_ax = sim["x"] / D_E
rho_ax = moving_average(sim["ax_DensityMEAN"] / RHO_J, x_ax, 0.05)
sim_pk = refined_peaks(x_ax, rho_ax, N_PEAKS)
exp_pk = refined_peaks(exp_ax["x"], exp_ax["rho"], N_PEAKS)
warp = PchipInterpolator(np.concatenate(([0.0], sim_pk)),
                         np.concatenate(([0.0], exp_pk)))
end_slope = float(warp.derivative()(sim_pk[-1]))
x_warped = np.where(x_ax <= sim_pk[-1],
                    warp(np.minimum(x_ax, sim_pk[-1])),
                    exp_pk[-1] + end_slope * (x_ax - sim_pk[-1]))


def inv_warp(x_exp_coord):
    """P^{-1}: exp coordinate -> sim coordinate (monotone)."""
    return float(np.interp(x_exp_coord, x_warped, x_ax))


# ---------------- data -------------------------------------------------
f5 = np.load(HERE.parent / "panda_fig5a_mj144_EPAPS.npz", allow_pickle=True)
RHOJ_EXP = float(f5["rhoj"])
fm = np.load(HERE.parent / "panda_m142den_field_EPAPS.npz")
xm, rm, Fm = fm["x"], fm["r"], fm["rho"]
sly = np.load(HERE / "SIM_L5_rhoMEAN_z0.npz")
slz = np.load(HERE / "SIM_L5_rhoMEAN_y0.npz")
xs3, rr3 = sly["x"], sly["y"]
cutY = sly["rho"] / RHO_J
cutZ = slz["rho"] / RHO_J
rpos = rr3 >= 0
rfold = rr3[rpos]


def azim_mean(j):
    return 0.25 * (cutY[j][rpos] + cutY[j][::-1][rpos]
                   + cutZ[j][rpos] + cutZ[j][::-1][rpos])


def rms_at(jj, rho42):
    return np.sqrt(np.mean((np.interp(np.abs(rm), rfold, azim_mean(jj))
                            - rho42)**2))


STATIONS = [0.6, 0.75, 0.9, 1.05, 1.25, 1.55, 2.0, 3.05]
OFF = 0.05
SEARCH = 0.3

print("station  col_x   x_sample=P^-1  shift    RMS_nom  RMS_warp  RMS_free")
rows = []
for stn in STATIONS:
    icol = np.argmin(np.abs(xm - stn))
    xcol = xm[icol]; rho42 = Fm[icol]
    j_nom = np.argmin(np.abs(xs3 - xcol))
    xw = inv_warp(xcol)
    j_w = np.argmin(np.abs(xs3 - xw))
    cand = np.where(np.abs(xs3 - xcol) <= SEARCH + 1e-9)[0]
    err = [rms_at(jj, rho42) for jj in cand]
    r_nom, r_w, r_free = rms_at(j_nom, rho42), rms_at(j_w, rho42), min(err)
    rows.append((stn, xcol, xw, j_w))
    print("%6.2f   %5.2f     %5.2f       %+0.3f    %.3f    %.3f     %.3f"
          % (stn, xcol, xw, xw - xcol, r_nom, r_w, r_free))

# ---------------- figure ----------------------------------------------
plt.rcParams.update({
    "font.size": 9, "axes.linewidth": 0.6,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
})
fig, axs = plt.subplots(2, 4, figsize=(7.2, 4.15), sharex=True, sharey=True)
for k, ((stn, xcol, xw, j_w), ax) in enumerate(zip(rows, axs.flat)):
    ax.set_axisbelow(True)
    ax.grid(True, ls="--", lw=0.5, color="0.87")
    jsel = np.abs(xs3 - xw) <= OFF + 1e-9
    stack = np.vstack([cutY[jsel], cutZ[jsel]])
    ax.fill_between(rr3, stack.min(0), stack.max(0), color=C3D, alpha=0.15,
                    lw=0, zorder=2)
    ax.plot(rr3, cutY[j_w], color=C3D, lw=1.1, ls="--", zorder=4)
    ax.plot(rr3, cutZ[j_w], color=C2D, lw=1.1, ls="--", zorder=4)
    e = f5["x%s" % str(stn)]
    ax.plot(e[:, 0], e[:, 1] / RHOJ_EXP, "o", color="0.70", ms=2.6,
            mfc="none", mew=0.7, zorder=4)
    ax.plot(-e[:, 0], e[:, 1] / RHOJ_EXP, "o", color="0.70", ms=2.6,
            mfc="none", mew=0.7, zorder=4)
    icol = np.argmin(np.abs(xm - stn))
    ax.plot(rm, Fm[icol], "o", color="k", ms=3.2, mew=0, zorder=5)
    ax.text(0.05, 0.06,
            "(%s) $x/D_e=%g$\n$P^{-1}\\!: %0.2f$" % (chr(97 + k), stn, xw),
            transform=ax.transAxes, ha="left", va="bottom", fontsize=7.5,
            zorder=6,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.75,
                      pad=1.2))
    ax.set_xlim(-0.85, 0.85)
    ax.set_xticks([-0.5, 0.0, 0.5])
    ax.tick_params(length=2.5)
for ax in axs[1]:
    ax.set_xlabel("$r/D_e$", labelpad=2)
for ax in axs[:, 0]:
    ax.set_ylabel("$\\langle\\rho\\rangle/\\rho_j$", labelpad=2)
axs[0, 0].set_ylim(0.3, 1.5)
handles = [
    Line2D([], [], marker="o", ls="", color="k", ms=3.2,
           label="Exp. $M_j=1.42$ survey"),
    Line2D([], [], marker="o", ls="", color="0.70", mfc="none", mew=0.7,
           ms=2.6, label="Exp. $M_j=1.44$ traverse (mirrored)"),
    Line2D([], [], color=C3D, lw=1.1, ls="--",
           label="3D L5, $y$-cut at $P^{-1}(x_s)$"),
    Line2D([], [], color=C2D, lw=1.1, ls="--",
           label="3D L5, $z$-cut at $P^{-1}(x_s)$"),
    Patch(facecolor=C3D, alpha=0.15, label="$\\pm0.05\\,D_e$ envelope"),
]
fig.legend(handles=handles, fontsize=6.5, frameon=False, ncol=5,
           loc="upper center", bbox_to_anchor=(0.5, 1.0),
           handlelength=1.3, columnspacing=0.8, handletextpad=0.4)
fig.tight_layout(w_pad=0.5, h_pad=0.5, rect=[0, 0, 1, 0.955])
out = HERE / "diagnostic_radial_warp_sampled"
fig.savefig(out.with_suffix(".png"), dpi=300)
fig.savefig(out.with_suffix(".pdf"))
print("wrote", out.with_suffix(".png").name)
