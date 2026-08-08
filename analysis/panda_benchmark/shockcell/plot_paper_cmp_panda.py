#!/usr/bin/env python3
"""Experiment-vs-computation mean-density comparison (selected grids).
Fig 1 (paperfig_cmp_axis_rho): centreline <rho>/rho_j — primary reference =
cycle mean of the 37 phase-conditioned M_j=1.42 profiles (m142den, filled
circles).
Two finest grids of each family only: Cartesian Levels 4/5 (solid), RZ
Levels 5/6 (dashed),
finest level emphasised; every curve individually identified in the legend.
Fig 2 (paperfig_cmp_radial_rho): Cartesian Level 5 radial profiles at eight stations,
full +-r extent. Primary experiment = cycle mean of the M_j=1.42 density
experiment (m142den, 37 phase-conditioned fields, two-sided
r = -0.63..+0.70, nearest
measurement column |Dx| <= 0.04); secondary = the M_j=1.44 experiment
(fig5a, one-sided,
mirrored, light). Simulation: two independent diametral cuts, along y
(z=0 plane, red dashed) and along z (y=0 plane, blue dashed), unfolded;
shaded band = min-max envelope of both cuts over |x - x_station| <=
0.05 D_e. The M_j=1.42 data are distributed already normalised by rho_j;
the M_j=1.44 measurements are normalised by their rho_j = 1.62."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from pathlib import Path

plt.rcParams.update({
    "font.size": 9, "axes.linewidth": 0.6,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
})
HERE = Path(__file__).parent
D = 0.0254
RHO_J = 1.6413
C3D, C2D = "#b3251e", "#1f4e9c"

def boxfilt(a, x, wD):
    w = max(3, int(round(wD / (x[1] - x[0]))))
    if w % 2 == 0:
        w += 1
    out = np.convolve(a, np.ones(w) / w, mode="same")
    out[:w//2] = a[:w//2]; out[-(w//2):] = a[-(w//2):]
    return out

# ---------- Fig 1: centreline ----------
xd = np.load(HERE.parent / "panda_m142den_axis_EPAPS.npz")
a3 = np.load(HERE / "AX3T_L5.npz", allow_pickle=True)
a3b = np.load(HERE / "AX3T_L4.npz", allow_pickle=True)
a2 = np.load(HERE / "AXT_L6.npz", allow_pickle=True)
a2b = np.load(HERE / "AXT_L5.npz", allow_pickle=True)
C3DL, C2DL = "#e39c94", "#8fa8d1"       # light tints for second-finest

fig, ax = plt.subplots(figsize=(7.0, 2.9))
ax.set_axisbelow(True)
ax.grid(True, ls="--", lw=0.5, color="0.87")
def axline(d, xkey, rkey, **kw):
    x = d[xkey] / D
    r = d[rkey] / RHO_J
    m = (x > 0.05) & (x <= 10.0)
    ax.plot(x[m], boxfilt(r, x, 0.05)[m], **kw)
# second-finest grids, individually identified
axline(a3b, "x", "ax_DensityMEAN", color=C3DL, lw=0.9, zorder=3)
axline(a2b, "z", "stat_DensityMEAN", color=C2DL, lw=0.9, ls="--", zorder=3)
# finest grids
axline(a3, "x", "ax_DensityMEAN", color=C3D, lw=1.5, zorder=4)
axline(a2, "z", "stat_DensityMEAN", color=C2D, lw=1.2, ls="--", zorder=4)
# primary experimental reference: matched condition M_j = 1.42
ax.plot(xd["x"], xd["rho"], "o", color="k", ms=4.0, zorder=5)
ax.set_xlim(0, 10); ax.set_ylim(0.3, 1.55)
ax.set_xlabel("$x/D_e$", labelpad=2)
ax.set_ylabel("$\\langle\\rho\\rangle/\\rho_j$", labelpad=2)
ax.tick_params(length=2.5)
blank = Line2D([], [], ls="", marker="", label=" ")
handles = [
    Line2D([], [], marker="o", ls="", color="k", ms=4.0,
           label="$M_j=1.42$ experiment"),
    blank,
    Line2D([], [], color=C2D, lw=1.2, ls="--", label="RZ Level 6"),
    Line2D([], [], color=C2DL, lw=0.9, ls="--", label="RZ Level 5"),
    Line2D([], [], color=C3D, lw=1.5, label="Cartesian Level 5"),
    Line2D([], [], color=C3DL, lw=0.9, label="Cartesian Level 4"),
]
ax.legend(handles=handles, fontsize=6.8, frameon=False, ncol=3,
          loc="upper right", handlelength=1.6, columnspacing=1.2)
fig.tight_layout()
fig.savefig(HERE / "paperfig_cmp_axis_rho.png", dpi=300)
fig.savefig(HERE / "paperfig_cmp_axis_rho.pdf")
print("wrote paperfig_cmp_axis_rho")

# ---------- Fig 1b: near-field/downstream broken-axis comparison ----------
fig, (ax_left, ax_right) = plt.subplots(
    1, 2, figsize=(7.2, 3.2), sharey=True,
    gridspec_kw={"width_ratios": [3.15, 1.75], "wspace": 0.04},
)

def broken_line(axis, data, xkey, rkey, xlo, xhi, **kwargs):
    x = data[xkey] / D
    rho = data[rkey] / RHO_J
    mask = (x >= xlo) & (x <= xhi)
    axis.plot(x[mask], boxfilt(rho, x, 0.05)[mask], **kwargs)

for axis in (ax_left, ax_right):
    axis.set_axisbelow(True)
    axis.grid(True, ls="--", lw=0.5, color="0.87")
    broken_line(axis, a3b, "x", "ax_DensityMEAN",
                0.0 if axis is ax_left else 20.0,
                8.0 if axis is ax_left else 32.0,
                color=C3DL, lw=1.05, zorder=3)
    broken_line(axis, a2b, "z", "stat_DensityMEAN",
                0.0 if axis is ax_left else 20.0,
                8.0 if axis is ax_left else 32.0,
                color=C2DL, lw=1.05, ls="--", zorder=3)
    broken_line(axis, a3, "x", "ax_DensityMEAN",
                0.0 if axis is ax_left else 20.0,
                8.0 if axis is ax_left else 32.0,
                color=C3D, lw=1.6, zorder=4)
    broken_line(axis, a2, "z", "stat_DensityMEAN",
                0.0 if axis is ax_left else 20.0,
                8.0 if axis is ax_left else 32.0,
                color=C2D, lw=1.35, ls="--", zorder=4)
    axis.tick_params(length=3.0, labelsize=9.5)

ax_left.plot(xd["x"], xd["rho"], "o", color="k", ms=4.5, zorder=5)
ax_left.set_xlim(0, 8)
ax_right.set_xlim(20, 32)
ax_left.set_xticks(np.arange(0, 8, 1))
ax_right.set_xticks([20, 24, 28, 32])
ax_left.set_ylim(0.3, 1.55)
ax_left.set_ylabel(r"$\langle\rho\rangle/\rho_j$", labelpad=2,
                   fontsize=10.5)
fig.supxlabel(r"$x/D_e$", y=0.01, fontsize=10.5)

ax_left.spines["right"].set_visible(False)
ax_right.spines["left"].set_visible(False)
ax_left.tick_params(right=False)
ax_right.tick_params(left=False, labelleft=False)
break_size = 0.015
break_style = dict(color="k", clip_on=False, lw=1.0)
ax_left.plot((1-break_size, 1+break_size), (-break_size, +break_size),
             transform=ax_left.transAxes, **break_style)
ax_left.plot((1-break_size, 1+break_size), (1-break_size, 1+break_size),
             transform=ax_left.transAxes, **break_style)
ax_right.plot((-break_size, +break_size), (-break_size, +break_size),
              transform=ax_right.transAxes, **break_style)
ax_right.plot((-break_size, +break_size), (1-break_size, 1+break_size),
              transform=ax_right.transAxes, **break_style)

broken_handles = [
    Line2D([], [], marker="o", ls="", color="k", ms=4.5,
           label=r"$M_j=1.42$ experiment"),
    Line2D([], [], color=C3D, lw=1.6, label="Cartesian Level 5"),
    Line2D([], [], color=C3DL, lw=1.05, label="Cartesian Level 4"),
    Line2D([], [], color=C2D, lw=1.35, ls="--", label="RZ Level 6"),
    Line2D([], [], color=C2DL, lw=1.05, ls="--", label="RZ Level 5"),
]
ax_right.legend(
    handles=broken_handles, fontsize=8.0, frameon=False, ncol=1,
    loc="upper right", handlelength=1.7, labelspacing=0.35,
    borderaxespad=0.35,
)
fig.subplots_adjust(left=0.09, right=0.985, bottom=0.18, top=0.97)
fig.savefig(HERE / "paperfig_cmp_axis_rho_break.png", dpi=300)
fig.savefig(HERE / "paperfig_cmp_axis_rho_break.pdf")
print("wrote paperfig_cmp_axis_rho_break")

# ---------- Fig 2: radial profiles (Cartesian only, both cuts, full +-r) ----------
f5 = np.load(HERE.parent / "panda_fig5a_mj144_EPAPS.npz", allow_pickle=True)
RHOJ_EXP = float(f5["rhoj"])
fm = np.load(HERE.parent / "panda_m142den_field_EPAPS.npz")
xm, rm, Fm = fm["x"], fm["r"], fm["rho"]        # already rho/rho_j
sly = np.load(HERE / "SIM_L5_rhoMEAN_z0.npz")   # cut along y (z = 0 plane)
slz = np.load(HERE / "SIM_L5_rhoMEAN_y0.npz")   # cut along z (y = 0 plane)
xs3, rr3 = sly["x"], sly["y"]                   # identical grids (verified)
cutY = sly["rho"] / RHO_J
cutZ = slz["rho"] / RHO_J
C3Z = "#e0862c"                                 # z-cut colour

STATIONS = [0.6, 0.75, 0.9, 1.05, 1.25, 1.55, 2.0, 3.05]
OFF = 0.05                              # envelope half-width [D_e]

fig, axs = plt.subplots(2, 4, figsize=(7.2, 4.15), sharex=True, sharey=True)
for k, (stn, ax) in enumerate(zip(STATIONS, axs.flat)):
    ax.set_axisbelow(True)
    ax.grid(True, ls="--", lw=0.5, color="0.87")
    e = f5["x%s" % str(stn)]
    rex = e[:, 0]; rhex = e[:, 1] / RHOJ_EXP    # measured side (r <= 0)
    icol = np.argmin(np.abs(xm - stn))          # nearest measurement column
    xcol = xm[icol]; rho42 = Fm[icol]
    # min-max envelope of BOTH cuts over |x - stn| <= OFF, signed r
    jsel = np.abs(xs3 - stn) <= OFF + 1e-9
    stack = np.vstack([cutY[jsel], cutZ[jsel]])
    ax.fill_between(rr3, stack.min(0), stack.max(0), color=C3D, alpha=0.15,
                    lw=0, zorder=2)
    # nominal station: the two diametral cuts, unfolded; identical dash
    # period, half-period phase offset -> alternating red/blue segments
    # where the two curves coincide
    j = np.argmin(np.abs(xs3 - stn))
    ax.plot(rr3, cutY[j], color=C3D, lw=1.1, ls=(0, (4, 4)), zorder=4)
    ax.plot(rr3, cutZ[j], color=C2D, lw=1.1, ls=(4, (4, 4)), zorder=4)
    # Experiments at M_j=1.42 (two-sided, primary) and M_j=1.44
    # (measured on one side and mirrored, strongly de-emphasised)
    ax.plot(rex, rhex, "o", color="0.72", ms=2.5, mfc="none", mew=0.7,
            zorder=3)
    ax.plot(-rex, rhex, "o", color="0.72", ms=2.5, mfc="none", mew=0.7,
            zorder=3)
    ax.plot(rm, rho42, "o", color="k", ms=3.2, mew=0, zorder=5)
    ax.text(0.05, 0.06, "$x/D_e=%g$" % stn,
            transform=ax.transAxes, ha="left", va="bottom", fontsize=7.5,
            zorder=6,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.75,
                      pad=1.2))
    ax.set_xlim(-0.85, 0.85)
    ax.set_xticks([-0.5, 0.0, 0.5])
    ax.tick_params(length=2.5)
for ax in axs[1]:
    ax.set_xlabel("$y/D_e$", labelpad=2)
for ax in axs[:, 0]:
    ax.set_ylabel("$\\langle\\rho\\rangle/\\rho_j$", labelpad=2)
axs[0, 0].set_ylim(0.3, 1.5)
handles = [
    Line2D([], [], marker="o", ls="", color="k", ms=3.2,
           label="$M_j=1.42$ experiment"),
    Line2D([], [], marker="o", ls="", color="0.72", mfc="none", mew=0.7,
           ms=2.5, label="$M_j=1.44$ experiment (mirrored)"),
    Line2D([], [], color=C3D, lw=1.1, ls=(0, (4, 4)),
           label="$y$-cut"),
    Line2D([], [], color=C2D, lw=1.1, ls=(4, (4, 4)),
           label="$z$-cut"),
    Patch(facecolor=C3D, alpha=0.15, label="$x\\pm0.05\\,D_e$ envelope"),
]
fig.legend(handles=handles, fontsize=8.0, frameon=False, ncol=5,
           loc="upper center", bbox_to_anchor=(0.5, 1.0),
           handlelength=1.3, columnspacing=0.8, handletextpad=0.4)
fig.tight_layout(w_pad=0.5, h_pad=0.5, rect=[0, 0, 1, 0.955])
fig.savefig(HERE / "paperfig_cmp_radial_rho.png", dpi=300)
fig.savefig(HERE / "paperfig_cmp_radial_rho.pdf")
print("wrote paperfig_cmp_radial_rho")

# ---------- Fig 3: radial profiles, axisymmetric radial--axial version ----------
# Same composition as Fig 2 with the RZ Level 6 solution: one profile per
# station, mirrored about the axis (exact by axisymmetry); envelope from
# the stored axial-offset profiles within +-OFF.
rt = np.load(HERE / "RADT_L6.npz", allow_pickle=True)
st2 = list(rt["stations"])
offs2 = np.asarray(rt["offsets"], dtype=float)
i0 = list(offs2).index(0.0)
osel2 = np.abs(offs2) <= OFF + 1e-9
rr2 = rt["rr"]                                  # already r/D, one-sided
rfull = np.concatenate([-rr2[::-1], rr2])

def mirror(p):
    return np.concatenate([p[::-1], p])

fig, axs = plt.subplots(2, 4, figsize=(7.2, 4.15), sharex=True, sharey=True)
for k, (stn, ax) in enumerate(zip(STATIONS, axs.flat)):
    ax.set_axisbelow(True)
    ax.grid(True, ls="--", lw=0.5, color="0.87")
    e = f5["x%s" % str(stn)]
    rex = e[:, 0]; rhex = e[:, 1] / RHOJ_EXP
    icol = np.argmin(np.abs(xm - stn))
    js = st2.index(stn)
    profs = np.array([boxfilt(p / RHO_J, rr2, 0.02)
                      for p in rt["rad_DensityMEAN"][js][osel2]])
    ax.fill_between(rfull, mirror(profs.min(0)), mirror(profs.max(0)),
                    color=C2D, alpha=0.15, lw=0, zorder=2)
    p2 = boxfilt(rt["rad_DensityMEAN"][js, i0] / RHO_J, rr2, 0.02)
    ax.plot(rfull, mirror(p2), color=C2D, lw=1.1, ls="--", zorder=4)
    ax.plot(rex, rhex, "o", color="0.72", ms=2.5, mfc="none", mew=0.7,
            zorder=3)
    ax.plot(-rex, rhex, "o", color="0.72", ms=2.5, mfc="none", mew=0.7,
            zorder=3)
    ax.plot(rm, Fm[icol], "o", color="k", ms=3.2, mew=0, zorder=5)
    ax.text(0.05, 0.06, "$x/D_e=%g$" % stn,
            transform=ax.transAxes, ha="left", va="bottom", fontsize=7.5,
            zorder=6,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.75,
                      pad=1.2))
    ax.set_xlim(-0.85, 0.85)
    ax.set_xticks([-0.5, 0.0, 0.5])
    ax.tick_params(length=2.5)
for ax in axs[1]:
    ax.set_xlabel("$y/D_e$", labelpad=2)
for ax in axs[:, 0]:
    ax.set_ylabel("$\\langle\\rho\\rangle/\\rho_j$", labelpad=2)
axs[0, 0].set_ylim(0.3, 1.5)
handles = [
    Line2D([], [], marker="o", ls="", color="k", ms=3.2,
           label="$M_j=1.42$ experiment"),
    Line2D([], [], marker="o", ls="", color="0.72", mfc="none", mew=0.7,
           ms=2.5, label="$M_j=1.44$ experiment (mirrored)"),
    Line2D([], [], color=C2D, lw=1.1, ls="--",
           label="RZ Level 6"),
    Patch(facecolor=C2D, alpha=0.15, label="$x\\pm0.05\\,D_e$ envelope"),
]
fig.legend(handles=handles, fontsize=8.0, frameon=False, ncol=4,
           loc="upper center", bbox_to_anchor=(0.5, 1.0),
           handlelength=1.3, columnspacing=0.9, handletextpad=0.4)
fig.tight_layout(w_pad=0.5, h_pad=0.5, rect=[0, 0, 1, 0.955])
fig.savefig(HERE / "paperfig_cmp_radial_rho_2d.png", dpi=300)
fig.savefig(HERE / "paperfig_cmp_radial_rho_2d.pdf")
print("wrote paperfig_cmp_radial_rho_2d")
