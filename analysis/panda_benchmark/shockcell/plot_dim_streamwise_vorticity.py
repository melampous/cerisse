#!/usr/bin/env python3
"""Streamwise-vorticity cross-sections for Thesis Section 3.7.

The saved y-z planes are composite Level 5 data sampled on a uniform D_e/256
array.  Coarser AMR values are repeated on this fine array.  Each derivative
therefore uses a centred stride that matches the local AMR level.  Streamwise
vorticity is evaluated from the in-plane velocity components,

    omega_x = d(u_z)/dy - d(u_y)/dz.

The local AMR level is reconstructed from the saved BoxArray bounds.  A value
is masked if any point in its centred derivative stencil belongs to another
level.  Four representative stations are used.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import distance_transform_edt


HERE = Path(__file__).parent
D_E = 0.0254
U_J = 414.2
STATIONS = (0.60, 1.25, 2.00, 3.05)
R_MAX = 0.72

plt.rcParams.update({
    "font.size": 9,
    "axes.linewidth": 0.6,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
})


def station_key(station):
    return "x" + str(float(station)).replace(".", "p")


def finest_level_map(box_data, y_dim, z_dim, x_dim):
    """Return the finest saved AMR level covering each y-z sample point."""
    y_phys = y_dim * D_E
    z_phys = z_dim * D_E
    y_grid, z_grid = np.meshgrid(y_phys, z_phys, indexing="ij")
    level = np.full(y_grid.shape, -1, dtype=np.int8)
    x_phys = x_dim * D_E
    tol = 1.0e-12

    for lev in range(6):
        for box in box_data[f"L5_lev{lev}"]:
            if not (box[0] - tol <= x_phys <= box[1] + tol):
                continue
            inside = (
                (y_grid >= box[2] - tol)
                & (y_grid <= box[3] + tol)
                & (z_grid >= box[4] - tol)
                & (z_grid <= box[5] + tol)
            )
            level[inside] = np.maximum(level[inside], lev)

    if np.any(level < 0):
        raise RuntimeError(f"AMR coverage is incomplete at x/D_e={x_dim:g}")
    return level


def local_level_vorticity(u_y, u_z, level, dy, dz):
    """Evaluate omega_x with a centred stride matched to the local level."""
    omega = np.full(u_y.shape, np.nan, dtype=np.float64)
    valid_all = np.zeros(u_y.shape, dtype=bool)
    finest = 5

    for lev in np.unique(level):
        stride = 2 ** (finest - int(lev))
        centre = np.s_[stride:-stride, stride:-stride]
        local_level = level[centre]
        valid = local_level == lev
        valid &= level[:-2 * stride, stride:-stride] == lev
        valid &= level[2 * stride:, stride:-stride] == lev
        valid &= level[stride:-stride, :-2 * stride] == lev
        valid &= level[stride:-stride, 2 * stride:] == lev

        duz_dy = (
            u_z[2 * stride:, stride:-stride]
            - u_z[:-2 * stride, stride:-stride]
        ) / (2.0 * stride * dy)
        duy_dz = (
            u_y[stride:-stride, 2 * stride:]
            - u_y[stride:-stride, :-2 * stride]
        ) / (2.0 * stride * dz)
        local_omega = duz_dy - duy_dz

        omega_view = omega[centre]
        omega_view[valid] = local_omega[valid]
        omega[centre] = omega_view
        valid_view = valid_all[centre]
        valid_view[valid] = True
        valid_all[centre] = valid_view

    return omega, valid_all


def fill_interface_mask_for_display(field):
    """Fill narrow derivative-stencil gaps using the nearest valid value.

    The AMR-aware derivative remains unchanged.  This operation is applied
    only to the displayed array so that level interfaces do not dominate the
    qualitative cross-section comparison.
    """
    invalid = ~np.isfinite(field)
    if not np.any(invalid):
        return field
    nearest = distance_transform_edt(
        invalid,
        return_distances=False,
        return_indices=True,
    )
    return field[tuple(nearest)]


planes = np.load(HERE / "SIM_L5_yz8_multi.npz", allow_pickle=True)
boxes = np.load(HERE / "GRID_boxes_3d_sens.npz")
structure = np.load(HERE / "SIM_L5_yz8_struct.npz")
field_names = list(planes["want"])
y = planes["y"].astype(np.float64)
z = planes["z"].astype(np.float64)
dy = float(planes["dy"])
dz = float(planes["dz"])

iy = field_names.index("y_velocity")
iz = field_names.index("z_velocity")

omega_x = []
levels = []
masks = []
streamwise_fraction = []
for station in STATIONS:
    values = planes[station_key(station)]
    u_y = values[iy].astype(np.float64)
    u_z = values[iz].astype(np.float64)
    lev = finest_level_map(boxes, y, z, station)
    vort, valid = local_level_vorticity(u_y, u_z, lev, dy, dz)
    vort *= D_E / U_J
    suffix = str(float(station)).replace(".", "p")
    omega_magnitude = structure[f"om_{suffix}"].astype(np.float64)
    omega_magnitude *= D_E / U_J
    radial_mask = y[:, None] ** 2 + z[None, :] ** 2 <= 0.8 ** 2
    analysis_mask = valid & radial_mask & np.isfinite(omega_magnitude)
    denominator = np.sum(omega_magnitude[analysis_mask] ** 2)
    fraction = np.sum(vort[analysis_mask] ** 2) / denominator
    omega_x.append(vort)
    levels.append(lev)
    masks.append(~valid)
    streamwise_fraction.append(fraction)

plot_radius = y[:, None] ** 2 + z[None, :] ** 2 <= 0.8 ** 2
all_values = np.concatenate([
    np.abs(field[np.isfinite(field) & plot_radius]).ravel()
    for field in omega_x
])
limit = 5.0 * np.ceil(np.nanpercentile(all_values, 99.5) / 5.0)

cmap = plt.get_cmap("RdBu_r").copy()
fig, axes = plt.subplots(
    1, 4, figsize=(8.2, 2.25), sharex=True, sharey=True
)
theta = np.linspace(0.0, 2.0 * np.pi, 361)

image = None
for panel, (station, field, ax) in enumerate(
    zip(STATIONS, omega_x, axes)
):
    display_field = fill_interface_mask_for_display(field)
    image = ax.imshow(
        display_field,
        origin="lower",
        extent=[z[0], z[-1], y[0], y[-1]],
        vmin=-limit,
        vmax=limit,
        cmap=cmap,
        aspect="equal",
        interpolation="nearest",
        rasterized=True,
    )
    ax.plot(
        0.5 * np.cos(theta),
        0.5 * np.sin(theta),
        color="0.20",
        linewidth=0.65,
        linestyle=(0, (2.0, 1.8)),
    )
    ax.text(
        0.045,
        0.95,
        f"({chr(97 + panel)}) $x/D_e={station:g}$",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7.4,
        bbox={
            "facecolor": "white",
            "edgecolor": "none",
            "alpha": 0.78,
            "pad": 1.1,
        },
    )
    ax.set_xlim(-R_MAX, R_MAX)
    ax.set_ylim(-R_MAX, R_MAX)
    ax.set_xticks([-0.5, 0.0, 0.5])
    ax.set_yticks([-0.5, 0.0, 0.5])
    ax.tick_params(length=2.4, labelsize=7.7)

for ax in axes:
    ax.set_xlabel("$z/D_e$", labelpad=1.5)
axes[0].set_ylabel("$y/D_e$", labelpad=1.0)

fig.subplots_adjust(
    left=0.055,
    right=0.995,
    top=0.985,
    bottom=0.255,
    wspace=0.06,
)
color_axis = fig.add_axes([0.32, 0.065, 0.36, 0.035])
colorbar = fig.colorbar(image, cax=color_axis, orientation="horizontal")
colorbar.set_label("$\\omega_x D_e/U_j$", fontsize=8.2, labelpad=1.0)
colorbar.ax.xaxis.set_label_position("top")
colorbar.set_ticks([-25.0, -12.5, 0.0, 12.5, 25.0])
colorbar.ax.tick_params(length=2.2, labelsize=7.6)

output = HERE / "paperfig_dim_streamwise_vorticity"
fig.savefig(output.with_suffix(".png"), dpi=300)
fig.savefig(output.with_suffix(".pdf"), dpi=300)
print(f"wrote {output.name}.png/pdf")
print(f"common colour limit: +/-{limit:g}")
for station, level, mask in zip(STATIONS, levels, masks):
    available = ",".join(str(item) for item in np.unique(level))
    masked = 100.0 * np.count_nonzero(mask) / mask.size
    print(
        f"x/D_e={station:g}: levels {available}; "
        f"interface mask {masked:.2f}%"
    )
for station, fraction in zip(STATIONS, streamwise_fraction):
    print(f"x/D_e={station:g}: chi_x={fraction:.4f}")
