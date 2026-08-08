#!/usr/bin/env python3
"""Instantaneous cross-sectional vorticity-layer fragmentation diagnostics.

The analysis uses the t = 10.16 ms y-z planes from the 3D Cartesian L5
calculation.  All planes are area-averaged to D_e/128 before comparison so
that the five stations use the same sampling width.  The analysed field is

    W = |omega| D_e/U_j.

Three quantities are reported:

  F_NA      non-axisymmetric fraction of the azimuthal W spectrum;
  N_eff     area-weighted effective number of connected high-W regions;
  F_max     area fraction occupied by the largest retained region.

The connected-component reference threshold is 0.40 W_99 at each station.
The range produced by 0.35--0.45 W_99 is shown to expose threshold
sensitivity.  Components smaller than 25 cells on the common D_e/128 grid are
rejected.
"""

from pathlib import Path
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import label, map_coordinates


HERE = Path(__file__).parent
D_E = 0.0254
U_J = 414.2
SCALE_VORTICITY = D_E / U_J

STATIONS = (0.60, 0.90, 1.05, 1.25, 1.55)
DISPLAY_STATIONS = (0.60, 1.25, 1.55)
R_INNER = 0.30
R_OUTER = 0.80
M_MAX = 32
M_CUTOFFS = (24, 32, 48)
KAPPAS = (0.35, 0.40, 0.45)
REFERENCE_KAPPA = 0.40
MINIMUM_CELLS = 25

plt.rcParams.update({
    "font.size": 8.5,
    "axes.linewidth": 0.6,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "legend.frameon": False,
})


def station_suffix(station):
    return str(float(station)).replace(".", "p")


def area_average_2x2(field):
    """Area-average a D_e/256 plane onto the common D_e/128 grid."""
    ny = field.shape[0] // 2 * 2
    nz = field.shape[1] // 2 * 2
    return field[:ny, :nz].reshape(ny // 2, 2, nz // 2, 2).mean(axis=(1, 3))


def area_average_coordinates(coordinate):
    n = coordinate.size // 2 * 2
    return coordinate[:n].reshape(n // 2, 2).mean(axis=1)


def azimuthal_power(field, y, z):
    """Return the spectrum of the radially integrated angular signal."""
    radius = np.linspace(R_INNER, R_OUTER, 128)
    theta = np.linspace(0.0, 2.0 * np.pi, 512, endpoint=False)
    rr, tt = np.meshgrid(radius, theta, indexing="ij")
    yy = rr * np.cos(tt)
    zz = rr * np.sin(tt)
    coordinates = np.array([
        (yy - y[0]) / (y[1] - y[0]),
        (zz - z[0]) / (z[1] - z[0]),
    ])
    polar = map_coordinates(field, coordinates, order=1, mode="nearest")
    angular_signal = np.trapz(polar * radius[:, None], radius, axis=0)
    fourier = np.fft.rfft(angular_signal) / theta.size
    power = np.abs(fourier) ** 2
    power[1:-1] *= 2.0
    return power


def component_metrics(field, annulus, kappa, cell_width):
    """Return N_eff, F_max, d_eff/D_e, and the retained component count."""
    w99 = np.percentile(field[annulus], 99.0)
    binary = annulus & (field >= kappa * w99)
    component_labels, _ = label(binary, structure=np.ones((3, 3), dtype=int))
    areas = np.bincount(component_labels.ravel())[1:]
    areas = areas[areas >= MINIMUM_CELLS].astype(np.float64)
    if areas.size == 0:
        return np.nan, np.nan, np.nan, 0, w99
    total_area = np.sum(areas)
    effective_number = total_area ** 2 / np.sum(areas ** 2)
    largest_fraction = np.max(areas) / total_area
    effective_area = np.sum(areas ** 2) / total_area
    effective_diameter = 2.0 * np.sqrt(effective_area / np.pi) * cell_width
    return (
        effective_number,
        largest_fraction,
        effective_diameter,
        int(areas.size),
        w99,
    )


structure = np.load(HERE / "SIM_L5_yz8_struct.npz")
planes = np.load(HERE / "SIM_L5_yz8_multi.npz", allow_pickle=True)
y = area_average_coordinates(planes["y"].astype(np.float64))
z = area_average_coordinates(planes["z"].astype(np.float64))
dy = float(y[1] - y[0])
dz = float(z[1] - z[0])

if not np.isclose(dy, 1.0 / 128.0) or not np.isclose(dz, 1.0 / 128.0):
    raise RuntimeError(
        f"common cross-sectional spacing is ({dy:g}, {dz:g}), not D_e/128"
    )

radial_coordinate = np.sqrt(y[:, None] ** 2 + z[None, :] ** 2)
annulus = (
    (radial_coordinate >= R_INNER)
    & (radial_coordinate <= R_OUTER)
)

fields = {}
nonaxisymmetric_fraction = {cutoff: [] for cutoff in M_CUTOFFS}
effective_number = {kappa: [] for kappa in KAPPAS}
largest_fraction = {kappa: [] for kappa in KAPPAS}
effective_diameter = {kappa: [] for kappa in KAPPAS}
retained_count = {kappa: [] for kappa in KAPPAS}
w99_by_station = {}

for station in STATIONS:
    suffix = station_suffix(station)
    vorticity = structure[f"om_{suffix}"].astype(np.float64)
    w_field = area_average_2x2(vorticity * SCALE_VORTICITY)
    fields[station] = w_field

    power = azimuthal_power(w_field, y, z)
    for cutoff in M_CUTOFFS:
        retained_power = power[:cutoff + 1]
        fraction = 1.0 - retained_power[0] / np.sum(retained_power)
        nonaxisymmetric_fraction[cutoff].append(float(fraction))

    for kappa in KAPPAS:
        n_eff, f_max, d_eff, count, w99 = component_metrics(
            w_field,
            annulus,
            kappa,
            dy,
        )
        effective_number[kappa].append(n_eff)
        largest_fraction[kappa].append(f_max)
        effective_diameter[kappa].append(d_eff)
        retained_count[kappa].append(count)
        w99_by_station[station] = w99

for cutoff in M_CUTOFFS:
    nonaxisymmetric_fraction[cutoff] = np.asarray(
        nonaxisymmetric_fraction[cutoff]
    )
for kappa in KAPPAS:
    effective_number[kappa] = np.asarray(effective_number[kappa])
    largest_fraction[kappa] = np.asarray(largest_fraction[kappa])
    effective_diameter[kappa] = np.asarray(effective_diameter[kappa])
    retained_count[kappa] = np.asarray(retained_count[kappa])


def sensitivity_range(values, keys, reference_key):
    stack = np.vstack([values[key] for key in keys])
    reference = values[reference_key]
    lower = reference - np.nanmin(stack, axis=0)
    upper = np.nanmax(stack, axis=0) - reference
    return reference, np.vstack((lower, upper))


f_na_reference, f_na_error = sensitivity_range(
    nonaxisymmetric_fraction,
    M_CUTOFFS,
    M_MAX,
)
n_eff_reference, n_eff_error = sensitivity_range(
    effective_number,
    KAPPAS,
    REFERENCE_KAPPA,
)
f_max_reference, f_max_error = sensitivity_range(
    largest_fraction,
    KAPPAS,
    REFERENCE_KAPPA,
)

# Save the values used in the figure for direct audit and later text updates.
csv_path = HERE / "instantaneous_vortex_fragmentation_metrics.csv"
with csv_path.open("w", newline="", encoding="utf-8") as stream:
    writer = csv.writer(stream)
    writer.writerow([
        "x_over_De",
        "F_NA_m24",
        "F_NA_m32",
        "F_NA_m48",
        "N_eff_kappa035",
        "N_eff_kappa040",
        "N_eff_kappa045",
        "F_max_kappa035",
        "F_max_kappa040",
        "F_max_kappa045",
        "d_eff_over_De_kappa040",
        "N_components_kappa040",
    ])
    for index, station in enumerate(STATIONS):
        writer.writerow([
            f"{station:.3f}",
            f"{nonaxisymmetric_fraction[24][index]:.8f}",
            f"{nonaxisymmetric_fraction[32][index]:.8f}",
            f"{nonaxisymmetric_fraction[48][index]:.8f}",
            f"{effective_number[0.35][index]:.8f}",
            f"{effective_number[0.40][index]:.8f}",
            f"{effective_number[0.45][index]:.8f}",
            f"{largest_fraction[0.35][index]:.8f}",
            f"{largest_fraction[0.40][index]:.8f}",
            f"{largest_fraction[0.45][index]:.8f}",
            f"{effective_diameter[0.40][index]:.8f}",
            int(retained_count[0.40][index]),
        ])

fig, axes = plt.subplots(2, 3, figsize=(7.25, 4.65))
image = None
theta_lip = np.linspace(0.0, 2.0 * np.pi, 361)

for panel, (station, axis) in enumerate(zip(DISPLAY_STATIONS, axes[0])):
    field = fields[station]
    image = axis.imshow(
        field,
        origin="lower",
        extent=[z[0], z[-1], y[0], y[-1]],
        vmin=0.0,
        vmax=40.0,
        cmap="magma",
        interpolation="nearest",
        rasterized=True,
        aspect="equal",
    )
    axis.contour(
        z,
        y,
        field,
        levels=[REFERENCE_KAPPA * w99_by_station[station]],
        colors="#65f0c0",
        linewidths=0.62,
    )
    axis.plot(
        0.5 * np.cos(theta_lip),
        0.5 * np.sin(theta_lip),
        color="0.78",
        linewidth=0.55,
        linestyle=":",
    )
    axis.text(
        0.045,
        0.95,
        f"({chr(97 + panel)}) $x/D_e={station:g}$",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=7.7,
        bbox={
            "facecolor": "white",
            "edgecolor": "none",
            "alpha": 0.80,
            "pad": 1.1,
        },
    )
    axis.set_xlim(-0.82, 0.82)
    axis.set_ylim(-0.82, 0.82)
    axis.set_xticks([-0.5, 0.0, 0.5])
    axis.set_yticks([-0.5, 0.0, 0.5])
    axis.tick_params(length=2.4, labelsize=7.7)
    axis.set_xlabel("$z/D_e$", labelpad=1.3)

axes[0, 0].set_ylabel("$y/D_e$", labelpad=1.0)

blue = "#2864b0"
orange = "#d55e00"
green = "#18845b"

axis = axes[1, 0]
axis.errorbar(
    STATIONS,
    f_na_reference,
    yerr=f_na_error,
    color=blue,
    marker="o",
    markersize=4.1,
    markeredgecolor="white",
    markeredgewidth=0.5,
    linewidth=1.25,
    capsize=2.1,
    elinewidth=0.75,
)
axis.set_ylabel("$F_{\\mathrm{NA}}$", labelpad=2.0)
axis.set_ylim(0.0, 0.24)
axis.set_title("(d) Non-axisymmetric fraction", fontsize=8.1, pad=4)

axis = axes[1, 1]
axis.errorbar(
    STATIONS,
    n_eff_reference,
    yerr=n_eff_error,
    color=orange,
    marker="s",
    markersize=4.0,
    markeredgecolor="white",
    markeredgewidth=0.5,
    linewidth=1.25,
    capsize=2.1,
    elinewidth=0.75,
)
axis.set_ylabel("$N_{\\mathrm{eff}}$", labelpad=2.0)
axis.set_ylim(0.0, max(34.0, np.nanmax(n_eff_reference + n_eff_error[1]) * 1.08))
axis.set_title("(e) Effective component number", fontsize=8.1, pad=4)

axis = axes[1, 2]
axis.errorbar(
    STATIONS,
    f_max_reference,
    yerr=f_max_error,
    color=green,
    marker="^",
    markersize=4.3,
    markeredgecolor="white",
    markeredgewidth=0.5,
    linewidth=1.25,
    capsize=2.1,
    elinewidth=0.75,
)
axis.set_ylabel("$F_{\\max}$", labelpad=2.0)
axis.set_ylim(0.0, 1.05)
axis.set_title("(f) Largest-component fraction", fontsize=8.1, pad=4)

for axis in axes[1]:
    axis.set_xlabel("$x/D_e$", labelpad=1.5)
    axis.set_xlim(0.55, 1.60)
    axis.set_xticks([0.6, 0.9, 1.2, 1.5])
    axis.grid(True, color="0.88", linestyle="--", linewidth=0.48)
    axis.tick_params(length=2.4, labelsize=7.7)
    axis.set_axisbelow(True)

fig.subplots_adjust(
    left=0.075,
    right=0.885,
    top=0.965,
    bottom=0.105,
    wspace=0.34,
    hspace=0.34,
)
cbar_axis = fig.add_axes([0.905, 0.605, 0.018, 0.285])
colourbar = fig.colorbar(image, cax=cbar_axis)
colourbar.set_label(
    r"$W=|\mathbf{\omega}|D_e/U_j$",
    fontsize=8.1,
    labelpad=3,
)
colourbar.ax.tick_params(length=2.3, labelsize=7.5)

output_stem = HERE / "paperfig_instantaneous_vortex_fragmentation"
fig.savefig(output_stem.with_suffix(".png"), dpi=300)
fig.savefig(output_stem.with_suffix(".pdf"), dpi=300)

print(f"wrote {output_stem.name}.png/pdf")
print(f"wrote {csv_path.name}")
for index, station in enumerate(STATIONS):
    print(
        f"x/D_e={station:4.2f}: F_NA={f_na_reference[index]:.3f}, "
        f"N_eff={n_eff_reference[index]:.2f}, "
        f"F_max={f_max_reference[index]:.3f}, "
        f"d_eff/D_e={effective_diameter[REFERENCE_KAPPA][index]:.3f}, "
        f"N={retained_count[REFERENCE_KAPPA][index]:d}"
    )
