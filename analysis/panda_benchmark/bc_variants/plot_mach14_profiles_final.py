#!/usr/bin/env python3
"""Plot the mean-field Mach number for the four nozzle-exit profiles.

The left column uses the nozzle diameter and the right column uses the
mass-flow-equivalent exit diameter.  Both coordinates are rescaled in the
right column.  Corresponding centreline Mach maxima are linked between
adjacent rows to show the accumulated phase shift.

The downstream extent of the mean supersonic core is obtained from the
four-connected component of M_bar >= 1 that intersects the nozzle exit
(the first available x plane within |y| <= 0.5 D_e).  Its terminal x
coordinate is refined by linear interpolation to M_bar = 1 across the
downstream-facing edge of the component.  This definition retains off-axis
extensions of the connected core that a centreline-only criterion misses.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import ConnectionPatch
from scipy import ndimage
from scipy.signal import find_peaks


HERE = Path(__file__).resolve().parent
X_MAX = 14.0

# Values used in Table 3.2 of the thesis.
CD_EFF = {
    "baseline": 0.8789,
    "tophat": 1.0000,
    "walltanh": 0.9755,
    "shift100": 0.9513,
}

CASES = [
    ("baseline", "Baseline shifted tanh"),
    ("tophat", "Top hat"),
    ("walltanh", "Wall-attached tanh"),
    ("shift100", "Thin shifted tanh"),
]
plt.rcParams.update(
    {
        "font.size": 13.2,
        "axes.labelsize": 13.2,
        "axes.titlesize": 14.0,
        "xtick.labelsize": 11.5,
        "ytick.labelsize": 11.5,
        "axes.linewidth": 0.7,
        "xtick.direction": "in",
        "ytick.direction": "in",
    }
)


def centreline_mach(data):
    """Return x and the centreline Mach number from the two nearest y cells."""
    x, y, mach = data["x"], data["y"], data["M"]
    j_centre = np.argsort(np.abs(y))[:2]
    return x, mach[:, j_centre].mean(axis=1)


def mach_maxima(x, mach, n_max=14):
    """Identify the successive resolved centreline Mach-number maxima."""
    width = max(3, int(0.06 / (x[1] - x[0])))
    filtered = np.convolve(mach, np.ones(width) / width, mode="same")
    distance = int(0.5 / (x[1] - x[0]))

    near, _ = find_peaks(filtered, distance=distance, prominence=0.03)
    near = [index for index in near if 0.2 < x[index] <= 8.0]

    far, _ = find_peaks(filtered, distance=distance, prominence=0.01)
    last_near = x[near[-1]] if near else 0.2
    far = [index for index in far if last_near + 0.4 < x[index] < X_MAX + 0.5]

    return [float(x[index]) for index in (near + far)[:n_max]]


def connected_supersonic_core(field, nozzle_radius=0.5):
    """Return the nozzle-connected M_bar >= 1 mask and its terminal point.

    Four-neighbour connectivity is used so that two regions touching only at
    a grid-cell corner are not treated as physically connected.  The nozzle
    seed is the supersonic part of the first available axial plane within the
    geometric exit radius.  The terminal x coordinate is interpolated across
    all downstream-facing component edges, and the farthest crossing is
    retained.
    """
    x, y, mach = field["x"], field["y"], field["M"]
    supersonic = np.isfinite(mach) & (mach >= 1.0)
    structure = ndimage.generate_binary_structure(rank=2, connectivity=1)
    labels, _ = ndimage.label(supersonic, structure=structure)

    exit_seed = (np.abs(y) <= nozzle_radius) & supersonic[0, :]
    seed_labels = np.unique(labels[0, exit_seed])
    seed_labels = seed_labels[seed_labels > 0]
    if seed_labels.size == 0:
        raise RuntimeError("No supersonic component intersects the nozzle exit")

    connected = np.isin(labels, seed_labels)
    i_connected, j_connected = np.where(connected)
    if i_connected.size == 0:
        raise RuntimeError("The nozzle-connected supersonic component is empty")

    last = np.argmax(x[i_connected])
    x_end = float(x[i_connected[last]])
    y_end = float(y[j_connected[last]])

    # Refine the downstream M_bar=1 boundary from cell-centred values.
    for i, j in zip(i_connected, j_connected):
        if i + 1 >= x.size or connected[i + 1, j]:
            continue
        if not (mach[i, j] >= 1.0 and mach[i + 1, j] < 1.0):
            continue
        fraction = (mach[i, j] - 1.0) / (mach[i, j] - mach[i + 1, j])
        crossing = float(x[i] + fraction * (x[i + 1] - x[i]))
        if crossing > x_end:
            x_end = crossing
            y_end = float(y[j])

    return connected, x_end, y_end


data = {}
for key, _ in CASES:
    field = np.load(HERE / f"MACH16_{key}.npz")
    x_axis, mach_axis = centreline_mach(field)
    core_mask, core_x, core_y = connected_supersonic_core(field)
    data[key] = {
        "field": field,
        "maxima": mach_maxima(x_axis, mach_axis),
        "axis_mach": mach_axis,
        "core_mask": core_mask,
        "core_x": core_x,
        "core_y": core_y,
    }

n_peaks = min(len(data[key]["maxima"]) for key, _ in CASES)
peak_colours = plt.get_cmap("cool")

# A compact source canvas avoids excessive down-scaling when inserted at
# text width and keeps the final tick and annotation sizes near 7 pt.
fig, axes = plt.subplots(
    4,
    2,
    figsize=(10.2, 4.05),
    sharex="col",
    sharey=True,
)
image = None

for column, rescaled in enumerate((False, True)):
    for row, (key, label) in enumerate(CASES):
        ax = axes[row, column]
        field = data[key]["field"]
        x, y, mach = field["x"], field["y"], field["M"]
        scale = 1.0 / np.sqrt(CD_EFF[key]) if rescaled else 1.0
        x_plot = x * scale
        y_plot = y * scale

        image = ax.imshow(
            mach.T,
            origin="lower",
            extent=[x_plot[0], x_plot[-1], y_plot[0], y_plot[-1]],
            vmin=0.0,
            vmax=2.0,
            cmap="turbo",
            aspect="equal",
            interpolation="nearest",
            rasterized=True,
        )
        ax.contour(
            x_plot,
            y_plot,
            mach.T,
            levels=[1.0],
            colors="white",
            linewidths=0.6,
            alpha=0.9,
        )

        for index in range(n_peaks):
            x_peak = data[key]["maxima"][index] * scale
            ax.plot(
                x_peak,
                0.0,
                marker="v",
                linestyle="none",
                color=peak_colours(index / max(n_peaks - 1, 1)),
                markersize=5.8,
                markeredgecolor="black",
                markeredgewidth=0.55,
                zorder=6,
            )

        x_core = data[key]["core_x"] * scale
        y_core = data[key]["core_y"] * scale
        ax.plot(
            x_core,
            y_core,
            marker="*",
            linestyle="none",
            color="white",
            markersize=12,
            markeredgecolor="black",
            markeredgewidth=0.65,
            zorder=7,
        )
        ax.axvline(x_core, color="white", linewidth=0.75, linestyle="--")
        ax.text(
            x_core + 0.14,
            0.98,
            f"{x_core:.1f}",
            fontsize=11.0,
            color="white",
            va="top",
            bbox={
                "facecolor": "0.15",
                "edgecolor": "none",
                "alpha": 0.60,
                "boxstyle": "round,pad=0.16",
            },
        )

        ax.set_xlim(0.0, X_MAX)
        ax.set_ylim(-1.2, 1.2)
        ax.set_yticks((-1, 0, 1))
        ax.tick_params(
            length=3.0,
            width=0.7,
            top=False,
            right=False,
            labelbottom=(row == len(CASES) - 1),
        )
        if column == 0:
            ax.set_ylabel(r"$r/D_e$", fontsize=10.5, labelpad=4)
        if row == len(CASES) - 1:
            ax.set_xlabel(
                r"$x/D_e$" if not rescaled else r"$x/D_{\mathrm{eff},k}$"
            )

    # Link the same centreline maximum between adjacent profile rows.
    for index in range(n_peaks):
        colour = peak_colours(index / max(n_peaks - 1, 1))
        for row in range(len(CASES) - 1):
            first = CASES[row][0]
            second = CASES[row + 1][0]
            first_scale = 1.0 / np.sqrt(CD_EFF[first]) if rescaled else 1.0
            second_scale = 1.0 / np.sqrt(CD_EFF[second]) if rescaled else 1.0
            first_peak = data[first]["maxima"][index] * first_scale
            second_peak = data[second]["maxima"][index] * second_scale
            if first_peak > X_MAX or second_peak > X_MAX:
                continue
            fig.add_artist(
                ConnectionPatch(
                    xyA=(first_peak, -1.2),
                    coordsA=axes[row, column].transData,
                    xyB=(second_peak, 1.2),
                    coordsB=axes[row + 1, column].transData,
                    color=colour,
                    linewidth=1.15,
                    alpha=0.9,
                    zorder=10,
                    clip_on=False,
                )
            )

axes[0, 0].set_title(r"Original coordinate, $x/D_e$")
axes[0, 1].set_title(
    r"Rescaled coordinates, $(x,r)/D_{\mathrm{eff},k}$"
)

colourbar_axis = fig.add_axes([0.924, 0.255, 0.015, 0.49])
colourbar = fig.colorbar(image, cax=colourbar_axis)
colourbar.set_label(r"$M_{\overline{\mathbf{q}}}$", labelpad=5)
colourbar.ax.tick_params(labelsize=11.5, length=2.5)

fig.subplots_adjust(
    left=0.120,
    right=0.905,
    top=0.940,
    bottom=0.090,
    hspace=0.03,
    wspace=0.11,
)

for extension in ("png", "pdf"):
    fig.savefig(
        HERE / f"mach14_profiles_final.{extension}",
        dpi=240,
        bbox_inches="tight",
        pad_inches=0.03,
    )

print("saved mach14_profiles_final.png/pdf")
for key, label in CASES:
    print(
        f"{label:22s} maxima: "
        + " ".join(f"{value:.2f}" for value in data[key]["maxima"][:n_peaks])
    )
    scale = 1.0 / np.sqrt(CD_EFF[key])
    print(
        f"{'':22s} connected core: "
        f"x/D_e={data[key]['core_x']:.4f}, "
        f"r/D_e={data[key]['core_y']:.4f}, "
        f"x/D_eff,k={data[key]['core_x'] * scale:.4f}"
    )
