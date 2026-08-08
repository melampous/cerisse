#!/usr/bin/env python3
"""Figure 3.15: centreline mean density for four nozzle-exit profiles.

The numerical curves and the mean-density profiles used for feature extraction
are smoothed over the common Chapter 3 width of 0.05 D_e.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy.signal import find_peaks


HERE = Path(__file__).resolve().parent
RHO_J = 1.6413
X_LIMITS = (0.72, 7.25)

STYLES = {
    "baseline": {
        "label": "Baseline shifted tanh",
        "colour": "#aab0b7",
        "linestyle": "-",
        "linewidth": 2.2,
    },
    "tophat": {
        "label": "Top hat",
        "colour": "#0b6ef5",
        "linestyle": (0, (5.5, 2.2)),
        "linewidth": 1.7,
    },
    "walltanh": {
        "label": "Wall-attached tanh",
        "colour": "#2ca02c",
        "linestyle": (0, (3.4, 1.35, 1.0, 1.35)),
        "linewidth": 1.7,
    },
    "shift100": {
        "label": "Thin shifted tanh",
        "colour": "#ff6a00",
        "linestyle": (0, (1.0, 1.35)),
        "linewidth": 1.8,
    },
}


def box_filter(values, x, width):
    """Moving box filter of a requested dimensional width in D_e."""
    count = max(3, int(round(width / (x[1] - x[0]))))
    if count % 2 == 0:
        count += 1
    result = np.convolve(values, np.ones(count) / count, mode="same")
    half = count // 2
    result[:half] = values[:half]
    result[-half:] = values[-half:]
    return result


def quadratic_peak(x, values, index):
    """Return the location and value of a three-point quadratic maximum."""
    if index <= 0 or index >= values.size - 1:
        return float(x[index]), float(values[index])
    curvature = values[index - 1] - 2.0 * values[index] + values[index + 1]
    if curvature == 0.0:
        return float(x[index]), float(values[index])
    offset = 0.5 * (values[index - 1] - values[index + 1]) / curvature
    offset = float(np.clip(offset, -1.0, 1.0))
    dx = x[1] - x[0]
    location = x[index] + offset * dx
    value = values[index] - 0.25 * (
        values[index - 1] - values[index + 1]
    ) * offset
    return float(location), float(value)


def closure_location(x, density):
    """Density-gradient first-cell closure measure used in Chapter 3."""
    expansion = np.flatnonzero((x >= 0.35) & (x <= 1.0))
    first_minimum = expansion[np.argmin(density[expansion])]
    recompression = np.flatnonzero(
        (np.arange(x.size) > first_minimum) & (x <= 1.7)
    )
    first_maximum = recompression[np.argmax(density[recompression])]
    gradient = np.gradient(density, x)
    candidates = np.arange(first_minimum, first_maximum + 1)
    selected = candidates[np.argmax(gradient[candidates])]
    return quadratic_peak(x, gradient, selected)[0]


def density_maxima(x, density, x_fc, count=6):
    """Successive mean-density maxima downstream of x_fc."""
    minimum_distance = max(1, int(round(0.5 / (x[1] - x[0]))))
    indices, _ = find_peaks(
        np.where(x > x_fc, density, -np.inf),
        prominence=0.02,
        distance=minimum_distance,
    )
    peaks = [quadratic_peak(x, density, i) for i in indices[:count]]
    return np.asarray(peaks, dtype=float)


def load_simulation(case):
    data = np.load(HERE / f"AXIS8_{case}.npz")
    x = data["x"]
    raw = data["DensityMEAN"] / RHO_J
    displayed = box_filter(raw, x, width=0.05)
    feature_curve = displayed
    x_fc = closure_location(x, raw)
    peaks = density_maxima(x, feature_curve, x_fc)
    # Plot peak amplitudes on the displayed curve while retaining the
    # density-based feature locations selected from the feature curve.
    peaks[:, 1] = np.interp(peaks[:, 0], x, displayed)
    return x, displayed, peaks


def load_experiment():
    data = np.load(HERE.parent / "panda_m142den_axis_EPAPS.npz")
    x = data["x"]
    density = data["rho"]
    x_fc = closure_location(x, density)
    peaks = density_maxima(x, density, x_fc)
    return x, density, peaks


def format_axis(ax):
    ax.set_axisbelow(True)
    ax.grid(True, color="0.86", linewidth=0.55, linestyle="--")
    ax.tick_params(direction="in", top=True, right=True, length=3.2, width=0.75)


plt.rcParams.update(
    {
        "font.size": 11.0,
        "axes.labelsize": 11.5,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "axes.linewidth": 0.75,
        "legend.fontsize": 9.2,
        "mathtext.fontset": "dejavusans",
    }
)

experiment_x, experiment_density, experiment_peaks = load_experiment()
simulations = {case: load_simulation(case) for case in STYLES}

fig = plt.figure(figsize=(9.0, 6.15))
grid = fig.add_gridspec(
    2,
    2,
    height_ratios=(1.85, 1.0),
    hspace=0.37,
    wspace=0.28,
    left=0.095,
    right=0.985,
    bottom=0.105,
    top=0.885,
)
ax_a = fig.add_subplot(grid[0, :])
ax_b = fig.add_subplot(grid[1, 0])
ax_c = fig.add_subplot(grid[1, 1])
for axis in (ax_a, ax_b, ax_c):
    format_axis(axis)

# Panel (a): centreline profiles and density-based feature markers.
ax_a.plot(
    experiment_x,
    experiment_density,
    linestyle="none",
    marker="o",
    markersize=4.0,
    markerfacecolor="white",
    markeredgecolor="black",
    markeredgewidth=0.85,
    zorder=5,
)
ax_a.plot(
    experiment_peaks[:, 0],
    experiment_peaks[:, 1],
    linestyle="none",
    marker="o",
    markersize=5.2,
    markerfacecolor="black",
    markeredgecolor="black",
    zorder=7,
)

for case, style in STYLES.items():
    x, density, peaks = simulations[case]
    mask = (x >= X_LIMITS[0]) & (x <= X_LIMITS[1])
    ax_a.plot(
        x[mask],
        density[mask],
        color=style["colour"],
        linestyle=style["linestyle"],
        linewidth=style["linewidth"],
        zorder=3,
    )
    ax_a.plot(
        peaks[:, 0],
        peaks[:, 1],
        linestyle="none",
        marker="v",
        markersize=4.6,
        markerfacecolor=style["colour"],
        markeredgecolor="black",
        markeredgewidth=0.45,
        zorder=6,
    )

for index, location in enumerate(experiment_peaks[:, 0], start=1):
    ax_a.text(
        location,
        1.455,
        rf"$n={index}$",
        ha="center",
        va="bottom",
        fontsize=9.2,
        color="0.32",
    )

ax_a.set_xlim(*X_LIMITS)
ax_a.set_ylim(0.35, 1.52)
ax_a.set_xlabel(r"$x/D_e$")
ax_a.set_ylabel(r"$\langle\rho\rangle/\rho_j$ on the axis")
ax_a.text(0.012, 0.955, "(a)", transform=ax_a.transAxes, va="top", fontsize=11)

# Panel (b): amplitudes at the successive mean-density maxima.
indices = np.arange(1, experiment_peaks.shape[0] + 1)
ax_b.plot(
    indices,
    experiment_peaks[:, 1],
    color="black",
    linewidth=1.25,
    marker="o",
    markersize=4.7,
    markerfacecolor="white",
    markeredgewidth=0.9,
    zorder=5,
)
for case, style in STYLES.items():
    peaks = simulations[case][2]
    ax_b.plot(
        indices[: peaks.shape[0]],
        peaks[:, 1],
        color=style["colour"],
        linestyle=style["linestyle"],
        linewidth=style["linewidth"],
        marker="s",
        markersize=4.0,
        markerfacecolor=style["colour"],
        markeredgewidth=0.0,
        zorder=4,
    )
ax_b.set_xlim(0.75, 6.25)
ax_b.set_ylim(0.84, 1.42)
ax_b.set_xticks(indices)
ax_b.set_xlabel(r"index $n$ of successive mean-density maxima")
ax_b.set_ylabel(r"$\langle\rho\rangle_{\max,n}/\rho_j$")
ax_b.text(0.025, 0.08, "(b)", transform=ax_b.transAxes, fontsize=11)

# Panel (c): axial offsets of the computed maxima from the experiment.
ax_c.axhline(0.0, color="0.35", linewidth=0.8, zorder=2)
for case, style in STYLES.items():
    peaks = simulations[case][2]
    count = min(peaks.shape[0], experiment_peaks.shape[0])
    offset = peaks[:count, 0] - experiment_peaks[:count, 0]
    ax_c.plot(
        indices[:count],
        offset,
        color=style["colour"],
        linestyle=style["linestyle"],
        linewidth=style["linewidth"],
        marker="s",
        markersize=4.0,
        markerfacecolor=style["colour"],
        markeredgewidth=0.0,
        zorder=4,
    )
ax_c.set_xlim(0.75, 6.25)
ax_c.set_ylim(-0.33, 0.54)
ax_c.set_xticks(indices)
ax_c.set_xlabel(r"index $n$ of successive mean-density maxima")
ax_c.set_ylabel(r"$(x_n-x_n^{\mathrm{exp}})/D_e$")
ax_c.text(0.025, 0.94, "(c)", transform=ax_c.transAxes, va="top", fontsize=11)

legend_handles = [
    Line2D(
        [],
        [],
        linestyle="none",
        marker="o",
        markersize=4.5,
        markerfacecolor="white",
        markeredgecolor="black",
        label=r"Experiment ($M_j=1.42$)",
    )
]
legend_handles.extend(
    Line2D(
        [],
        [],
        color=style["colour"],
        linestyle=style["linestyle"],
        linewidth=style["linewidth"],
        label=style["label"],
    )
    for style in STYLES.values()
)
fig.legend(
    handles=legend_handles,
    loc="upper center",
    bbox_to_anchor=(0.54, 0.985),
    ncol=5,
    frameon=False,
    handlelength=2.5,
    handletextpad=0.55,
    columnspacing=1.35,
)

for extension in ("pdf", "png"):
    fig.savefig(
        HERE / f"axis_4profile_v2.{extension}",
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.04,
    )

print("Density maxima used in Figure 3.15")
print("Experiment:", np.round(experiment_peaks[:, 0], 5))
for case, style in STYLES.items():
    print(f"{style['label']:<22s}", np.round(simulations[case][2][:, 0], 5))
print("wrote axis_4profile_v2.pdf and axis_4profile_v2.png")
