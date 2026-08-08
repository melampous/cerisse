#!/usr/bin/env python3
"""Reproducible audit of the current Mj=1.42 Panda comparison data.

This script uses only the local NPZ files in this directory.  It deliberately
separates numerical-scheme differences, finite-window variability, spatial
dimensionality, and the near-axis density deficit.  It does not access AWS or
modify any simulation output.
"""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "current_data_analysis"
OUT.mkdir(exist_ok=True)

D = 0.0254
RHO_J_SIM = 1.6413
P_AMB = 99_780.0
R_AIR = 287.05
GAMMA = 1.4
WENO_STATS_START = 6.5e-3
TENO_STATS_START = 7.2e-3
THREED_STATS_START = 6.5e-3
THREED_STATS_END = 10.16e-3

RADIAL_STATIONS_8 = [0.60, 0.75, 0.90, 1.05, 1.25, 1.55, 2.00, 3.05]
RADIAL_STATIONS_30 = np.array(
    [
        0.60,
        0.75,
        0.90,
        1.05,
        1.25,
        1.40,
        1.55,
        1.70,
        1.85,
        2.00,
        2.15,
        2.30,
        2.45,
        2.60,
        2.75,
        2.90,
        3.05,
        3.20,
        3.35,
        3.50,
        3.75,
        4.00,
        4.25,
        4.50,
        4.75,
        5.00,
        5.25,
        5.50,
        5.75,
        6.00,
    ]
)

# The intervals are bounded by midpoints between the same-Mach experimental
# compression maxima.  This avoids treating the nozzle-exit shoulder near
# x/D=0.4 as a shock-cell maximum.
PEAK_INTERVALS = [
    (0.95, 1.865),
    (1.865, 3.14),
    (3.14, 4.30),
    (4.30, 5.34),
    (5.34, 6.30),
]


@dataclass
class Extremum:
    case: str
    cell: int
    peak_x: float
    peak_value: float
    preceding_trough_x: float
    preceding_trough_value: float
    peak_to_trough: float
    prominence: float
    spacing_from_previous: float | None
    resolved: bool


def load(name: str):
    return np.load(ROOT / name, allow_pickle=True)


def rms(a: np.ndarray, b: np.ndarray) -> float:
    valid = np.isfinite(a) & np.isfinite(b)
    return float(np.sqrt(np.mean((a[valid] - b[valid]) ** 2)))


def block_mean(
    mean_a: np.ndarray,
    time_a: float,
    mean_b: np.ndarray,
    time_b: float,
    stats_start: float = WENO_STATS_START,
) -> np.ndarray:
    """Recover the mean over [time_a, time_b] from two cumulative means."""

    if not (stats_start < time_a < time_b):
        raise ValueError("invalid cumulative-mean times")
    weight_a = time_a - stats_start
    weight_b = time_b - stats_start
    return (mean_b * weight_b - mean_a * weight_a) / (time_b - time_a)


def smooth_profile(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    grid = np.arange(0.50, 6.501, 0.005)
    values = np.interp(grid, x, y)
    values = savgol_filter(values, window_length=11, polyorder=3)
    return grid, values


def extract_cell_extrema(case: str, x: np.ndarray, y: np.ndarray) -> list[Extremum]:
    grid, values = smooth_profile(x, y)
    peak_indices: list[int] = []
    boundary_hits: list[bool] = []

    for lo, hi in PEAK_INTERVALS:
        in_window = np.flatnonzero((grid >= lo) & (grid <= hi))
        index = int(in_window[np.argmax(values[in_window])])
        peak_indices.append(index)
        boundary_hits.append(min(grid[index] - lo, hi - grid[index]) < 0.03)

    # A preceding expansion minimum is defined before each compression peak.
    trough_indices: list[int] = []
    left = 0.55
    for peak_index in peak_indices:
        right = grid[peak_index]
        in_window = np.flatnonzero((grid >= left) & (grid <= right))
        trough_indices.append(int(in_window[np.argmin(values[in_window])]))
        left = right

    # A following minimum is needed for a local prominence estimate.
    following_troughs: list[int] = []
    for i, peak_index in enumerate(peak_indices):
        right = grid[peak_indices[i + 1]] if i + 1 < len(peak_indices) else 6.50
        in_window = np.flatnonzero((grid >= grid[peak_index]) & (grid <= right))
        following_troughs.append(int(in_window[np.argmin(values[in_window])]))

    rows: list[Extremum] = []
    previous_resolved_peak: float | None = None
    for i, (peak_index, trough_index, following_index) in enumerate(
        zip(peak_indices, trough_indices, following_troughs)
    ):
        peak = float(values[peak_index])
        trough = float(values[trough_index])
        peak_to_trough = peak - trough
        prominence = peak - 0.5 * (trough + float(values[following_index]))
        resolved = bool(
            not boundary_hits[i] and peak_to_trough >= 0.08 and prominence >= 0.04
        )
        spacing = None
        if resolved and previous_resolved_peak is not None:
            spacing = float(grid[peak_index] - previous_resolved_peak)
        if resolved:
            previous_resolved_peak = float(grid[peak_index])
        rows.append(
            Extremum(
                case=case,
                cell=i + 1,
                peak_x=float(grid[peak_index]),
                peak_value=peak,
                preceding_trough_x=float(grid[trough_index]),
                preceding_trough_value=trough,
                peak_to_trough=peak_to_trough,
                prominence=prominence,
                spacing_from_previous=spacing,
                resolved=resolved,
            )
        )
    return rows


def experimental_radial_profile(exp5, station: float) -> tuple[np.ndarray, np.ndarray]:
    key = next(
        key
        for key in exp5.files
        if key.startswith("x") and abs(float(key[1:]) - station) < 1.0e-8
    )
    data = exp5[key]
    radius = np.abs(data[:, 0])
    density = data[:, 1] / float(exp5["rhoj"])
    order = np.argsort(radius)
    return radius[order], density[order]


def radial_profile_rms(
    radius: np.ndarray,
    profile: np.ndarray,
    exp_radius: np.ndarray,
    exp_density: np.ndarray,
) -> float:
    valid = np.isfinite(profile)
    interpolated = np.interp(exp_radius, radius[valid], profile[valid])
    return rms(interpolated, exp_density)


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def scheme_and_window_analysis() -> dict:
    teno_ax = load("cmp_teno16500_axstrip.npz")
    weno_a_ax = load("cmp_weno14750_axstrip.npz")
    weno_b_ax = load("cmp_weno16500_axstrip.npz")
    z = teno_ax["zz"]
    teno_axis = teno_ax["rho"][:, 0] / RHO_J_SIM
    weno_block_axis = block_mean(
        weno_a_ax["rho"][:, 0],
        float(weno_a_ax["time"]),
        weno_b_ax["rho"][:, 0],
        float(weno_b_ax["time"]),
    ) / RHO_J_SIM

    cumulative_axis: list[tuple[float, np.ndarray, str]] = [
        (
            float(weno_a_ax["time"]),
            weno_a_ax["rho"][:, 0] / RHO_J_SIM,
            "plt14750",
        ),
        (
            float(weno_b_ax["time"]),
            weno_b_ax["rho"][:, 0] / RHO_J_SIM,
            "plt16500",
        ),
    ]
    for plot in ["17500", "18500", "20500", "23500", "27500", "32080"]:
        data = load(f"win_plt{plot}.npz")
        cumulative_axis.append(
            (float(data["time"]), data["axis"] / RHO_J_SIM, f"plt{plot}")
        )

    adjacent_blocks: list[tuple[float, float, np.ndarray]] = []
    for (time_a, mean_a, _), (time_b, mean_b, _) in zip(
        cumulative_axis[:-1], cumulative_axis[1:]
    ):
        adjacent_blocks.append(
            (time_a, time_b, block_mean(mean_a, time_a, mean_b, time_b))
        )

    final_axis = cumulative_axis[-1][1]
    axial_region = (z >= 0.90) & (z <= 6.00)
    scheme_axis_rms = rms(teno_axis[axial_region], weno_block_axis[axial_region])
    temporal_axis_rms = rms(
        adjacent_blocks[0][2][axial_region], adjacent_blocks[1][2][axial_region]
    )

    convergence_rows: list[dict] = []
    convergence_x = (z >= 0.60) & (z <= 5.00)
    for time, profile, name in cumulative_axis:
        convergence_rows.append(
            {
                "endpoint": name,
                "time_ms": time * 1.0e3,
                "cumulative_window_ms": (time - WENO_STATS_START) * 1.0e3,
                "axis_rms_to_final_x0p6_5": rms(
                    profile[convergence_x], final_axis[convergence_x]
                ),
            }
        )
    write_csv(
        OUT / "weno_cumulative_convergence.csv",
        convergence_rows,
        list(convergence_rows[0]),
    )

    teno_rad = load("cmp_teno16500_radscan.npz")
    weno_a_rad = load("cmp_weno14750_radscan.npz")
    weno_b_rad = load("cmp_weno16500_radscan.npz")
    weno_full_rad = load("m142_2d_radscan_wide.npz")
    next_rad = load("win_plt17500.npz")
    exp5 = load("panda_fig5a_mj144_EPAPS.npz")
    radius = teno_rad["rr"]
    compare_radius = radius <= 1.20

    scheme_rows: list[dict] = []
    for station in RADIAL_STATIONS_8:
        key = f"s{station:.2f}_o+0.00"
        teno_profile = teno_rad[key] / RHO_J_SIM
        weno_block = block_mean(
            weno_a_rad[key],
            float(weno_a_rad["time"]),
            weno_b_rad[key],
            float(weno_b_rad["time"]),
        ) / RHO_J_SIM
        next_block = block_mean(
            weno_b_rad[key] / RHO_J_SIM,
            float(weno_b_rad["time"]),
            next_rad[f"s{station:.2f}"] / RHO_J_SIM,
            float(next_rad["time"]),
        )
        exp_radius, exp_density = experimental_radial_profile(exp5, station)
        scheme_rows.append(
            {
                "x_over_D": station,
                "teno_weno_direct_rms_r0_1p2": rms(
                    teno_profile[compare_radius], weno_block[compare_radius]
                ),
                "weno_neighbor_block_rms_r0_1p2": rms(
                    weno_block[compare_radius], next_block[compare_radius]
                ),
                "teno_rms_vs_exp_mj1p44": radial_profile_rms(
                    radius, teno_profile, exp_radius, exp_density
                ),
                "weno_block_rms_vs_exp_mj1p44": radial_profile_rms(
                    radius, weno_block, exp_radius, exp_density
                ),
                "weno_full_rms_vs_exp_mj1p44": radial_profile_rms(
                    radius,
                    weno_full_rad[key] / RHO_J_SIM,
                    exp_radius,
                    exp_density,
                ),
            }
        )
    write_csv(OUT / "scheme_radial_metrics.csv", scheme_rows, list(scheme_rows[0]))

    exp_axis = load("panda_m142den_axis_EPAPS.npz")
    scheme_extrema: list[Extremum] = []
    scheme_extrema += extract_cell_extrema(
        "experiment_m142_phase_cycle_mean", exp_axis["x"], exp_axis["rho"]
    )
    scheme_extrema += extract_cell_extrema("TENO5_short_mean", z, teno_axis)
    scheme_extrema += extract_cell_extrema(
        "WENOZ5_matched_block", z, weno_block_axis
    )
    scheme_extrema += extract_cell_extrema("WENOZ5_final_mean", z, final_axis)

    # Extrema of every adjacent WENO block quantify how much a short record can
    # move the apparent shock-cell peaks without changing the scheme.
    block_extrema: list[list[Extremum]] = []
    for index, (time_a, time_b, profile) in enumerate(adjacent_blocks):
        rows = extract_cell_extrema(f"WENO_block_{index}", z, profile)
        block_extrema.append(rows)
        scheme_extrema += rows

    write_csv(
        OUT / "scheme_shock_cell_metrics.csv",
        [asdict(row) for row in scheme_extrema],
        list(asdict(scheme_extrema[0])),
    )

    fig, axes = plt.subplots(2, 2, figsize=(14.0, 9.2))
    ax = axes[0, 0]
    ax.plot(z, final_axis, color="0.70", lw=1.8, label="WENO-Z5, 9.25 ms mean")
    ax.plot(
        z,
        weno_block_axis,
        color="0.20",
        lw=1.7,
        ls="--",
        label="WENO-Z5, 0.849 ms block",
    )
    ax.plot(z, teno_axis, color="#0066b3", lw=2.0, label="TENO5, 0.889 ms mean")
    ax.plot(
        exp_axis["x"],
        exp_axis["rho"],
        "o",
        ms=4.0,
        mfc="white",
        mec="black",
        label="Panda, Mj=1.42 phase-cycle mean",
    )
    ax.set(xlim=(0.5, 5.2), ylim=(0.3, 1.45), xlabel="x/D", ylabel=r"$\bar{\rho}/\rho_j$")
    ax.set_title("(a) Scheme comparison: short records are not converged")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8.3)

    ax = axes[0, 1]
    stations = np.array([row["x_over_D"] for row in scheme_rows])
    direct = np.array([row["teno_weno_direct_rms_r0_1p2"] for row in scheme_rows])
    temporal = np.array([row["weno_neighbor_block_rms_r0_1p2"] for row in scheme_rows])
    width = 0.055
    ax.bar(stations - width / 2, direct, width=width, color="#0066b3", label="TENO-WENO")
    ax.bar(
        stations + width / 2,
        temporal,
        width=width,
        color="#d55e00",
        label="WENO adjacent blocks",
    )
    ax.set(xlabel="x/D", ylabel="profile RMS difference, r/D <= 1.2")
    ax.set_title("(b) Scheme difference is below finite-block variability")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(fontsize=8.5)

    ax = axes[1, 0]
    windows = np.array([row["cumulative_window_ms"] for row in convergence_rows])
    convergence = np.array(
        [row["axis_rms_to_final_x0p6_5"] for row in convergence_rows]
    )
    ax.plot(windows[:-1], convergence[:-1], "o-", color="#009e73", lw=1.8)
    ax.axvline(
        (float(teno_ax["time"]) - TENO_STATS_START) * 1.0e3,
        color="#0066b3",
        ls="--",
        lw=1.2,
        label="TENO record length",
    )
    ax.set(
        xscale="log",
        xlabel="WENO cumulative averaging window (ms)",
        ylabel="axis RMS to 9.25 ms mean, 0.6 <= x/D <= 5",
    )
    ax.set_title("(c) WENO mean convergence")
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=8.5)

    ax = axes[1, 1]
    mid_times = np.array([(a + b) * 0.5e3 for a, b, _ in adjacent_blocks])
    lengths = np.array([(b - a) * 1.0e3 for a, b, _ in adjacent_blocks])
    peak2 = np.array([rows[1].peak_x for rows in block_extrema])
    amp2 = np.array([rows[1].peak_to_trough for rows in block_extrema])
    scatter = ax.scatter(
        mid_times,
        peak2,
        c=lengths,
        cmap="viridis",
        s=55,
        edgecolor="black",
        linewidth=0.5,
    )
    ax.axhline(
        extract_cell_extrema("exp", exp_axis["x"], exp_axis["rho"])[1].peak_x,
        color="black",
        ls=":",
        lw=1.2,
        label="experiment",
    )
    ax.plot(
        [0.5 * (TENO_STATS_START + float(teno_ax["time"])) * 1.0e3],
        [extract_cell_extrema("teno", z, teno_axis)[1].peak_x],
        marker="*",
        ms=13,
        color="#0066b3",
        mec="black",
        label="TENO short mean",
    )
    ax.set(xlabel="block midpoint time (ms)", ylabel="second compression maximum x/D")
    ax.set_title("(d) Apparent peak location varies between WENO blocks")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8.0)
    colorbar = fig.colorbar(scatter, ax=ax, pad=0.02)
    colorbar.set_label("block length (ms)")

    fig.suptitle("Mj=1.42 scheme and finite-window audit", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(OUT / "scheme_and_window_audit.png", dpi=180)
    plt.close(fig)

    return {
        "teno_window_ms": (float(teno_ax["time"]) - TENO_STATS_START) * 1.0e3,
        "weno_matched_block_ms": (
            float(weno_b_ax["time"]) - float(weno_a_ax["time"])
        )
        * 1.0e3,
        "scheme_axis_rms_x0p9_6": scheme_axis_rms,
        "neighbor_weno_block_axis_rms_x0p9_6": temporal_axis_rms,
        "median_scheme_radial_rms": float(np.median(direct)),
        "median_neighbor_block_radial_rms": float(np.median(temporal)),
        "weno_3p553ms_axis_rms_to_final": next(
            row["axis_rms_to_final_x0p6_5"]
            for row in convergence_rows
            if row["endpoint"] == "plt20500"
        ),
    }


def dimensionality_analysis() -> dict:
    exp_axis = load("panda_m142den_axis_EPAPS.npz")
    two_d = load("win_plt20500.npz")
    three_d = load("m142_3d_axrings.npz")
    x2 = two_d["zz"]
    y2 = two_d["axis"] / RHO_J_SIM
    x3 = three_d["xx"]
    y3 = three_d["ring00"] / RHO_J_SIM
    xe = exp_axis["x"]
    ye = exp_axis["rho"]

    extrema_rows: list[Extremum] = []
    extrema_by_case = {
        "experiment_m142_phase_cycle_mean": extract_cell_extrema(
            "experiment_m142_phase_cycle_mean", xe, ye
        ),
        "2D_WENO_3p553ms_mean": extract_cell_extrema(
            "2D_WENO_3p553ms_mean", x2, y2
        ),
        "3D_WENO_3p660ms_near_axis": extract_cell_extrema(
            "3D_WENO_3p660ms_near_axis", x3, y3
        ),
    }
    for rows in extrema_by_case.values():
        extrema_rows += rows
    write_csv(
        OUT / "dimensionality_shock_cell_metrics.csv",
        [asdict(row) for row in extrema_rows],
        list(asdict(extrema_rows[0])),
    )

    axis_ranges = [(0.60, 1.50), (1.50, 3.00), (3.00, 5.00), (0.60, 5.00)]
    axis_error_rows: list[dict] = []
    for name, x, y in [("2D", x2, y2), ("3D", x3, y3)]:
        for lo, hi in axis_ranges:
            selected = (xe >= lo) & (xe <= hi)
            axis_error_rows.append(
                {
                    "case": name,
                    "x_min_over_D": lo,
                    "x_max_over_D": hi,
                    "axis_rms_vs_same_mach_experiment": rms(
                        np.interp(xe[selected], x, y), ye[selected]
                    ),
                }
            )
    write_csv(
        OUT / "dimensionality_axis_rms.csv",
        axis_error_rows,
        list(axis_error_rows[0]),
    )

    all2 = load("all30_2d.npz")
    all3 = load("all30_3d.npz")
    exp5 = load("panda_fig5a_mj144_EPAPS.npz")
    radial_rows: list[dict] = []
    for station in RADIAL_STATIONS_30:
        exp_radius, exp_density = experimental_radial_profile(exp5, float(station))
        profile2 = all2[f"s{station:.2f}"] / RHO_J_SIM
        profile3 = all3[f"s{station:.2f}"] / RHO_J_SIM
        radial_rows.append(
            {
                "x_over_D": station,
                "two_d_rms_vs_exp_mj1p44": radial_profile_rms(
                    all2["rr"], profile2, exp_radius, exp_density
                ),
                "three_d_rms_vs_exp_mj1p44": radial_profile_rms(
                    all3["rr"], profile3, exp_radius, exp_density
                ),
            }
        )
    write_csv(
        OUT / "dimensionality_radial_metrics_30stations.csv",
        radial_rows,
        list(radial_rows[0]),
    )

    region_masks = {
        "first_cell_x_le_1p25": RADIAL_STATIONS_30 <= 1.25,
        "closure_transition_1p4_to_2p3": (RADIAL_STATIONS_30 >= 1.40)
        & (RADIAL_STATIONS_30 <= 2.30),
        "cells_2_to_4_x2p45_to_4": (RADIAL_STATIONS_30 >= 2.45)
        & (RADIAL_STATIONS_30 <= 4.00),
        "downstream_x_gt_4": RADIAL_STATIONS_30 > 4.00,
    }
    radial2 = np.array([row["two_d_rms_vs_exp_mj1p44"] for row in radial_rows])
    radial3 = np.array([row["three_d_rms_vs_exp_mj1p44"] for row in radial_rows])
    region_rows: list[dict] = []
    for region, selected in region_masks.items():
        region_rows.append(
            {
                "region": region,
                "station_count": int(np.count_nonzero(selected)),
                "two_d_mean_rms": float(np.mean(radial2[selected])),
                "three_d_mean_rms": float(np.mean(radial3[selected])),
                "two_d_median_rms": float(np.median(radial2[selected])),
                "three_d_median_rms": float(np.median(radial3[selected])),
            }
        )
    write_csv(OUT / "dimensionality_radial_regions.csv", region_rows, list(region_rows[0]))

    fig = plt.figure(figsize=(14.2, 9.3))
    grid = fig.add_gridspec(2, 2, height_ratios=(1.1, 1.0), hspace=0.30, wspace=0.25)
    ax = fig.add_subplot(grid[0, :])
    ax.plot(x2, y2, color="0.35", lw=1.9, ls="--", label="2D WENO-Z5, 3.553 ms mean")
    ax.plot(
        x3,
        y3,
        color="#0066b3",
        lw=2.0,
        label="3D WENO-Z5, 3.660 ms near-axis proxy",
    )
    ax.plot(xe, ye, "o", ms=4.2, mfc="white", mec="black", label="Panda, Mj=1.42")
    ax.set(xlim=(0.5, 6.2), ylim=(0.3, 1.48), xlabel="x/D", ylabel=r"$\bar{\rho}/\rho_j$")
    ax.set_title("(a) Axial mean: 3D preserves cells but develops an upstream phase drift")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8.8)

    ax = fig.add_subplot(grid[1, 0])
    experiment_extrema = extrema_by_case["experiment_m142_phase_cycle_mean"]
    cells = np.arange(1, 5)
    exp_positions = np.array([row.peak_x for row in experiment_extrema[:4]])
    exp_amplitudes = np.array([row.peak_to_trough for row in experiment_extrema[:4]])
    for label, color, marker in [
        ("2D_WENO_3p553ms_mean", "0.35", "s"),
        ("3D_WENO_3p660ms_near_axis", "#0066b3", "o"),
    ]:
        rows = extrema_by_case[label][:4]
        errors = np.array(
            [row.peak_x - exp_positions[i] if row.resolved else np.nan for i, row in enumerate(rows)]
        )
        ax.plot(cells, errors, marker=marker, color=color, lw=1.8, label=label.split("_")[0])
    ax.axhline(0.0, color="black", lw=1.0)
    ax.set(xticks=cells, xlabel="compression maximum number", ylabel="position error, x_sim/D - x_exp/D")
    ax.set_title("(b) After cell 1, compression maxima move progressively upstream in 3D")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8.5)

    ax = fig.add_subplot(grid[1, 1])
    for label, color, marker in [
        ("2D_WENO_3p553ms_mean", "0.35", "s"),
        ("3D_WENO_3p660ms_near_axis", "#0066b3", "o"),
    ]:
        rows = extrema_by_case[label][:4]
        ratios = np.array(
            [
                row.peak_to_trough / exp_amplitudes[i] if row.resolved else np.nan
                for i, row in enumerate(rows)
            ]
        )
        ax.plot(cells, ratios, marker=marker, color=color, lw=1.8, label=label.split("_")[0])
    ax.axhline(1.0, color="black", lw=1.0, ls=":")
    ax.set(xticks=cells, xlabel="shock cell", ylabel="peak-to-trough amplitude / experiment")
    ax.set_title("(c) 2D loses amplitude; 3D matches cells 2-3 but overshoots cell 4")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8.5)

    fig.suptitle("Mj=1.42 dimensionality audit (different native AMR resolutions)", fontsize=14)
    fig.savefig(OUT / "dimensionality_axis_audit.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    representative = [1.05, 1.55, 2.00, 2.75, 4.00, 5.00]
    fig, axes = plt.subplots(2, 3, figsize=(13.8, 7.6), sharex=True, sharey=True)
    for ax, station in zip(axes.flat, representative):
        exp_radius, exp_density = experimental_radial_profile(exp5, station)
        ax.plot(
            all2["rr"],
            all2[f"s{station:.2f}"] / RHO_J_SIM,
            color="0.35",
            lw=1.6,
            ls="--",
            label="2D",
        )
        ax.plot(
            all3["rr"],
            all3[f"s{station:.2f}"] / RHO_J_SIM,
            color="#0066b3",
            lw=1.8,
            label="3D four-arm mean",
        )
        ax.plot(exp_radius, exp_density, "o", ms=3.8, mfc="white", mec="black", label="Panda Mj=1.44")
        ax.set_title(f"x/D = {station:g}")
        ax.set(xlim=(0, 1.25), ylim=(0.2, 1.52))
        ax.grid(alpha=0.25)
    for ax in axes[-1, :]:
        ax.set_xlabel("r/D")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\bar{\rho}/\rho_j$")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        ncol=3,
        fontsize=9.0,
        bbox_to_anchor=(0.5, 0.995),
    )
    fig.suptitle(
        "Representative radial profiles; nominal-station RMS alone is phase-sensitive",
        fontsize=13,
        y=0.945,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(OUT / "dimensionality_radial_profiles.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(12.6, 4.8))
    ax.plot(RADIAL_STATIONS_30, radial2, "s--", color="0.35", lw=1.5, ms=4.5, label="2D")
    ax.plot(RADIAL_STATIONS_30, radial3, "o-", color="#0066b3", lw=1.7, ms=4.5, label="3D four-arm mean")
    ax.set(xlabel="nominal radial-profile station x/D", ylabel="profile RMS vs Panda Mj=1.44")
    ax.set_title("Raw radial-profile mismatch: useful locally, not a dimensionality ranking")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "dimensionality_radial_rms.png", dpi=180)
    plt.close(fig)

    return {
        "two_d_axis_rms_same_mach_x0p6_5": next(
            row["axis_rms_vs_same_mach_experiment"]
            for row in axis_error_rows
            if row["case"] == "2D" and row["x_min_over_D"] == 0.60 and row["x_max_over_D"] == 5.00
        ),
        "three_d_axis_rms_same_mach_x0p6_5": next(
            row["axis_rms_vs_same_mach_experiment"]
            for row in axis_error_rows
            if row["case"] == "3D" and row["x_min_over_D"] == 0.60 and row["x_max_over_D"] == 5.00
        ),
        "radial_regions": region_rows,
        "note": "3D data are a four-arm azimuthal sample and a composite AMR profile, not a full circular average at uniform D/256.",
    }


def axis_core_analysis() -> dict:
    fields = load("m142_2d_fields.npz")
    components = [str(value) for value in fields["comps"]]
    component_index = {name: i for i, name in enumerate(components)}
    data = fields["fields"].astype(np.float64)
    dx = float(fields["dxt"])
    x = (np.arange(data.shape[1]) + 0.5) * dx / D
    radius = (np.arange(data.shape[2]) + 0.5) * dx / D
    density = data[component_index["DensityMEAN"]]
    pressure = data[component_index["pressureMEAN"]]

    all2 = load("all30_2d.npz")
    all3 = load("all30_3d.npz")
    rows: list[dict] = []
    for station in RADIAL_STATIONS_30:
        j = int(np.argmin(np.abs(x - station)))
        rho_axis = float(density[j, 0])
        rho_r005 = float(np.interp(0.05, radius, density[j]))
        p_axis = float(pressure[j, 0])
        p_r005 = float(np.interp(0.05, radius, pressure[j]))
        temperature_axis_proxy = p_axis / (R_AIR * rho_axis)
        temperature_r005_proxy = p_r005 / (R_AIR * rho_r005)
        entropy_proxy = np.log(p_axis / p_r005) - GAMMA * np.log(rho_axis / rho_r005)

        profile2 = all2[f"s{station:.2f}"] / RHO_J_SIM
        profile3 = all3[f"s{station:.2f}"] / RHO_J_SIM
        contrast2 = float(np.interp(0.05, all2["rr"], profile2) - profile2[0])
        contrast3 = float(np.interp(0.05, all3["rr"], profile3) - profile3[0])
        rows.append(
            {
                "x_over_D": station,
                "two_d_density_axis_over_rhoj": rho_axis / RHO_J_SIM,
                "two_d_density_r0p05_over_rhoj": rho_r005 / RHO_J_SIM,
                "two_d_density_deficit_r0p05_minus_axis": (rho_r005 - rho_axis) / RHO_J_SIM,
                "three_d_density_deficit_r0p05_minus_axis": contrast3,
                "two_d_pressure_axis_over_pamb": p_axis / P_AMB,
                "two_d_pressure_r0p05_over_pamb": p_r005 / P_AMB,
                "two_d_temperature_axis_proxy_K": temperature_axis_proxy,
                "two_d_temperature_r0p05_proxy_K": temperature_r005_proxy,
                "two_d_temperature_proxy_excess_K": temperature_axis_proxy - temperature_r005_proxy,
                "two_d_mean_field_entropy_proxy": entropy_proxy,
                "two_d_density_contrast_from_D1024_profile": contrast2,
            }
        )
    write_csv(OUT / "axis_core_metrics.csv", rows, list(rows[0]))

    stations = np.array([row["x_over_D"] for row in rows])
    contrast2 = np.array([row["two_d_density_deficit_r0p05_minus_axis"] for row in rows])
    contrast3 = np.array([row["three_d_density_deficit_r0p05_minus_axis"] for row in rows])
    pressure_delta = np.array(
        [row["two_d_pressure_axis_over_pamb"] - row["two_d_pressure_r0p05_over_pamb"] for row in rows]
    )
    temperature_delta = np.array([row["two_d_temperature_proxy_excess_K"] for row in rows])
    entropy_proxy = np.array([row["two_d_mean_field_entropy_proxy"] for row in rows])

    fig, axes = plt.subplots(2, 2, figsize=(13.8, 8.8))
    plot_region = stations >= 1.05
    ax = axes[0, 0]
    ax.plot(
        stations[plot_region],
        contrast2[plot_region],
        "s--",
        color="0.35",
        lw=1.6,
        label="2D axisymmetric",
    )
    ax.plot(
        stations[plot_region],
        contrast3[plot_region],
        "o-",
        color="#0066b3",
        lw=1.7,
        label="3D four-arm sample average",
    )
    ax.axhline(0.0, color="black", lw=0.9)
    ax.set(xlabel="x/D", ylabel=r"$\bar{\rho}(r/D=0.05)-\bar{\rho}(axis)$, normalized")
    ax.set_title("(a) Persistent low-density core is much stronger in 2D")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8.5)

    ax = axes[0, 1]
    ax.plot(stations[plot_region], pressure_delta[plot_region], "o-", color="#d55e00", lw=1.6)
    ax.axhline(0.0, color="black", lw=0.9)
    ax.set(xlabel="x/D", ylabel=r"$\bar{p}_{axis}/p_a-\bar{p}_{r/D=0.05}/p_a$")
    ax.set_title("(b) The 2D pressure is nearly flat or slightly higher on axis")
    ax.grid(alpha=0.25)

    ax = axes[1, 0]
    ax.plot(stations[plot_region], temperature_delta[plot_region], "o-", color="#cc79a7", lw=1.7)
    ax.axhline(0.0, color="black", lw=0.9)
    ax.set(xlabel="x/D", ylabel="mean-field temperature proxy excess (K)")
    ax.set_title("(c) Ratio-of-means temperature proxy is elevated on axis")
    ax.grid(alpha=0.25)

    ax = axes[1, 1]
    ax.plot(stations[plot_region], entropy_proxy[plot_region], "o-", color="#009e73", lw=1.7)
    ax.axhline(0.0, color="black", lw=0.9)
    ax.set(
        xlabel="x/D",
        ylabel=r"$\ln(\bar{p}_{axis}/\bar{p}_{0.05})-\gamma\ln(\bar{\rho}_{axis}/\bar{\rho}_{0.05})$",
    )
    ax.set_title("(d) Positive mean-field entropy proxy; mechanism is not yet isolated")
    ax.grid(alpha=0.25)

    fig.suptitle("Near-axis density-deficit diagnostic", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(OUT / "axis_core_diagnostic.png", dpi=180)
    plt.close(fig)

    post_first_closure = stations >= 1.25
    downstream = stations >= 3.05
    return {
        "two_d_max_density_deficit_x_ge_1p25": float(np.max(contrast2[post_first_closure])),
        "two_d_median_density_deficit_x_ge_3p05": float(np.median(contrast2[downstream])),
        "three_d_median_density_deficit_x_ge_3p05": float(np.median(contrast3[downstream])),
        "two_d_max_temperature_proxy_excess_K_x_ge_1p25": float(
            np.max(temperature_delta[post_first_closure])
        ),
        "two_d_max_mean_field_entropy_proxy_x_ge_1p25": float(
            np.max(entropy_proxy[post_first_closure])
        ),
        "caveat": "pbar/(R*rhobar) and ln(pbar)-gamma*ln(rhobar) are mean-field proxies, not directly averaged thermodynamic quantities.",
    }


def spatial_sampling_analysis() -> dict:
    rings = load("m142_3d_axrings.npz")
    fig4 = load("panda_fig4_centerline_EPAPS.npz")
    experiment = fig4["M1.43"]
    exp_x = experiment[:, 0]
    exp_y = experiment[:, 1] / float(fig4["rhoj_M1.43"])
    selected = exp_x <= 5.0

    ring_values = []
    ring_counts = []
    for index in range(16):
        ring_values.append(rings[f"ring{index:02d}"] / RHO_J_SIM)
        ring_counts.append(int(rings[f"ring{index:02d}_n"]))
    ring_values = np.stack(ring_values)
    ring_counts = np.array(ring_counts)

    radius_rows: list[dict] = []
    for index in range(16):
        weights = ring_counts[: index + 1]
        disc = np.average(ring_values[: index + 1], axis=0, weights=weights)
        radius_rows.append(
            {
                "nominal_disc_radius_over_D": 0.005 + 0.01 * index,
                "rms_vs_fig4_mj1p43_x_le_5": rms(
                    np.interp(exp_x[selected], rings["xx"], disc), exp_y[selected]
                ),
            }
        )
    write_csv(OUT / "three_d_disc_radius_sensitivity.csv", radius_rows, list(radius_rows[0]))
    return {
        "axis_proxy_rms": radius_rows[0]["rms_vs_fig4_mj1p43_x_le_5"],
        "r0p15_disc_rms": radius_rows[-1]["rms_vs_fig4_mj1p43_x_le_5"],
        "interpretation": "This is a spatial sampling sensitivity test only. It does not establish temporal convergence.",
    }


def main() -> None:
    summary = {
        "provenance": {
            "simulation_density_normalization_kg_m3": RHO_J_SIM,
            "weno_statistics_start_ms": WENO_STATS_START * 1.0e3,
            "teno_statistics_start_ms": TENO_STATS_START * 1.0e3,
            "three_d_statistics_window_ms": (THREED_STATS_END - THREED_STATS_START) * 1.0e3,
            "axis_experiment": "Panda m142den, Mj=1.42, 37 phase bins averaged over phase",
            "radial_experiment": "Panda figure 5a, Mj=1.44",
        },
        "scheme_and_window": scheme_and_window_analysis(),
        "dimensionality": dimensionality_analysis(),
        "axis_core": axis_core_analysis(),
        "spatial_sampling": spatial_sampling_analysis(),
    }
    with (OUT / "summary_metrics.json").open("w", encoding="ascii") as stream:
        json.dump(summary, stream, indent=2)
        stream.write("\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
