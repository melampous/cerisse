#!/usr/bin/env python3
"""Compare the first eight R-Z rings in the Run165 axis diagnostics.

This is a read-only post-processing utility.  It deliberately separates two
signatures:

* a curved wave has a smoothly varying axial gradient locus as ring number
  increases;
* a grid-attached filament is a non-monotone radial extremum that persists at
  one ring number over a long axial interval.

The script writes reproducible CSV tables and compact figures under its own
analysis directory.  It never modifies a plotfile or production source.
"""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[3]
ANALYSIS = REPO / "analysis/skew_axis_diagnosis"
RUN165 = REPO / "exm/underexpanded_jet/2d/run165_analytic_exit_matrix"
OUTPUT = Path(__file__).resolve().parent

CASES = {
    "skew_baseline_500us": ANALYSIS / "baseline/plt01250",
    "skew_baseline_1500us": ANALYSIS / "baseline/plt03750",
    "skew_K1_no_axis_sync_500us": ANALYSIS / "fixed_collar_500us/plt01250",
    "skew_K1_no_axis_sync_1500us": ANALYSIS / "fixed_collar_fields/plt03750",
    "skew_K1_axis_pair_500us": ANALYSIS / "fixed_axis_pair_500us/plt01250",
    "skew_K1_parity_axis_500us": (
        ANALYSIS / "fixed_axis_parity_500us/plt01250"
    ),
    "hllc_muscl_500us": (
        RUN165 / "results_hllc_muscl_o2_rz_cfl02_500us/plt01178"
    ),
}

FIELDS = (
    "pressure",
    "density",
    "temperature",
    "radial_velocity",
    "mach",
)

LABELS = {
    "pressure": r"$p/p_\infty$",
    "density": r"$\rho/\rho_\infty$",
    "temperature": r"$T/T_\infty$",
    "radial_velocity": r"$u_r$ [m/s]",
    "mach": "Mach",
}

P_INF = 534.58067022560419
RHO_INF = 0.0287655648
T_INF = 64.72986748
NORMALIZER = {
    "pressure": P_INF,
    "density": RHO_INF,
    "temperature": T_INF,
    "radial_velocity": 1.0,
    "mach": 1.0,
}

NEAR_Z_LO = -0.080
NEAR_Z_HI = -0.005
FILAMENT_THRESHOLD = 0.02


def load_axis_module():
    path = ANALYSIS / "analyze_axis_rings.py"
    spec = importlib.util.spec_from_file_location("axis_ring_reader", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.patch_yt_cylindrical_readonly_edge()
    return module


def write_rows(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def longest_true_run(mask: np.ndarray, spacing: float) -> float:
    padded = np.concatenate(([False], mask.astype(bool), [False]))
    changes = np.diff(padded.astype(np.int8))
    starts = np.flatnonzero(changes == 1)
    ends = np.flatnonzero(changes == -1)
    if starts.size == 0:
        return 0.0
    return float(np.max(ends - starts) * spacing)


def radial_extremum_residual(values: np.ndarray, floor: float) -> np.ndarray:
    """Signed non-monotone residual for rings 1..6.

    A smooth monotone radial wave crossing can have a large first derivative,
    but it should not create a long-lived alternating local extremum at one
    fixed ring.  This second-difference-like diagnostic targets the latter.
    """

    left = values[:-2]
    centre = values[1:-1]
    right = values[2:]
    scale = np.maximum.reduce(
        (
            np.abs(left),
            np.abs(centre),
            np.abs(right),
            np.full_like(centre, floor),
        )
    )
    return (centre - 0.5 * (left + right)) / scale


def first8(module, plotfile: Path) -> dict:
    dataset, level, fields = module.load_fields(plotfile)
    nr, nz = fields["pressure"].shape
    if nr < 8:
        raise RuntimeError(f"{plotfile}: only {nr} radial rings")
    r_lo = float(dataset.domain_left_edge[0])
    r_hi = float(dataset.domain_right_edge[0])
    z_lo = float(dataset.domain_left_edge[1])
    z_hi = float(dataset.domain_right_edge[1])
    dr = (r_hi - r_lo) / nr
    dz = (z_hi - z_lo) / nz
    r = r_lo + (np.arange(nr) + 0.5) * dr
    z = z_lo + (np.arange(nz) + 0.5) * dz
    return {
        "dataset": dataset,
        "level": level,
        "r": r[:8],
        "z": z,
        "dr": dr,
        "dz": dz,
        "fields": {name: fields[name][:8].copy() for name in FIELDS},
        "sound_speed": fields["sound_speed"][:8].copy(),
        "axial_velocity": fields["axial_velocity"][:8].copy(),
    }


def export_raw(data: dict[str, dict]) -> None:
    rows: list[dict] = []
    for case, item in data.items():
        time_s = float(item["dataset"].current_time)
        for ring, radius in enumerate(item["r"]):
            for j, axial in enumerate(item["z"]):
                row = {
                    "case": case,
                    "time_s": f"{time_s:.17g}",
                    "level": item["level"],
                    "ring": ring,
                    "r_mm": f"{radius * 1.0e3:.12g}",
                    "z_mm": f"{axial * 1.0e3:.12g}",
                }
                for field in FIELDS:
                    row[field] = f"{item['fields'][field][ring, j]:.17g}"
                row["axial_velocity"] = (
                    f"{item['axial_velocity'][ring, j]:.17g}"
                )
                rows.append(row)
    write_rows(OUTPUT / "first8_values.csv", rows)


def export_filament_metrics(data: dict[str, dict]) -> None:
    rows: list[dict] = []
    floors = {
        "pressure": 1.0,
        "density": 1.0e-12,
        "temperature": 1.0e-12,
        "radial_velocity": 1.0,
        "mach": 1.0e-12,
    }
    for case, item in data.items():
        z = item["z"]
        active = (z >= NEAR_Z_LO) & (z <= NEAR_Z_HI)
        active_indices = np.flatnonzero(active)
        for field in FIELDS:
            residual = radial_extremum_residual(
                item["fields"][field], floors[field]
            )
            for local_ring in range(residual.shape[0]):
                ring = local_ring + 1
                magnitude = np.abs(residual[local_ring, active])
                peak_local = int(np.argmax(magnitude))
                peak_j = int(active_indices[peak_local])
                above = magnitude > FILAMENT_THRESHOLD
                rows.append(
                    {
                        "case": case,
                        "time_s": f"{float(item['dataset'].current_time):.17g}",
                        "field": field,
                        "ring": ring,
                        "r_mm": f"{item['r'][ring] * 1.0e3:.12g}",
                        "peak_abs_radial_extremum": f"{magnitude[peak_local]:.17g}",
                        "peak_z_mm": f"{z[peak_j] * 1.0e3:.12g}",
                        "p95_abs_radial_extremum": (
                            f"{np.quantile(magnitude, 0.95):.17g}"
                        ),
                        "fraction_over_2pct": f"{np.mean(above):.17g}",
                        "longest_over_2pct_mm": (
                            f"{longest_true_run(above, item['dz']) * 1.0e3:.17g}"
                        ),
                    }
                )
    write_rows(OUTPUT / "radial_filament_metrics.csv", rows)


def export_axis_metrics(module, data: dict[str, dict]) -> None:
    rows: list[dict] = []
    scalar_fields = (
        "pressure",
        "density",
        "temperature",
        "axial_velocity",
    )
    floors = {
        "pressure": 1.0,
        "density": 1.0e-12,
        "temperature": 1.0e-12,
        "axial_velocity": 1.0,
    }
    for case, item in data.items():
        z = item["z"]
        active = (z >= NEAR_Z_LO) & (z <= NEAR_Z_HI)
        active_indices = np.flatnonzero(active)
        for field in scalar_fields:
            values = (
                item["axial_velocity"]
                if field == "axial_velocity"
                else item["fields"][field]
            )
            _, residual, relative = module.scalar_metrics(values, floors[field])
            local = int(np.argmax(relative[active]))
            j = int(active_indices[local])
            rows.append(
                {
                    "case": case,
                    "time_s": f"{float(item['dataset'].current_time):.17g}",
                    "field": field,
                    "axis_defect_max": f"{relative[j]:.17g}",
                    "z_mm": f"{z[j] * 1.0e3:.12g}",
                    "signed_residual": f"{residual[j]:.17g}",
                }
            )
        ur_prediction = module.odd_axis_prediction(
            item["fields"]["radial_velocity"]
        )
        ur_residual = item["fields"]["radial_velocity"][0] - ur_prediction
        relative = np.abs(ur_residual) / np.maximum(item["sound_speed"][0], 1.0)
        local = int(np.argmax(relative[active]))
        j = int(active_indices[local])
        rows.append(
            {
                "case": case,
                "time_s": f"{float(item['dataset'].current_time):.17g}",
                "field": "radial_velocity_over_sound",
                "axis_defect_max": f"{relative[j]:.17g}",
                "z_mm": f"{z[j] * 1.0e3:.12g}",
                "signed_residual": f"{ur_residual[j]:.17g}",
            }
        )
    write_rows(OUTPUT / "axis_parity_metrics_near_exit.csv", rows)


def export_gradient_locus(data: dict[str, dict]) -> None:
    rows: list[dict] = []
    for case, item in data.items():
        z = item["z"]
        active = (z >= -0.045) & (z <= NEAR_Z_HI)
        active_indices = np.flatnonzero(active)
        for field in ("pressure", "mach"):
            for ring in range(8):
                values = item["fields"][field][ring]
                if field == "pressure":
                    values = np.log(np.maximum(values, 1.0e-300))
                gradient = np.abs(np.gradient(values, item["dz"]))
                local = int(np.argmax(gradient[active]))
                j = int(active_indices[local])
                rows.append(
                    {
                        "case": case,
                        "time_s": f"{float(item['dataset'].current_time):.17g}",
                        "field": field,
                        "ring": ring,
                        "r_mm": f"{item['r'][ring] * 1.0e3:.12g}",
                        "dominant_gradient_z_mm": f"{z[j] * 1.0e3:.12g}",
                        "gradient_magnitude_per_m": f"{gradient[j]:.17g}",
                    }
                )
    write_rows(OUTPUT / "near_exit_gradient_locus.csv", rows)


def export_time_stationarity(data: dict[str, dict]) -> None:
    pairs = (
        ("skew_baseline", "skew_baseline_500us", "skew_baseline_1500us"),
        (
            "skew_K1_no_axis_sync",
            "skew_K1_no_axis_sync_500us",
            "skew_K1_no_axis_sync_1500us",
        ),
    )
    rows: list[dict] = []
    for label, early_name, late_name in pairs:
        early = data[early_name]
        late = data[late_name]
        active = (early["z"] >= NEAR_Z_LO) & (early["z"] <= NEAR_Z_HI)
        for field in FIELDS:
            a = early["fields"][field][:, active]
            b = late["fields"][field][:, active]
            floor = 1.0 if field in ("pressure", "radial_velocity") else 1.0e-12
            scale = np.maximum.reduce(
                (np.abs(a), np.abs(b), np.full_like(a, floor))
            )
            relative = np.abs(b - a) / scale
            rows.append(
                {
                    "pair": label,
                    "field": field,
                    "max_abs_change": f"{np.max(np.abs(b-a)):.17g}",
                    "max_relative_change": f"{np.max(relative):.17g}",
                    "rms_relative_change": (
                        f"{np.sqrt(np.mean(relative * relative)):.17g}"
                    ),
                }
            )
    write_rows(OUTPUT / "stationarity_500_to_1500us.csv", rows)


def plot_profiles(data: dict[str, dict], names: list[str], output: Path) -> None:
    figure, axes = plt.subplots(
        len(FIELDS), len(names), figsize=(4.2 * len(names), 12.5),
        sharex=True, squeeze=False,
    )
    colors = plt.cm.viridis(np.linspace(0.03, 0.97, 8))
    for column, name in enumerate(names):
        item = data[name]
        z_mm = item["z"] * 1.0e3
        active = (z_mm >= -80.0) & (z_mm <= 0.5)
        for row, field in enumerate(FIELDS):
            axis = axes[row, column]
            values = item["fields"][field] / NORMALIZER[field]
            for ring in range(8):
                axis.plot(
                    z_mm[active], values[ring, active], color=colors[ring],
                    lw=1.15, label=f"r{ring} ({item['r'][ring]*1e3:.2f} mm)",
                )
            axis.grid(True, alpha=0.2)
            if column == 0:
                axis.set_ylabel(LABELS[field])
            if row == 0:
                axis.set_title(name.replace("_", "\n"), fontsize=10)
            if row == len(FIELDS) - 1:
                axis.set_xlabel("axial z [mm]")
    handles, labels = axes[0, -1].get_legend_handles_labels()
    figure.legend(
        handles, labels, loc="center right", bbox_to_anchor=(0.995, 0.5),
        fontsize=8,
    )
    figure.suptitle("Run165: first eight finest-level radial rings", y=0.995)
    figure.tight_layout(rect=(0, 0, 0.91, 0.985))
    figure.savefig(output, dpi=180)
    plt.close(figure)


def plot_filament_heatmaps(data: dict[str, dict]) -> None:
    names = [
        "skew_baseline_500us",
        "skew_K1_no_axis_sync_500us",
        "skew_K1_axis_pair_500us",
        "hllc_muscl_500us",
    ]
    fields = ("pressure", "density", "temperature", "mach")
    figure, axes = plt.subplots(
        len(fields), len(names), figsize=(15.5, 10.5), sharex=True,
        sharey=True, squeeze=False,
    )
    for column, name in enumerate(names):
        item = data[name]
        z_mm = item["z"] * 1.0e3
        active = (z_mm >= -80.0) & (z_mm <= -5.0)
        for row, field in enumerate(fields):
            residual = np.abs(radial_extremum_residual(
                item["fields"][field], 1.0 if field == "pressure" else 1.0e-12
            ))
            image = axes[row, column].imshow(
                residual[:, active], origin="lower", aspect="auto",
                extent=(z_mm[active][0], z_mm[active][-1], 0.5, 6.5),
                vmin=0.0, vmax=0.20, cmap="magma",
            )
            axes[row, column].set_yticks(range(1, 7))
            axes[row, column].grid(False)
            if column == 0:
                axes[row, column].set_ylabel(f"{field}\ncentre ring k")
            if row == 0:
                axes[row, column].set_title(name.replace("_", "\n"), fontsize=10)
            if row == len(fields) - 1:
                axes[row, column].set_xlabel("axial z [mm]")
    colorbar = figure.colorbar(image, ax=axes, fraction=0.02, pad=0.015)
    colorbar.set_label(r"$|q_k-(q_{k-1}+q_{k+1})/2|/$ local scale")
    figure.suptitle(
        "Grid-attached filament diagnostic: radial local-extremum residual",
        y=0.995,
    )
    figure.savefig(OUTPUT / "radial_filament_heatmaps_500us.png", dpi=180)
    plt.close(figure)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    missing = [str(path) for path in CASES.values() if not path.is_dir()]
    if missing:
        raise FileNotFoundError("missing plotfiles:\n" + "\n".join(missing))
    module = load_axis_module()
    data = {name: first8(module, path) for name, path in CASES.items()}
    export_raw(data)
    export_filament_metrics(data)
    export_axis_metrics(module, data)
    export_gradient_locus(data)
    export_time_stationarity(data)
    plot_profiles(
        data,
        [
            "skew_baseline_500us",
            "skew_K1_no_axis_sync_500us",
            "skew_K1_axis_pair_500us",
            "hllc_muscl_500us",
        ],
        OUTPUT / "first8_profiles_500us.png",
    )
    plot_profiles(
        data,
        ["skew_baseline_1500us", "skew_K1_no_axis_sync_1500us"],
        OUTPUT / "first8_profiles_1500us.png",
    )
    plot_filament_heatmaps(data)
    print(OUTPUT)


if __name__ == "__main__":
    main()
