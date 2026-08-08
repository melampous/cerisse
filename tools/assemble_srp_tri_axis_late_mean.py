#!/usr/bin/env python3
"""Assemble a non-uniform-time late mean from strict 4-axis snapshots."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatterMathtext
import numpy as np


GAMMA = 1.4
AXES = (
    ("body_axis", 0.275, 0.275, "body axis"),
    ("nozzle_0", 0.307, 0.275, "nozzle 0 axis"),
    ("nozzle_1", 0.259, 0.302, "nozzle 1 axis"),
    ("nozzle_2", 0.259, 0.248, "nozzle 2 axis"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot-root", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--steps", nargs="+", type=int,
                        default=[6000, 6500, 7000, 7500, 8000, 8500, 9000])
    parser.add_argument("--p-inf", type=float, default=574.56)
    parser.add_argument("--body-diameter", type=float, default=0.127)
    parser.add_argument("--body-centre-x", type=float, default=0.300)
    parser.add_argument("--nozzle-plane-x", type=float, default=0.299)
    parser.add_argument("--dpi", type=int, default=600)
    return parser.parse_args()


def trapezoidal_weights(times: np.ndarray) -> np.ndarray:
    times = np.asarray(times, dtype=np.float64)
    if times.ndim != 1 or times.size < 2 or np.any(np.diff(times) <= 0.0):
        raise ValueError("Snapshot times must be one-dimensional and strictly increasing")
    weights = np.empty_like(times)
    weights[0] = 0.5 * (times[1] - times[0])
    weights[-1] = 0.5 * (times[-1] - times[-2])
    weights[1:-1] = 0.5 * (times[2:] - times[:-2])
    weights /= times[-1] - times[0]
    if not np.isclose(np.sum(weights), 1.0, rtol=0.0, atol=5.0e-15):
        raise RuntimeError("Trapezoidal weights do not sum to one")
    return weights


def segment_ids(valid: np.ndarray) -> tuple[np.ndarray, list[tuple[int, int]]]:
    labels = np.full(valid.shape, -1, dtype=np.int32)
    padded = np.r_[False, valid, False]
    starts = np.flatnonzero(padded[1:] & ~padded[:-1])
    stops = np.flatnonzero(~padded[1:] & padded[:-1])
    segments = []
    for identifier, (start, stop) in enumerate(zip(starts, stops)):
        labels[start:stop] = identifier
        segments.append((int(start), int(stop)))
    return labels, segments


def segmented_abs_gradient(values, coords, segments):
    result = np.full(values.shape, np.nan, dtype=np.float64)
    for start, stop in segments:
        if stop - start >= 2:
            result[start:stop] = np.abs(
                np.gradient(values[start:stop], coords[start:stop], edge_order=1)
            )
    return result


def weighted_mean(stack: np.ndarray, weights: np.ndarray, valid: np.ndarray) -> np.ndarray:
    result = np.tensordot(weights, stack, axes=(0, 0))
    result = np.asarray(result, dtype=np.float64)
    result[~valid] = np.nan
    return result


def write_csv(path: Path, profile: dict[str, np.ndarray]) -> None:
    columns = list(profile)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        for row_index in range(len(profile["x_m"])):
            row = []
            for column in columns:
                value = profile[column][row_index]
                if np.issubdtype(profile[column].dtype, np.integer):
                    row.append(str(int(value)))
                else:
                    row.append("{:.12g}".format(float(value)))
            writer.writerow(row)


def finite_range(values):
    finite = values[np.isfinite(values)]
    return [float(np.min(finite)), float(np.max(finite))] if finite.size else [None, None]


def configure_matplotlib() -> None:
    matplotlib.rcParams.update({
        "font.family": "serif",
        "font.serif": ["STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": 8.5,
        "axes.labelsize": 9.0,
        "axes.titlesize": 9.0,
        "axes.linewidth": 0.65,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "legend.frameon": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def render_primary(profiles, path_stem, times, nozzle_x, p_inf, db, dpi):
    configure_matplotlib()
    colors = {
        "body_axis": "#111111", "nozzle_0": "#d62728",
        "nozzle_1": "#2474b7", "nozzle_2": "#269c3c",
    }
    labels = {
        name: "{} ($y={:.3f}$ m, $z={:.3f}$ m)".format(label, y, z)
        for name, y, z, label in AXES
    }
    fig, axes = plt.subplots(
        3, 1, figsize=(7.08, 7.5), sharex=True, constrained_layout=True,
        gridspec_kw={"height_ratios": [1.0, 1.08, 1.08]},
    )
    for name, profile in profiles.items():
        style = "-" if name == "body_axis" else "--"
        width = 1.35 if name == "body_axis" else 1.05
        axes[0].plot(profile["x_m"], profile["favre_Mach"], style,
                     color=colors[name], lw=width, label=labels[name])
        axes[1].plot(profile["x_m"], profile["mean_p_over_pinf"], style,
                     color=colors[name], lw=width)
        axes[2].plot(profile["x_m"], profile["abs_dmeanrho_dx_kg_m4"], style,
                     color=colors[name], lw=width)
    for index, axis in enumerate(axes):
        axis.axvline(nozzle_x, color="#666666", ls=":", lw=0.9,
                    label=r"nozzle exit plane, $x=0.299$ m" if index == 0 else None)
        axis.grid(True, which="both", color="#b8b8b8", alpha=0.24, lw=0.45)
    axes[0].set_ylabel(r"Favre Mach number, $\widetilde{M}$")
    axes[0].set_ylim(bottom=0.0)
    axes[0].legend(loc="upper right", fontsize=6.9, ncol=2, handlelength=2.4)
    axes[1].set_ylabel(r"Mean pressure, $\overline{p}/p_\infty$")
    axes[1].set_yscale("log")
    pressure_values = np.concatenate([
        profile["mean_p_over_pinf"][np.isfinite(profile["mean_p_over_pinf"])]
        for profile in profiles.values()
    ])
    if pressure_values.size:
        lower = 10.0 ** np.floor(np.log10(max(np.min(pressure_values), 1.0e-12)))
        upper = 10.0 ** np.ceil(np.log10(np.max(pressure_values)))
        axes[1].set_ylim(lower, upper)
    axes[1].yaxis.set_major_formatter(LogFormatterMathtext())
    right_axis = axes[1].secondary_yaxis(
        "right", functions=(lambda ratio: ratio * p_inf, lambda p: p / p_inf)
    )
    right_axis.set_ylabel(r"Mean pressure, $\overline{p}$ [Pa]")
    axes[2].set_ylabel(r"$|\partial \overline{\rho}/\partial x|$ [kg m$^{-4}$]")
    axes[2].set_yscale("log")
    axes[2].set_ylim(bottom=1.0e-3)
    axes[2].set_xlabel(r"Axial coordinate, $x$ [m]")
    top_axis = axes[0].secondary_xaxis(
        "top", functions=(lambda x: x / db, lambda x_over_db: x_over_db * db)
    )
    top_axis.set_xlabel(r"$x/D_b$")
    fig.suptitle(
        "SRP tri-nozzle late mean — trapezoidal $t={:.3f}$–${:.3f}$ ms".format(
            1.0e3 * times[0], 1.0e3 * times[-1]
        ),
        fontsize=10.0,
    )
    fig.savefig(path_stem.with_suffix(".png"), dpi=dpi)
    fig.savefig(path_stem.with_suffix(".pdf"))
    plt.close(fig)


def render_definition_comparison(profiles, output, times, dpi):
    configure_matplotlib()
    fig, axes = plt.subplots(2, 2, figsize=(7.08, 5.2), sharex=True,
                             constrained_layout=True)
    for axis, (name, y, z, label) in zip(axes.flat, AXES):
        profile = profiles[name]
        axis.plot(profile["x_m"], profile["favre_Mach"], color="#1f77b4",
                  lw=1.2, label=r"Favre-derived $\widetilde{M}$")
        axis.plot(profile["x_m"], profile["mean_instantaneous_Mach"],
                  color="#d62728", ls="--", lw=1.0,
                  label=r"$\overline{M(t)}$")
        axis.set_title("{}: $y={:.3f}$, $z={:.3f}$ m".format(label, y, z))
        axis.set_ylabel("Mach number")
        axis.grid(True, alpha=0.23, lw=0.45)
    for axis in axes[-1, :]:
        axis.set_xlabel("x [m]")
    axes[0, 0].legend(fontsize=7.0)
    fig.suptitle(
        "Late-mean definition check, $t={:.3f}$–${:.3f}$ ms".format(
            1.0e3 * times[0], 1.0e3 * times[-1]
        ), fontsize=9.5,
    )
    fig.savefig(output, dpi=dpi)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    snapshot_root = Path(args.snapshot_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    steps = sorted(set(args.steps))
    if len(steps) < 2:
        raise ValueError("At least two late snapshots are required")

    records = []
    for step in steps:
        path = snapshot_root / "plt{:05d}".format(step) / "axis_profiles_combined.npz"
        if not path.is_file():
            raise RuntimeError("Missing snapshot profile archive: {}".format(path))
        records.append(np.load(path, allow_pickle=False))
    try:
        times = np.asarray([float(record["simulation_time_s"]) for record in records])
        weights = trapezoidal_weights(times)
        profiles = {}
        metadata = {}
        combined = {
            "steps": np.asarray(steps, dtype=np.int64),
            "times_s": times,
            "trapezoidal_weights": weights,
            "p_inf_Pa": np.asarray(args.p_inf),
            "body_diameter_m": np.asarray(args.body_diameter),
            "body_centre_x_m": np.asarray(args.body_centre_x),
            "nozzle_exit_plane_x_m": np.asarray(args.nozzle_plane_x),
        }
        for name, y, z, label in AXES:
            def stack(field):
                return np.stack([record["{}__{}".format(name, field)]
                                 for record in records]).astype(np.float64)

            x = stack("x_m")[0]
            for record_index in range(1, len(records)):
                if not np.allclose(stack("x_m")[record_index], x,
                                   rtol=0.0, atol=1.0e-14):
                    raise RuntimeError("x-grid mismatch for {}".format(name))
            valid_stack = stack("fluid_valid").astype(bool)
            common_valid = np.all(valid_stack, axis=0)
            labels, segments = segment_ids(common_valid)

            rho_stack = stack("rho_kg_m3")
            p_stack = stack("p_Pa")
            temperature_stack = stack("T_K")
            u_stack = stack("u_m_s")
            v_stack = stack("v_m_s")
            w_stack = stack("w_m_s")
            mach_stack = stack("Mach")
            rho_mean = weighted_mean(rho_stack, weights, common_valid)
            pressure_mean = weighted_mean(p_stack, weights, common_valid)
            temperature_mean = weighted_mean(temperature_stack, weights, common_valid)
            u_mean = weighted_mean(u_stack, weights, common_valid)
            v_mean = weighted_mean(v_stack, weights, common_valid)
            w_mean = weighted_mean(w_stack, weights, common_valid)
            mach_mean = weighted_mean(mach_stack, weights, common_valid)
            mean_rho_u = weighted_mean(rho_stack * u_stack, weights, common_valid)
            mean_rho_v = weighted_mean(rho_stack * v_stack, weights, common_valid)
            mean_rho_w = weighted_mean(rho_stack * w_stack, weights, common_valid)
            mean_rho_T = weighted_mean(
                rho_stack * temperature_stack, weights, common_valid
            )
            tiny = np.finfo(np.float64).tiny
            favre_u = mean_rho_u / np.maximum(rho_mean, tiny)
            favre_v = mean_rho_v / np.maximum(rho_mean, tiny)
            favre_w = mean_rho_w / np.maximum(rho_mean, tiny)
            favre_temperature = mean_rho_T / np.maximum(rho_mean, tiny)
            mean_sound_speed = np.sqrt(
                GAMMA * pressure_mean / np.maximum(rho_mean, tiny)
            )
            favre_mach = np.sqrt(
                favre_u * favre_u + favre_v * favre_v + favre_w * favre_w
            ) / np.maximum(mean_sound_speed, tiny)
            for values in (favre_u, favre_v, favre_w, favre_temperature, favre_mach):
                values[~common_valid] = np.nan
            grad_mean_rho = segmented_abs_gradient(rho_mean, x, segments)
            mean_instantaneous_gradient = weighted_mean(
                stack("abs_drhodx_kg_m4"), weights, common_valid
            )
            sld_max = np.max(stack("sld"), axis=0)
            ghs_max = np.max(stack("ghs"), axis=0)
            source_level_min = np.min(stack("source_level"), axis=0).astype(np.uint8)

            profile = {
                "x_m": x,
                "x_over_Db": x / args.body_diameter,
                "X_body": (x - args.body_centre_x) / args.body_diameter,
                "favre_Mach": favre_mach,
                "mean_instantaneous_Mach": mach_mean,
                "mean_p_Pa": pressure_mean,
                "mean_p_over_pinf": pressure_mean / args.p_inf,
                "mean_rho_kg_m3": rho_mean,
                "mean_T_K": temperature_mean,
                "favre_T_K": favre_temperature,
                "mean_u_m_s": u_mean,
                "mean_v_m_s": v_mean,
                "mean_w_m_s": w_mean,
                "favre_u_m_s": favre_u,
                "favre_v_m_s": favre_v,
                "favre_w_m_s": favre_w,
                "mean_axis_velocity_m_s": u_mean.copy(),
                "favre_axis_velocity_m_s": favre_u.copy(),
                "abs_dmeanrho_dx_kg_m4": grad_mean_rho,
                "mean_abs_drhodx_kg_m4": mean_instantaneous_gradient,
                "sld_max": sld_max,
                "ghs_max": ghs_max,
                "fluid_valid_all_times": common_valid.astype(np.uint8),
                "fluid_segment_id": labels,
                "source_level_min": source_level_min,
            }
            profiles[name] = profile
            csv_name = "{}_late_mean_profile.csv".format(name)
            write_csv(output_dir / csv_name, profile)
            for field, values in profile.items():
                combined["{}__{}".format(name, field)] = values
            metadata[name] = {
                "label": label,
                "y_m": y,
                "z_m": z,
                "csv": csv_name,
                "valid_all_times_count": int(np.count_nonzero(common_valid)),
                "sample_count": int(len(x)),
                "fluid_segments_index": [[start, stop - 1]
                                           for start, stop in segments],
                "fluid_segments_x_m": [[float(x[start]), float(x[stop - 1])]
                                         for start, stop in segments],
                "favre_Mach_range": finite_range(favre_mach),
                "mean_pressure_Pa_range": finite_range(pressure_mean),
            }

        np.savez_compressed(output_dir / "axis_profiles_late_mean_combined.npz", **combined)
        primary_stem = output_dir / "axis_profiles_late_mean_Mach_pressure_schlieren"
        render_primary(
            profiles, primary_stem, times, args.nozzle_plane_x,
            args.p_inf, args.body_diameter, args.dpi,
        )
        render_definition_comparison(
            profiles, output_dir / "axis_profiles_late_mean_Mach_definitions.png",
            times, args.dpi,
        )
        summary = {
            "schema_version": 1,
            "steps": steps,
            "times_s": times.tolist(),
            "times_ms": (1.0e3 * times).tolist(),
            "trapezoidal_weights": weights.tolist(),
            "window_s": [float(times[0]), float(times[-1])],
            "window_ms": [float(1.0e3 * times[0]), float(1.0e3 * times[-1])],
            "window_duration_s": float(times[-1] - times[0]),
            "averaging": {
                "ordinary_means": "sum_i w_i q_i",
                "favre_velocity": "sum_i(w_i rho_i velocity_i) / sum_i(w_i rho_i)",
                "favre_temperature": "sum_i(w_i rho_i T_i) / sum_i(w_i rho_i)",
                "favre_Mach": (
                    "magnitude(Favre velocity) / sqrt(gamma*mean(p)/mean(rho))"
                ),
                "mean_instantaneous_Mach": "sum_i(w_i M_i)",
                "validity": (
                    "complete-record only: a sample must pass the strict four-corner "
                    "fluid test in every snapshot"
                ),
                "abs_dmeanrho_dx": (
                    "absolute derivative of time-mean density inside each common-fluid segment"
                ),
                "mean_abs_drhodx": "time mean of the instantaneous segmented derivatives",
            },
            "p_inf_Pa": args.p_inf,
            "body_diameter_m": args.body_diameter,
            "body_centre_x_m": args.body_centre_x,
            "nozzle_exit_plane_x_m": args.nozzle_plane_x,
            "axes": metadata,
            "outputs": {
                "combined_npz": "axis_profiles_late_mean_combined.npz",
                "primary_png": primary_stem.with_suffix(".png").name,
                "primary_pdf": primary_stem.with_suffix(".pdf").name,
                "Mach_definition_comparison_png": (
                    "axis_profiles_late_mean_Mach_definitions.png"
                ),
            },
            "interpretation_warning": (
                "This is a finite-window late-time diagnostic over seven snapshots, "
                "not a demonstrated statistically converged turbulent mean."
            ),
        }
        with (output_dir / "axis_profiles_late_mean_summary.json").open("w") as handle:
            json.dump(summary, handle, indent=2, sort_keys=True)
        print("Wrote late axis means to {}".format(output_dir), flush=True)
    finally:
        for record in records:
            record.close()


if __name__ == "__main__":
    main()
