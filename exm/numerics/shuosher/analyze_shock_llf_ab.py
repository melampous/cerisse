#!/usr/bin/env python3
"""Analyse a same-binary Shu--Osher ``afd_shock_llf`` A/B study."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shlex
from pathlib import Path

import numpy as np

try:
    import yt
except ImportError as exc:  # pragma: no cover - dependency diagnostic
    raise SystemExit(
        "analyze_shock_llf_ab.py requires yt (python3 -m pip install yt)"
    ) from exc


MODES = ("shock_llf_0", "shock_llf_1")
FIELDS = ("Density", "Xmom", "Energy", "pressure", "x_velocity")
GAMMA = 1.4
SMOOTHNESS_THRESHOLD = 0.08
SHOCK_PRESSURE_JUMP = 0.01
SHOCK_COMPRESSION = 0.001
ENTROPY_WINDOW = (-1.5, 1.5)
SHOCK_WINDOW = (2.1, 2.6)


def latest_plotfile(path: Path) -> Path:
    if path.is_dir() and (path / "Header").is_file():
        return path
    pattern = re.compile(r"^plt(\d+)$")
    candidates = [
        item
        for item in path.glob("plt*")
        if pattern.fullmatch(item.name) and (item / "Header").is_file()
    ]
    if not candidates:
        candidates = [
            item
            for item in path.glob("**/plt*")
            if pattern.fullmatch(item.name) and (item / "Header").is_file()
        ]
    if not candidates:
        raise FileNotFoundError(f"no AMReX plotfile found below {path}")
    return max(
        candidates,
        key=lambda item: (int(pattern.fullmatch(item.name)[1]), str(item)),
    )


def field_key(dataset, name: str):
    wanted = name.lower()
    for key in dataset.field_list:
        if key[1].lower() == wanted:
            return key
    raise KeyError(f"field {name!r} is unavailable; fields are {dataset.field_list}")


def read_case(path: Path) -> dict:
    plotfile = latest_plotfile(path)
    dataset = yt.load(str(plotfile))
    if int(dataset.index.max_level) != 0:
        raise ValueError(f"A/B analysis requires a uniform level-0 result: {plotfile}")
    data = dataset.all_data()
    x = np.asarray(data[("index", "x")].d, dtype=float)
    order = np.argsort(x)
    arrays = {
        name: np.asarray(data[field_key(dataset, name)].d, dtype=float)[order]
        for name in FIELDS
    }
    return {
        "plotfile": plotfile,
        "dataset": dataset,
        "x": x[order],
        "fields": arrays,
    }


def advance_seconds(path: Path) -> float:
    log = path / "run.log"
    if not log.is_file():
        return float("nan")
    match = re.search(
        r"Run Time advance\s*=\s*([0-9.eE+-]+)",
        log.read_text(errors="replace"),
    )
    return float(match.group(1)) if match else float("nan")


def final_step(plotfile: Path) -> int:
    match = re.fullmatch(r"plt(\d+)", plotfile.name)
    return int(match[1]) if match else -1


def require_same_grid(first: dict, second: dict, cells: int) -> None:
    x0 = first["x"]
    x1 = second["x"]
    if x0.size != cells or x1.size != cells:
        raise ValueError(
            f"N={cells} result sizes are {x0.size} and {x1.size}, not {cells}"
        )
    if not np.array_equal(x0, x1):
        raise ValueError(f"shock_llf=0 and 1 use different grids at N={cells}")
    time0 = float(first["dataset"].current_time)
    time1 = float(second["dataset"].current_time)
    if time0 != time1:
        raise ValueError(
            f"shock_llf=0 and 1 have different final times at N={cells}: "
            f"{time0} and {time1}"
        )
    for item in (first, second):
        dataset = item["dataset"]
        if np.asarray(dataset.domain_dimensions, dtype=int)[0] != cells:
            raise ValueError(
                f"plotfile {item['plotfile']} reports dimensions "
                f"{dataset.domain_dimensions}, expected N={cells}"
            )


def reference_density(reference: dict, x: np.ndarray) -> np.ndarray:
    x_ref = reference["x"]
    rho_ref = reference["fields"]["Density"]
    if x_ref.size <= x.size:
        raise ValueError("the diagnostic reference must be finer than every A/B case")
    tolerance = 32.0 * np.finfo(float).eps * max(1.0, np.max(np.abs(x_ref)))
    if np.min(x) < np.min(x_ref) - tolerance or np.max(x) > np.max(x_ref) + tolerance:
        raise ValueError("case cell centres extend beyond the reference cell centres")
    return np.interp(x, x_ref, rho_ref)


def norm_row(prefix: str, difference: np.ndarray) -> dict[str, float]:
    return {
        f"{prefix}_L1": float(np.mean(np.abs(difference))),
        f"{prefix}_L2": float(np.sqrt(np.mean(difference * difference))),
        f"{prefix}_Linf": float(np.max(np.abs(difference))),
    }


def selected(values: np.ndarray, x: np.ndarray, window: tuple[float, float]):
    mask = (x >= window[0]) & (x <= window[1])
    if not np.any(mask):
        raise ValueError(f"no cell centres lie in window {window}")
    return values[mask]


def shock_metrics(case: dict) -> dict[str, float]:
    """Return reproducible density-gradient measures of the leading shock."""
    x = case["x"]
    rho = case["fields"]["Density"]
    face_x = 0.5 * (x[:-1] + x[1:])
    gradients = np.abs(np.diff(rho) / np.diff(x))
    mask = (face_x >= SHOCK_WINDOW[0]) & (face_x <= SHOCK_WINDOW[1])
    candidates = np.flatnonzero(mask)
    if candidates.size == 0:
        raise ValueError(f"no density-gradient faces lie in {SHOCK_WINDOW}")
    index = int(candidates[np.argmax(gradients[candidates])])
    shock_x = float(face_x[index])
    dx = float(x[1] - x[0])

    # Median states on either side suppress the post-shock entropy oscillation.
    left_mask = (x >= shock_x - 0.25) & (x <= shock_x - 0.10)
    right_mask = (x >= shock_x + 0.10) & (x <= shock_x + 0.25)
    if not np.any(left_mask) or not np.any(right_mask):
        raise ValueError("insufficient cells for the shock plateau estimates")
    rho_left = float(np.median(rho[left_mask]))
    rho_right = float(np.median(rho[right_mask]))
    jump = abs(rho_left - rho_right)
    maximum_gradient = float(gradients[index])
    equivalent_thickness = (
        jump / maximum_gradient if maximum_gradient > 0.0 else float("nan")
    )
    return {
        "shock_x_max_density_gradient": shock_x,
        "shock_density_gradient_max": maximum_gradient,
        "shock_density_left_median": rho_left,
        "shock_density_right_median": rho_right,
        "shock_density_jump": jump,
        "shock_equivalent_thickness": equivalent_thickness,
        "shock_equivalent_thickness_cells": equivalent_thickness / dx,
    }


def shock_gate_metrics(case: dict) -> dict[str, float | int]:
    """Reproduce the one-dimensional bulk pressure-compression LLF gate."""

    x = case["x"]
    density = case["fields"]["Density"]
    pressure = case["fields"]["pressure"]
    velocity = case["fields"]["x_velocity"]
    sound = np.sqrt(GAMMA * pressure / density)
    tiny = np.finfo(float).tiny
    face_x: list[float] = []
    broad_flags: list[bool] = []
    pressure_flags: list[bool] = []
    compression_flags: list[bool] = []
    fallback_flags: list[bool] = []

    for face in range(3, x.size - 2):
        ids = np.arange(face - 3, face + 3)
        broad = False
        for values, floor in (
            (density, float(np.max(np.abs(density[ids])))),
            (pressure, float(np.max(np.abs(pressure[ids])))),
            (velocity, float(np.max(np.abs(sound[ids])))),
        ):
            for m in range(4):
                a, b, c = values[ids[m : m + 3]]
                denominator = (
                    abs(a) + 2.0 * abs(b) + abs(c) + floor + tiny
                )
                if (
                    abs(a - 2.0 * b + c) / denominator
                    > SMOOTHNESS_THRESHOLD
                ):
                    broad = True

        pressure_jump = any(
            abs(pressure[ids[m + 1]] - pressure[ids[m]])
            / max(pressure[ids[m]], pressure[ids[m + 1]], tiny)
            > SHOCK_PRESSURE_JUMP
            for m in range(5)
        )
        compression = any(
            max(velocity[ids[m] - 1] - velocity[ids[m] + 1], 0.0)
            / (
                2.0
                * max(
                    sound[ids[m] - 1],
                    sound[ids[m]],
                    sound[ids[m] + 1],
                    tiny,
                )
            )
            > SHOCK_COMPRESSION
            for m in range(1, 5)
        )
        face_x.append(float(0.5 * (x[face - 1] + x[face])))
        broad_flags.append(broad)
        pressure_flags.append(pressure_jump)
        compression_flags.append(compression)
        fallback_flags.append(broad and pressure_jump and compression)

    locations = np.asarray(face_x)
    broad_array = np.asarray(broad_flags, dtype=bool)
    pressure_array = np.asarray(pressure_flags, dtype=bool)
    compression_array = np.asarray(compression_flags, dtype=bool)
    fallback_array = np.asarray(fallback_flags, dtype=bool)
    result: dict[str, float | int] = {
        "offline_broad_nonsmooth_faces": int(np.count_nonzero(broad_array)),
        "offline_pressure_jump_faces": int(np.count_nonzero(pressure_array)),
        "offline_compression_faces": int(
            np.count_nonzero(compression_array)
        ),
        "offline_pressure_compression_faces": int(
            np.count_nonzero(pressure_array & compression_array)
        ),
        "offline_llf_fallback_faces": int(
            np.count_nonzero(fallback_array)
        ),
    }
    for prefix, window in (
        ("entropy", ENTROPY_WINDOW),
        ("shock", SHOCK_WINDOW),
    ):
        selected_faces = (
            (locations >= window[0]) & (locations <= window[1])
        )
        result[f"{prefix}_offline_broad_nonsmooth_faces"] = int(
            np.count_nonzero(broad_array & selected_faces)
        )
        result[f"{prefix}_offline_pressure_jump_faces"] = int(
            np.count_nonzero(pressure_array & selected_faces)
        )
        result[f"{prefix}_offline_compression_faces"] = int(
            np.count_nonzero(compression_array & selected_faces)
        )
        result[f"{prefix}_offline_llf_fallback_faces"] = int(
            np.count_nonzero(fallback_array & selected_faces)
        )
    return result


def recorded_command(path: Path) -> list[str]:
    for line in path.read_text().splitlines():
        if line.startswith("command="):
            return shlex.split(line.removeprefix("command="))
    raise ValueError(f"no command= record in {path}")


def normalized_paired_command(command: list[str]) -> list[str]:
    normalized = []
    for token in command:
        if token.startswith("amr.plot_file="):
            normalized.append("amr.plot_file=<variant-output>/plt")
        elif token.startswith("cns.afd_shock_llf="):
            normalized.append("cns.afd_shock_llf=<variant>")
        else:
            normalized.append(token)
    return normalized


def verify_paired_commands(root: Path, cells_list: list[int]) -> dict:
    pair_records = {}
    reference_normalized = None
    for cells in cells_list:
        commands = {
            mode: recorded_command(root / mode / f"N{cells}" / "run_command.txt")
            for mode in MODES
        }
        executable_paths = [command[6] for command in commands.values()]
        executable_hashes = []
        for executable in executable_paths:
            digest = hashlib.sha256(Path(executable).read_bytes()).hexdigest()
            executable_hashes.append(digest)
        if len(set(executable_paths)) != 1 or len(set(executable_hashes)) != 1:
            raise ValueError(f"paired N={cells} runs did not use one executable")
        normalized = {
            mode: normalized_paired_command(command)
            for mode, command in commands.items()
        }
        if normalized[MODES[0]] != normalized[MODES[1]]:
            raise ValueError(
                f"paired N={cells} command lines differ by more than output path "
                "and cns.afd_shock_llf"
            )
        cross_grid = list(normalized[MODES[0]])
        cross_grid = [
            "amr.n_cell=<N>" if token.startswith("amr.n_cell=") else token
            for token in cross_grid
        ]
        if reference_normalized is None:
            reference_normalized = cross_grid
        elif cross_grid != reference_normalized:
            raise ValueError(
                f"N={cells} command differs from the common A/B configuration"
            )
        pair_records[f"N{cells}"] = {
            "same_executable_path": True,
            "same_executable_sha256": True,
            "executable_sha256": executable_hashes[0],
            "only_variant_argument": "cns.afd_shock_llf",
        }
    return {
        "all_pairs_passed": True,
        "all_grids_use_common_configuration": True,
        "pairs": pair_records,
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path, help="root of the completed A/B study")
    parser.add_argument("--cells", nargs="+", type=int, default=[200, 400, 800])
    parser.add_argument(
        "--reference",
        type=Path,
        help="optional fine-grid diagnostic reference, not an exact solution",
    )
    parser.add_argument("--per-case-csv", type=Path)
    parser.add_argument("--pairwise-csv", type=Path)
    parser.add_argument("--figure", type=Path)
    parser.add_argument("--png", type=Path)
    parser.add_argument("--metadata", type=Path)
    args = parser.parse_args()

    per_case_path = args.per_case_csv or args.root / "per_case_metrics.csv"
    pairwise_path = args.pairwise_csv or args.root / "pairwise_metrics.csv"
    figure_path = args.figure or args.root / "density_ab_comparison.pdf"
    png_path = args.png or args.root / "density_ab_comparison.png"
    metadata_path = args.metadata or args.root / "analysis_metadata.json"

    reference = read_case(args.reference) if args.reference else None
    command_invariants = verify_paired_commands(args.root, args.cells)
    cases: dict[tuple[str, int], dict] = {}
    per_case_rows: list[dict] = []
    pairwise_rows: list[dict] = []

    for cells in args.cells:
        for mode in MODES:
            case_dir = args.root / mode / f"N{cells}"
            case = read_case(case_dir)
            cases[(mode, cells)] = case
            rho = case["fields"]["Density"]
            pressure = case["fields"]["pressure"]
            row = {
                "mode": mode,
                "afd_shock_llf": int(mode[-1]),
                "cells": cells,
                "time": float(case["dataset"].current_time),
                "final_step": final_step(case["plotfile"]),
                "advance_seconds": advance_seconds(case_dir),
                "all_fields_finite": all(
                    bool(np.all(np.isfinite(values)))
                    for values in case["fields"].values()
                ),
                "density_min": float(np.min(rho)),
                "density_max": float(np.max(rho)),
                "density_TV": float(np.sum(np.abs(np.diff(rho)))),
                "pressure_min": float(np.min(pressure)),
                "pressure_max": float(np.max(pressure)),
                "afd_smoothness_threshold": SMOOTHNESS_THRESHOLD,
                "afd_shock_pressure_jump": SHOCK_PRESSURE_JUMP,
                "afd_shock_compression": SHOCK_COMPRESSION,
            }
            entropy_density = selected(rho, case["x"], ENTROPY_WINDOW)
            row.update(
                {
                    "entropy_window_lo": ENTROPY_WINDOW[0],
                    "entropy_window_hi": ENTROPY_WINDOW[1],
                    "entropy_density_min": float(np.min(entropy_density)),
                    "entropy_density_max": float(np.max(entropy_density)),
                    "entropy_density_peak_to_peak": float(np.ptp(entropy_density)),
                    "entropy_density_TV": float(
                        np.sum(np.abs(np.diff(entropy_density)))
                    ),
                }
            )
            row.update(shock_metrics(case))
            row.update(shock_gate_metrics(case))
            if reference is not None:
                reference_time = float(reference["dataset"].current_time)
                if float(case["dataset"].current_time) != reference_time:
                    raise ValueError(
                        f"time mismatch between {case['plotfile']} and "
                        f"{reference['plotfile']}"
                    )
                density_error = rho - reference_density(reference, case["x"])
                row.update(norm_row("density_reference", density_error))
                entropy_error = selected(
                    density_error, case["x"], ENTROPY_WINDOW
                )
                row.update(
                    norm_row("entropy_density_reference", entropy_error)
                )
            per_case_rows.append(row)

        off = cases[("shock_llf_0", cells)]
        on = cases[("shock_llf_1", cells)]
        require_same_grid(off, on, cells)
        pair_row: dict[str, float | int] = {
            "cells": cells,
            "time": float(off["dataset"].current_time),
        }
        for field in FIELDS:
            difference = on["fields"][field] - off["fields"][field]
            label = field.lower()
            pair_row.update(norm_row(label, difference))
            scale = max(
                float(np.max(np.abs(off["fields"][field]))),
                np.finfo(float).tiny,
            )
            pair_row[f"{label}_relative_Linf"] = (
                float(np.max(np.abs(difference))) / scale
            )
            pair_row[f"{label}_changed_cells_gt_1e-12"] = int(
                np.count_nonzero(np.abs(difference) > 1.0e-12)
            )
        density_difference = (
            on["fields"]["Density"] - off["fields"]["Density"]
        )
        pair_row.update(
            norm_row(
                "entropy_density",
                selected(density_difference, off["x"], ENTROPY_WINDOW),
            )
        )
        pair_row.update(
            norm_row(
                "shock_density",
                selected(density_difference, off["x"], SHOCK_WINDOW),
            )
        )
        maximum_index = int(np.argmax(np.abs(density_difference)))
        pair_row["density_max_difference_x"] = float(off["x"][maximum_index])
        off_shock = shock_metrics(off)
        on_shock = shock_metrics(on)
        pair_row["shock_position_difference"] = (
            on_shock["shock_x_max_density_gradient"]
            - off_shock["shock_x_max_density_gradient"]
        )
        pair_row["shock_thickness_difference"] = (
            on_shock["shock_equivalent_thickness"]
            - off_shock["shock_equivalent_thickness"]
        )
        pairwise_rows.append(pair_row)

    write_csv(per_case_path, per_case_rows)
    write_csv(pairwise_path, pairwise_rows)

    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(
        len(args.cells),
        2,
        figsize=(7.2, 2.35 * len(args.cells)),
        sharex="col",
        squeeze=False,
        gridspec_kw={"width_ratios": [2.25, 1.0]},
    )
    for row_index, cells in enumerate(args.cells):
        off = cases[("shock_llf_0", cells)]
        on = cases[("shock_llf_1", cells)]
        x = off["x"]
        rho_off = off["fields"]["Density"]
        rho_on = on["fields"]["Density"]
        profile_axis = axes[row_index, 0]
        shock_axis = axes[row_index, 1]
        if reference is not None:
            profile_axis.plot(
                reference["x"],
                reference["fields"]["Density"],
                color="0.72",
                linewidth=0.75,
                label="N=6400 diagnostic reference" if row_index == 0 else None,
            )
        profile_axis.plot(
            x,
            rho_off,
            color="#0072B2",
            linewidth=1.0,
            label=r"$\mathtt{afd\_shock\_llf}=0$" if row_index == 0 else None,
        )
        profile_axis.plot(
            x,
            rho_on,
            color="#D55E00",
            linewidth=0.9,
            linestyle="--",
            label=r"$\mathtt{afd\_shock\_llf}=1$" if row_index == 0 else None,
        )
        shock_axis.plot(
            x,
            rho_off,
            color="#0072B2",
            linewidth=1.0,
        )
        shock_axis.plot(
            x,
            rho_on,
            color="#D55E00",
            linewidth=0.9,
            linestyle="--",
        )
        if reference is not None:
            shock_axis.plot(
                reference["x"],
                reference["fields"]["Density"],
                color="0.72",
                linewidth=0.75,
            )
        off_shock = shock_metrics(off)
        on_shock = shock_metrics(on)
        shock_axis.axvline(
            off_shock["shock_x_max_density_gradient"],
            color="#0072B2",
            linewidth=0.6,
            alpha=0.7,
        )
        shock_axis.axvline(
            on_shock["shock_x_max_density_gradient"],
            color="#D55E00",
            linewidth=0.6,
            linestyle="--",
            alpha=0.7,
        )
        profile_axis.set_ylabel(rf"$\rho$, $N={cells}$")
        shock_axis.set_ylabel(r"$\rho$")
        profile_axis.set_xlim(*ENTROPY_WINDOW)
        shock_axis.set_xlim(*SHOCK_WINDOW)
        profile_axis.grid(alpha=0.18, linewidth=0.4)
        shock_axis.grid(alpha=0.18, linewidth=0.4)
    axes[-1, 0].set_xlabel(r"$x$")
    axes[-1, 1].set_xlabel(r"$x$")
    axes[0, 0].set_title("Post-shock entropy-wave region")
    axes[0, 1].set_title("Leading-shock region")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="upper center",
        ncol=len(labels),
        frameon=False,
        bbox_to_anchor=(0.5, 1.0),
    )
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.955))
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(figure_path, dpi=300)
    figure.savefig(png_path, dpi=300)
    plt.close(figure)

    metadata = {
        "root": str(args.root.resolve()),
        "modes": list(MODES),
        "cells": args.cells,
        "fields": list(FIELDS),
        "entropy_window": list(ENTROPY_WINDOW),
        "shock_window": list(SHOCK_WINDOW),
        "shock_position_definition": (
            "Face location of the maximum absolute density gradient in the "
            "shock window."
        ),
        "shock_thickness_definition": (
            "Median density jump across the shock divided by the maximum "
            "absolute density gradient."
        ),
        "offline_shock_gate": {
            "broad_smoothness_threshold": SMOOTHNESS_THRESHOLD,
            "relative_pressure_jump_threshold": SHOCK_PRESSURE_JUMP,
            "normal_compression_threshold": SHOCK_COMPRESSION,
            "semantics": (
                "The LLF fallback requires the broad curvature sensor, "
                "the pressure-jump gate, and the compression gate."
            ),
        },
        "paired_command_invariants": command_invariants,
        "reference": str(latest_plotfile(args.reference).resolve())
        if args.reference
        else None,
        "reference_semantics": (
            "Fine-grid numerical profile used only as a common diagnostic; "
            "it is not an exact solution."
            if args.reference
            else None
        ),
        "per_case_csv": str(per_case_path.resolve()),
        "pairwise_csv": str(pairwise_path.resolve()),
        "figure": str(figure_path.resolve()),
        "png": str(png_path.resolve()),
    }
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")

    for row in pairwise_rows:
        print(
            f"N={int(row['cells']):4d}  "
            f"density L1={row['density_L1']:.6e}  "
            f"L2={row['density_L2']:.6e}  "
            f"Linf={row['density_Linf']:.6e}  "
            f"x(max)={row['density_max_difference_x']:.6f}"
        )
    print(f"wrote {per_case_path}")
    print(f"wrote {pairwise_path}")
    print(f"wrote {figure_path}")
    print(f"wrote {png_path}")
    print(f"wrote {metadata_path}")


if __name__ == "__main__":
    yt.set_log_level(40)
    main()
