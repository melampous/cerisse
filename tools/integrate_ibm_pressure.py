#!/usr/bin/env python3
"""Audit and integrate pressure from Cerisse IBM surface VTP files.

The VTP connectivity orientation defines the body-outward normal.  The force
reported here is the pressure force on the body,

    F_p = - integral_S (p - p_ref) n dS.

Two-node cells are treated as 2-D edges per unit span; cells with three or more
nodes are triangulated about their first vertex.  Multiple files also emit
area-weighted pressure comparisons when their meshes are identical.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np


def numbers(node: ET.Element | None, dtype=float) -> np.ndarray:
    if node is None:
        raise ValueError("missing VTP DataArray")
    return np.fromstring(node.text or "", sep=" ", dtype=dtype)


def named_array(parent: ET.Element | None, name: str) -> ET.Element:
    if parent is None:
        raise ValueError(f"missing VTP section containing {name}")
    for node in parent.findall("DataArray"):
        if node.attrib.get("Name") == name:
            return node
    raise ValueError(f"missing VTP DataArray {name}")


def points_array(parent: ET.Element | None) -> ET.Element:
    """Return the VTK points array, accepting legacy unnamed arrays."""
    if parent is None:
        raise ValueError("missing VTP Points section")
    arrays = parent.findall("DataArray")
    for node in arrays:
        if node.attrib.get("Name") == "Points":
            return node
    if len(arrays) == 1:
        return arrays[0]
    raise ValueError("missing unambiguous VTP points DataArray")


def read_vtp(path: Path) -> dict[str, object]:
    piece = ET.parse(path).getroot().find(".//Piece")
    if piece is None:
        raise ValueError(f"{path}: missing PolyData/Piece")

    points = numbers(points_array(piece.find("Points"))).reshape(-1, 3)
    topology = piece.find("Polys")
    if topology is None:
        topology = piece.find("Lines")
    if topology is None:
        topology = piece.find("Verts")
    connectivity = numbers(named_array(topology, "connectivity"), np.int64)
    offsets = numbers(named_array(topology, "offsets"), np.int64)
    starts = np.concatenate((np.array([0], dtype=np.int64), offsets[:-1]))
    cells = [connectivity[start:stop] for start, stop in zip(starts, offsets)]

    fields: dict[str, np.ndarray] = {}
    cell_data = piece.find("CellData")
    if cell_data is not None:
        for node in cell_data.findall("DataArray"):
            name = node.attrib.get("Name")
            if not name:
                continue
            dtype = np.int64 if node.attrib.get("type", "").startswith("Int") else float
            fields[name] = numbers(node, dtype)

    metadata: dict[str, float | list[float]] = {}
    field_data = piece.find("FieldData")
    if field_data is not None:
        for node in field_data.findall("DataArray"):
            name = node.attrib.get("Name")
            if not name:
                continue
            value = numbers(node)
            if value.size == 1:
                metadata[name] = float(value[0])
            elif value.size:
                metadata[name] = value.tolist()

    if "Pressure" not in fields:
        raise ValueError(f"{path}: missing CellData/Pressure")
    if len(cells) != fields["Pressure"].size:
        raise ValueError(
            f"{path}: {len(cells)} cells but {fields['Pressure'].size} pressures"
        )
    return {
        "path": path,
        "points": points,
        "cells": cells,
        "connectivity": connectivity,
        "offsets": offsets,
        "fields": fields,
        "metadata": metadata,
    }


def oriented_measures(surface: dict[str, object]) -> tuple[np.ndarray, np.ndarray]:
    points = np.asarray(surface["points"])
    fields = surface["fields"]
    direct_normal = fields.get("FluidOutwardNormal")
    direct_measure = fields.get("QuadratureMeasure")
    if direct_normal is not None and direct_measure is not None:
        fluid_normal = np.asarray(direct_normal, dtype=float).reshape(-1, 3)
        measure = np.asarray(direct_measure, dtype=float)
        if fluid_normal.shape[0] != measure.size:
            raise ValueError("wall quadrature normal/measure size mismatch")
        if points.shape[0] != measure.size:
            raise ValueError("wall quadrature point/measure size mismatch")
        # The generic integration path uses body-outward area vectors and
        # applies -p*n.  Cut-control VTP stores the opposite convention:
        # the outward normal of the fluid control volume.
        return -fluid_normal * measure[:, None], points
    vectors = []
    centroids = []
    for cell in surface["cells"]:
        vertices = points[np.asarray(cell, dtype=np.int64)]
        if vertices.shape[0] == 2:
            edge = vertices[1] - vertices[0]
            vectors.append(np.array([edge[1], -edge[0], 0.0]))
        elif vertices.shape[0] >= 3:
            area_vector = np.zeros(3)
            for index in range(1, vertices.shape[0] - 1):
                area_vector += 0.5 * np.cross(
                    vertices[index] - vertices[0],
                    vertices[index + 1] - vertices[0],
                )
            vectors.append(area_vector)
        else:
            raise ValueError("surface cell has fewer than two vertices")
        centroids.append(vertices.mean(axis=0))
    return np.asarray(vectors), np.asarray(centroids)


def weighted_l2(values: np.ndarray, weight: np.ndarray) -> float:
    return math.sqrt(float(np.sum(weight * values * values) / np.sum(weight)))


def selected_cells(surface: dict[str, object], args: argparse.Namespace) -> np.ndarray:
    pressure = np.asarray(surface["fields"]["Pressure"])
    selected = np.ones(pressure.size, dtype=bool)
    quality = surface["fields"].get("IP_quality")
    if quality is not None and not args.include_invalid_ip:
        selected &= np.asarray(quality) >= 0
    return selected


def diagnostic_summary(
    values: np.ndarray, area: np.ndarray, selected: np.ndarray
) -> dict[str, object]:
    valid = selected & np.isfinite(values) & (values >= 0.0)
    if not np.any(valid):
        return {"applicable_cells": 0}
    selected_values = values[valid]
    selected_area = area[valid]
    return {
        "applicable_cells": int(np.count_nonzero(valid)),
        "min": float(np.min(selected_values)),
        "max": float(np.max(selected_values)),
        "area_weighted_mean": float(
            np.sum(selected_area * selected_values) / np.sum(selected_area)
        ),
        "applicable_measure": float(np.sum(selected_area)),
    }


def analyze(surface: dict[str, object], args: argparse.Namespace) -> dict[str, object]:
    vectors, centroids = oriented_measures(surface)
    area = np.linalg.norm(vectors, axis=1)
    pressure = np.asarray(surface["fields"]["Pressure"], dtype=float)
    selected = selected_cells(surface, args)
    finite_geometry = (
        np.isfinite(vectors).all(axis=1)
        & np.isfinite(centroids).all(axis=1)
        & np.isfinite(area)
    )
    positive_measure = finite_geometry & (area > 0.0)
    invalid_geometry_count = int(np.count_nonzero(~finite_geometry))
    nonpositive_measure_count = int(
        np.count_nonzero(finite_geometry & (area <= 0.0))
    )
    selected &= positive_measure
    selected_pressure = pressure[selected]
    selected_area = area[selected]
    selected_vectors = vectors[selected]
    finite = np.isfinite(selected_pressure)
    if not np.any(selected):
        raise ValueError(f"{surface['path']}: no selected surface cells")
    finite_area = selected_area[finite]
    finite_pressure = selected_pressure[finite]
    finite_vectors = selected_vectors[finite]
    pressure_force = -np.sum(
        (finite_pressure - args.p_ref)[:, None] * finite_vectors, axis=0
    )
    total_measure = float(np.sum(area[positive_measure]))
    selected_measure = float(np.sum(selected_area))
    closure = np.sum(vectors[positive_measure], axis=0)
    failures: list[str] = []
    nonfinite_count = int(np.count_nonzero(~finite))
    nonpositive_count = int(np.count_nonzero(finite & (selected_pressure <= 0.0)))
    closure_relative = float(
        np.linalg.norm(closure) / max(total_measure, np.finfo(float).tiny)
    )
    if nonfinite_count:
        failures.append("selected surface pressure contains nonfinite values")
    if invalid_geometry_count:
        failures.append("surface contains nonfinite geometry")
    if nonpositive_measure_count:
        failures.append("surface contains non-positive element measure")
    if nonpositive_count and not args.allow_nonpositive_pressure:
        failures.append("selected surface pressure contains nonpositive values")
    if closure_relative > args.normal_closure_tolerance:
        failures.append("surface-normal closure tolerance exceeded")

    result: dict[str, object] = {
        "file": str(surface["path"]),
        "cells": int(pressure.size),
        "selected_cells": int(np.count_nonzero(selected)),
        "total_measure": total_measure,
        "selected_measure": selected_measure,
        "geometry": {
            "nonfinite_cells": invalid_geometry_count,
            "nonpositive_measure_cells": nonpositive_measure_count,
        },
        "normal_closure_vector": closure.tolist(),
        "normal_closure_relative": closure_relative,
        "pressure": {
            "finite": int(np.count_nonzero(finite)),
            "nonfinite": nonfinite_count,
            "nonpositive": nonpositive_count,
            "min": float(np.min(finite_pressure)) if finite.any() else math.nan,
            "max": float(np.max(finite_pressure)) if finite.any() else math.nan,
            "area_weighted_mean": float(
                np.sum(finite_area * finite_pressure) / np.sum(finite_area)
            ) if finite.any() else math.nan,
        },
        "pressure_force": pressure_force.tolist(),
        "audit_failures": failures,
    }

    finite_selected = selected & np.isfinite(pressure)
    if np.any(finite_selected):
        finite_pressure = pressure[finite_selected]
        finite_area = area[finite_selected]
        finite_centroids = centroids[finite_selected]
        stagnation_index = int(np.argmax(finite_pressure))
        below_reference = finite_pressure < args.p_ref
        result["stagnation_pressure"] = {
            "value": float(finite_pressure[stagnation_index]),
            "centroid": finite_centroids[stagnation_index].tolist(),
            "ratio_to_reference": (
                float(finite_pressure[stagnation_index] / args.p_ref)
                if args.p_ref != 0.0 else None
            ),
        }
        result["pressure_below_reference"] = {
            "minimum_minus_reference": float(np.min(finite_pressure - args.p_ref)),
            "cells_below_reference": int(np.count_nonzero(below_reference)),
            "cell_fraction": float(np.mean(below_reference)),
            "measure_fraction": float(
                np.sum(finite_area[below_reference]) / np.sum(finite_area)
            ),
        }

    quality = surface["fields"].get("IP_quality")
    if quality is not None:
        values, counts = np.unique(quality.astype(int), return_counts=True)
        result["ip_quality_counts"] = {
            str(int(value)): int(count) for value, count in zip(values, counts)
        }
        result["invalid_ip_cells_excluded"] = int(
            np.count_nonzero(np.asarray(quality) < 0) if not args.include_invalid_ip else 0
        )
        invalid_fraction = float(np.mean(np.asarray(quality) < 0))
        result["invalid_ip_cell_fraction"] = invalid_fraction
        if (args.max_invalid_ip_fraction is not None and
                invalid_fraction > args.max_invalid_ip_fraction):
            failures.append("invalid image-point fraction exceeds threshold")

    fit_order = surface["fields"].get("IP_fit_order")
    if fit_order is not None:
        fit_values = np.asarray(fit_order, dtype=int)
        values, counts = np.unique(fit_values, return_counts=True)
        applicable = selected & (fit_values >= 0)
        applicable_measure = float(np.sum(area[applicable]))
        measure_fraction = {
            str(order): (
                float(np.sum(area[applicable & (fit_values == order)])) /
                applicable_measure
                if applicable_measure > 0.0 else None
            )
            for order in (0, 1, 2)
        }
        result["ip_fit_order"] = {
            "counts": {
                str(int(value)): int(count)
                for value, count in zip(values, counts)
            },
            "measure_fraction": measure_fraction,
            "constant_fit_measure_fraction": measure_fraction["0"],
        }
        constant_fraction = measure_fraction["0"]
        if (args.max_constant_fit_measure_fraction is not None and
                constant_fraction is not None and
                constant_fraction > args.max_constant_fit_measure_fraction):
            failures.append("constant image-point fit measure fraction exceeds threshold")

    shock = surface["fields"].get("PressureShockSensor")
    legacy_compression = surface["fields"].get("PressureCompressionSensor")
    legacy_jump = surface["fields"].get("PressureJumpSensor")
    fallback = surface["fields"].get("PressureFallbackFraction")
    if (shock is not None or legacy_compression is not None or
            legacy_jump is not None or fallback is not None):
        diagnostics: dict[str, object] = {}
        if shock is not None:
            diagnostics["shock_sensor"] = diagnostic_summary(
                np.asarray(shock, dtype=float), area, selected
            )
        elif legacy_compression is not None:
            diagnostics["legacy_compression_sensor"] = diagnostic_summary(
                np.asarray(legacy_compression, dtype=float), area, selected
            )
        elif legacy_jump is not None:
            diagnostics["legacy_symmetric_jump_sensor"] = diagnostic_summary(
                np.asarray(legacy_jump, dtype=float), area, selected
            )
        if fallback is not None:
            fallback_values = np.asarray(fallback, dtype=float)
            fallback_summary = diagnostic_summary(fallback_values, area, selected)
            applicable = selected & np.isfinite(fallback_values) & (fallback_values >= 0.0)
            fallback_summary["triggered_cells"] = int(
                np.count_nonzero(applicable & (fallback_values > 1.0e-12))
            )
            fallback_summary["full_fallback_cells"] = int(
                np.count_nonzero(applicable & (fallback_values >= 1.0 - 1.0e-12))
            )
            applicable_measure = float(np.sum(area[applicable]))
            triggered = applicable & (fallback_values > 1.0e-12)
            full = applicable & (fallback_values >= 1.0 - 1.0e-12)
            fallback_summary["triggered_cell_fraction"] = (
                float(np.count_nonzero(triggered) / np.count_nonzero(applicable))
                if np.any(applicable) else None
            )
            fallback_summary["triggered_measure_fraction"] = (
                float(np.sum(area[triggered]) / applicable_measure)
                if applicable_measure > 0.0 else None
            )
            fallback_summary["full_fallback_measure_fraction"] = (
                float(np.sum(area[full]) / applicable_measure)
                if applicable_measure > 0.0 else None
            )
            diagnostics["fallback_fraction"] = fallback_summary
            triggered_measure_fraction = fallback_summary[
                "triggered_measure_fraction"
            ]
            if (args.max_fallback_measure_fraction is not None and
                    triggered_measure_fraction is not None and
                    triggered_measure_fraction > args.max_fallback_measure_fraction):
                failures.append("pressure fallback measure fraction exceeds threshold")
        result["pressure_closure_diagnostics"] = diagnostics

    if args.reference_area is not None:
        if args.rho_ref is None or args.u_ref is None:
            raise ValueError("--reference-area requires --rho-ref and --u-ref")
        dynamic_pressure = 0.5 * args.rho_ref * args.u_ref**2
        result["dynamic_pressure"] = dynamic_pressure
        result["reference_area"] = args.reference_area
        result["pressure_force_coefficient"] = (
            pressure_force / (dynamic_pressure * args.reference_area)
        ).tolist()
        cp = (pressure - args.p_ref) / dynamic_pressure
        selected_cp = cp[finite_selected]
        selected_cp_area = area[finite_selected]
        if selected_cp.size == 0:
            raise ValueError(f"{surface['path']}: no finite pressure for Cp")
        cp_stagnation_index = int(np.argmax(selected_cp))
        normals = vectors[finite_selected] / area[finite_selected, None]
        flow_direction = unit_vector(args.flow_direction, "--flow-direction")
        upstream_alignment = -normals @ flow_direction
        max_alignment = float(np.max(upstream_alignment))
        upstream_cap = upstream_alignment >= max_alignment - args.upstream_cap_width
        cap_area = selected_cp_area[upstream_cap]
        cap_pressure = pressure[finite_selected][upstream_cap]
        cap_cp = selected_cp[upstream_cap]
        result["cp"] = {
            "min": float(np.min(selected_cp)),
            "max": float(np.max(selected_cp)),
            "area_weighted_mean": float(
                np.sum(selected_cp_area * selected_cp) / np.sum(selected_cp_area)
            ),
            "stagnation": float(selected_cp[cp_stagnation_index]),
            "stagnation_centroid": centroids[finite_selected][cp_stagnation_index].tolist(),
            "undershoot_below_freestream": float(min(0.0, np.min(selected_cp))),
        }
        result["upstream_stagnation_cap"] = {
            "alignment_max": max_alignment,
            "alignment_width": args.upstream_cap_width,
            "cells": int(np.count_nonzero(upstream_cap)),
            "measure": float(np.sum(cap_area)),
            "pressure_area_weighted_mean": float(
                np.sum(cap_area * cap_pressure) / np.sum(cap_area)
            ),
            "pressure_max": float(np.max(cap_pressure)),
            "cp_area_weighted_mean": float(
                np.sum(cap_area * cap_cp) / np.sum(cap_area)
            ),
            "cp_max": float(np.max(cap_cp)),
        }

        drag_direction = flow_direction
        lift_direction = unit_vector(args.lift_direction, "--lift-direction")
        force_coefficient = np.asarray(result["pressure_force_coefficient"])
        result["pressure_coefficients"] = {
            "drag": float(np.dot(force_coefficient, drag_direction)),
            "lift": float(np.dot(force_coefficient, lift_direction)),
            "absolute_symmetry_lift_error": float(
                abs(np.dot(force_coefficient, lift_direction))
            ),
        }
    result["passed"] = not failures
    return result


def unit_vector(values: list[float], option: str) -> np.ndarray:
    vector = np.asarray(values, dtype=float)
    magnitude = float(np.linalg.norm(vector))
    if not np.isfinite(magnitude) or magnitude <= 0.0:
        raise ValueError(f"{option} must define a finite nonzero vector")
    return vector / magnitude


def write_surface_profile(
    surface: dict[str, object], args: argparse.Namespace, output: Path
) -> None:
    vectors, centroids = oriented_measures(surface)
    area = np.linalg.norm(vectors, axis=1)
    normals = vectors / np.maximum(area[:, None], np.finfo(float).tiny)
    pressure = np.asarray(surface["fields"]["Pressure"], dtype=float)
    selected = selected_cells(surface, args)

    dynamic_pressure = None
    if args.rho_ref is not None and args.u_ref is not None:
        dynamic_pressure = 0.5 * args.rho_ref * args.u_ref**2

    fields = surface["fields"]
    shock = fields.get("PressureShockSensor")
    fallback = fields.get("PressureFallbackFraction")
    quality = fields.get("IP_quality")
    fit_order = fields.get("IP_fit_order")
    level = fields.get("Level")

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow([
            "cell", "selected", "x", "y", "z", "nx", "ny", "nz",
            "measure", "pressure", "cp", "pressure_shock_sensor",
            "pressure_fallback_fraction", "ip_quality", "ip_fit_order", "level",
        ])
        for index in range(pressure.size):
            cp = (
                (pressure[index] - args.p_ref) / dynamic_pressure
                if dynamic_pressure is not None and dynamic_pressure > 0.0
                else math.nan
            )
            writer.writerow([
                index,
                int(selected[index]),
                *centroids[index].tolist(),
                *normals[index].tolist(),
                area[index],
                pressure[index],
                cp,
                shock[index] if shock is not None else math.nan,
                fallback[index] if fallback is not None else math.nan,
                quality[index] if quality is not None else -1,
                fit_order[index] if fit_order is not None else -1,
                level[index] if level is not None else -1,
            ])


def same_mesh(first: dict[str, object], second: dict[str, object]) -> bool:
    return (
        np.array_equal(first["connectivity"], second["connectivity"])
        and np.array_equal(first["offsets"], second["offsets"])
        and np.allclose(first["points"], second["points"], rtol=0.0, atol=1.0e-14)
    )


def compare(
    first: dict[str, object], second: dict[str, object],
    first_result: dict[str, object], second_result: dict[str, object],
    args: argparse.Namespace,
) -> dict[str, object]:
    if not same_mesh(first, second):
        return {"comparable": False, "reason": "surface meshes differ"}
    vectors, _ = oriented_measures(first)
    area = np.linalg.norm(vectors, axis=1)
    p0 = np.asarray(first["fields"]["Pressure"], dtype=float)
    p1 = np.asarray(second["fields"]["Pressure"], dtype=float)
    selected = selected_cells(first, args) & selected_cells(second, args)
    vectors = vectors[selected]
    area = area[selected]
    p0 = p0[selected]
    p1 = p1[selected]
    delta = p1 - p0
    force0 = np.asarray(first_result["pressure_force"])
    force1 = np.asarray(second_result["pressure_force"])
    result: dict[str, object] = {
        "comparable": True,
        "selected_cells": int(np.count_nonzero(selected)),
        "pressure_delta_area_weighted_L2": weighted_l2(delta, area),
        "pressure_delta_Linf": float(np.max(np.abs(delta))),
        "pressure_relative_L2": weighted_l2(delta, area) / weighted_l2(p0, area),
        "pressure_relative_Linf": float(np.max(np.abs(delta)) / np.max(np.abs(p0))),
        "pressure_force_delta": (force1 - force0).tolist(),
    }
    if "pressure_force_coefficient" in first_result:
        coeff0 = np.asarray(first_result["pressure_force_coefficient"])
        coeff1 = np.asarray(second_result["pressure_force_coefficient"])
        result["pressure_force_coefficient_delta"] = (coeff1 - coeff0).tolist()
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="+", type=Path)
    parser.add_argument("--p-ref", type=float, default=0.0)
    parser.add_argument("--rho-ref", type=float)
    parser.add_argument("--u-ref", type=float)
    parser.add_argument("--reference-area", type=float,
                        help="Reference length in 2-D or area in 3-D.")
    parser.add_argument(
        "--flow-direction", type=float, nargs=3, default=[1.0, 0.0, 0.0],
        metavar=("DX", "DY", "DZ"),
        help="Freestream direction used to project pressure drag (default: +x).",
    )
    parser.add_argument(
        "--lift-direction", type=float, nargs=3, default=[0.0, 1.0, 0.0],
        metavar=("LX", "LY", "LZ"),
        help="Symmetry/lift direction used to project pressure lift (default: +y).",
    )
    parser.add_argument("--include-invalid-ip", action="store_true",
                        help="Include faces with IP_quality < 0 (excluded by default).")
    parser.add_argument("--allow-nonpositive-pressure", action="store_true")
    parser.add_argument("--normal-closure-tolerance", type=float, default=1.0e-6)
    parser.add_argument("--max-invalid-ip-fraction", type=float)
    parser.add_argument("--max-fallback-measure-fraction", type=float)
    parser.add_argument("--max-constant-fit-measure-fraction", type=float)
    parser.add_argument(
        "--upstream-cap-width", type=float, default=0.02,
        help=(
            "Include faces whose opposing-flow normal alignment is within this "
            "amount of the surface maximum (default: 0.02)."
        ),
    )
    parser.add_argument("--require-same-mesh", action="store_true",
                        help="Fail unless every VTP uses the identical fixed surface mesh.")
    parser.add_argument("--labels", nargs="*",
                        help="Optional labels matching the input VTP files.")
    parser.add_argument(
        "--profile", type=Path,
        help="Write a per-surface-cell CSV profile (requires exactly one VTP).",
    )
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()

    if args.profile is not None and len(args.files) != 1:
        parser.error("--profile requires exactly one input VTP")
    if args.labels is not None and len(args.labels) != len(args.files):
        parser.error("--labels count must match input files")
    if args.reference_area is not None and args.reference_area <= 0.0:
        parser.error("--reference-area must be positive")
    if args.rho_ref is not None and args.rho_ref <= 0.0:
        parser.error("--rho-ref must be positive")
    if args.u_ref is not None and args.u_ref <= 0.0:
        parser.error("--u-ref must be positive")
    if not 0.0 < args.upstream_cap_width <= 2.0:
        parser.error("--upstream-cap-width must satisfy 0 < value <= 2")
    for option in (
        "max_invalid_ip_fraction", "max_fallback_measure_fraction",
        "max_constant_fit_measure_fraction",
    ):
        value = getattr(args, option)
        if value is not None and not 0.0 <= value <= 1.0:
            parser.error(
                f"--{option.replace('_', '-')} must satisfy 0 <= value <= 1"
            )

    surfaces = [read_vtp(path) for path in args.files]
    results = [analyze(surface, args) for surface in surfaces]
    labels = args.labels or [path.stem for path in args.files]
    for label, result in zip(labels, results):
        result["label"] = label
    if args.profile is not None:
        write_surface_profile(surfaces[0], args, args.profile)
    mesh_matches_first = [same_mesh(surfaces[0], surface) for surface in surfaces]
    mesh_contract = {
        "matches_first": mesh_matches_first,
        "all_identical": all(mesh_matches_first),
    }
    if args.require_same_mesh and not mesh_contract["all_identical"]:
        for result, matches in zip(results, mesh_matches_first):
            if not matches:
                result["audit_failures"].append("surface mesh differs from first file")
                result["passed"] = False
    report: dict[str, object] = {
        "surfaces": results,
        "fixed_mesh_contract": mesh_contract,
    }
    if len(surfaces) == 2:
        report["second_minus_first"] = compare(
            surfaces[0], surfaces[1], results[0], results[1], args
        )
    elif len(surfaces) > 2:
        report["each_minus_finest"] = [
            compare(surface, surfaces[-1], result, results[-1], args)
            for surface, result in zip(surfaces[:-1], results[:-1])
        ]

    report["passed"] = all(result["passed"] for result in results)

    output = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.write_text(output)
    print(output, end="")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
