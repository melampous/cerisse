#!/usr/bin/env python3
"""Audit and integrate pressure from Cerisse IBM surface VTP files.

The VTP connectivity orientation defines the body-outward normal.  The force
reported here is the pressure force on the body,

    F_p = - integral_S (p - p_ref) n dS.

Two-node cells are treated as 2-D edges per unit span; cells with three or more
nodes are triangulated about their first vertex.  Passing two files also emits
an area-weighted pressure A/B comparison when their meshes are identical.
"""

from __future__ import annotations

import argparse
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


def read_vtp(path: Path) -> dict[str, object]:
    piece = ET.parse(path).getroot().find(".//Piece")
    if piece is None:
        raise ValueError(f"{path}: missing PolyData/Piece")

    points = numbers(named_array(piece.find("Points"), "Points")).reshape(-1, 3)
    polys = piece.find("Polys")
    connectivity = numbers(named_array(polys, "connectivity"), np.int64)
    offsets = numbers(named_array(polys, "offsets"), np.int64)
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
    }


def oriented_measures(surface: dict[str, object]) -> tuple[np.ndarray, np.ndarray]:
    points = np.asarray(surface["points"])
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
    }


def analyze(surface: dict[str, object], args: argparse.Namespace) -> dict[str, object]:
    vectors, _ = oriented_measures(surface)
    area = np.linalg.norm(vectors, axis=1)
    pressure = np.asarray(surface["fields"]["Pressure"], dtype=float)
    selected = selected_cells(surface, args)
    selected_pressure = pressure[selected]
    selected_area = area[selected]
    selected_vectors = vectors[selected]
    finite = np.isfinite(selected_pressure)
    pressure_force = -np.sum(
        (selected_pressure - args.p_ref)[:, None] * selected_vectors, axis=0
    )
    total_measure = float(np.sum(selected_area))
    closure = np.sum(selected_vectors, axis=0)

    result: dict[str, object] = {
        "file": str(surface["path"]),
        "cells": int(pressure.size),
        "selected_cells": int(np.count_nonzero(selected)),
        "total_measure": total_measure,
        "normal_closure_vector": closure.tolist(),
        "normal_closure_relative": float(np.linalg.norm(closure) / total_measure),
        "pressure": {
            "finite": int(np.count_nonzero(finite)),
            "nonfinite": int(np.count_nonzero(~finite)),
            "nonpositive": int(np.count_nonzero(finite & (selected_pressure <= 0.0))),
            "min": float(np.min(selected_pressure[finite])) if finite.any() else math.nan,
            "max": float(np.max(selected_pressure[finite])) if finite.any() else math.nan,
            "area_weighted_mean": float(
                np.sum(selected_area * selected_pressure) / np.sum(selected_area)
            ),
        },
        "pressure_force": pressure_force.tolist(),
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
            diagnostics["fallback_fraction"] = fallback_summary
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
    return result


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
    parser.add_argument("--include-invalid-ip", action="store_true",
                        help="Include faces with IP_quality < 0 (excluded by default).")
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()

    surfaces = [read_vtp(path) for path in args.files]
    results = [analyze(surface, args) for surface in surfaces]
    report: dict[str, object] = {"surfaces": results}
    if len(surfaces) == 2:
        report["second_minus_first"] = compare(
            surfaces[0], surfaces[1], results[0], results[1], args
        )
    elif len(surfaces) > 2:
        report["comparisons"] = [
            compare(surfaces[0], surface, results[0], result, args)
            for surface, result in zip(surfaces[1:], results[1:])
        ]

    output = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.write_text(output)
    print(output, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
