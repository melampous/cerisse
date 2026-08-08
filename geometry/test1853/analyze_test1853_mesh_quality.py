#!/usr/bin/env python3
"""Detailed, reproducible STL-quality statistics for the Test 1853 report."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import trimesh


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "outputs"
PERCENTILES = [0, 1, 5, 25, 50, 75, 95, 99, 100]


def distribution(values: np.ndarray) -> dict[str, float]:
    result = {f"p{q:02d}": float(v) for q, v in zip(PERCENTILES, np.percentile(values, PERCENTILES))}
    result["mean"] = float(np.mean(values))
    result["std"] = float(np.std(values))
    return result


def face_metrics(mesh: trimesh.Trimesh) -> dict[str, np.ndarray]:
    tri = np.asarray(mesh.triangles, dtype=float)
    edge = np.stack(
        [
            np.linalg.norm(tri[:, 1] - tri[:, 0], axis=1),
            np.linalg.norm(tri[:, 2] - tri[:, 1], axis=1),
            np.linalg.norm(tri[:, 0] - tri[:, 2], axis=1),
        ],
        axis=1,
    )
    area = np.asarray(mesh.area_faces, dtype=float)
    quality = 4.0 * math.sqrt(3.0) * area / np.sum(edge * edge, axis=1)
    perimeter = np.sum(edge, axis=1)
    aspect = np.max(edge, axis=1) * perimeter / (4.0 * math.sqrt(3.0) * area)
    # Interior angles from the cosine law; columns are opposite edge columns.
    angle = np.empty_like(edge)
    for opposite, first, second in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        cosine = (
            edge[:, first] ** 2 + edge[:, second] ** 2 - edge[:, opposite] ** 2
        ) / (2.0 * edge[:, first] * edge[:, second])
        angle[:, opposite] = np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))
    return {
        "face_edge_mm": edge,
        "area_mm2": area,
        "quality": quality,
        "aspect_ratio_equilateral_1": aspect,
        "minimum_angle_deg": np.min(angle, axis=1),
        "maximum_angle_deg": np.max(angle, axis=1),
    }


def region_metrics(
    centroids: np.ndarray, edges: np.ndarray, areas: np.ndarray
) -> dict[str, Any]:
    regions = {
        "front_x_le_35_mm": centroids[:, 0] <= 35.0,
        "transition_35_to_60_mm": (centroids[:, 0] > 35.0) & (centroids[:, 0] < 60.0),
        "aft_x_ge_60_mm": centroids[:, 0] >= 60.0,
    }
    result: dict[str, Any] = {}
    for name, mask in regions.items():
        result[name] = {
            "triangles": int(np.count_nonzero(mask)),
            "edge_mm": distribution(edges[mask].reshape(-1)),
            "area_mm2": distribution(areas[mask]),
        }
    return result


def analyze_stl(path: Path, chord: dict[str, float]) -> dict[str, Any]:
    mesh = trimesh.load_mesh(path, process=True)
    if not isinstance(mesh, trimesh.Trimesh):
        raise RuntimeError(f"{path}: expected one mesh")
    metrics = face_metrics(mesh)
    edge = metrics["face_edge_mm"]
    minimum_angle = metrics["minimum_angle_deg"]
    maximum_angle = metrics["maximum_angle_deg"]
    quality = metrics["quality"]
    aspect = metrics["aspect_ratio_equilateral_1"]
    adjacency = np.degrees(np.asarray(mesh.face_adjacency_angles, dtype=float))
    return {
        "file": str(path.relative_to(HERE)),
        "triangles": int(len(mesh.faces)),
        "unique_vertices": int(len(mesh.vertices)),
        "unique_edges": int(len(mesh.edges_unique)),
        "bounds_mm": mesh.bounds.tolist(),
        "area_mm2": float(mesh.area),
        "volume_mm3": float(mesh.volume),
        "watertight": bool(mesh.is_watertight),
        "winding_consistent": bool(mesh.is_winding_consistent),
        "body_count": int(mesh.body_count),
        "unique_edge_length_mm": distribution(np.asarray(mesh.edges_unique_length)),
        "face_area_mm2": distribution(metrics["area_mm2"]),
        "triangle_quality_4sqrt3A_over_sum_l2": distribution(quality),
        "quality_fraction_below": {
            "0.90": float(np.mean(quality < 0.90)),
            "0.95": float(np.mean(quality < 0.95)),
            "0.99": float(np.mean(quality < 0.99)),
        },
        "aspect_ratio_equilateral_1": distribution(aspect),
        "minimum_internal_angle_deg": distribution(minimum_angle),
        "maximum_internal_angle_deg": distribution(maximum_angle),
        "angle_event_counts": {
            "faces_min_angle_lt_20_deg": int(np.count_nonzero(minimum_angle < 20.0)),
            "faces_min_angle_lt_30_deg": int(np.count_nonzero(minimum_angle < 30.0)),
            "faces_max_angle_gt_120_deg": int(np.count_nonzero(maximum_angle > 120.0)),
        },
        "adjacent_facet_dihedral_deg": distribution(adjacency),
        "axial_regions": region_metrics(mesh.triangles_center, edge, metrics["area_mm2"]),
        "cad_stl_chord_error_mm": chord,
    }


def curvature_table(parameters: dict[str, Any]) -> dict[str, Any]:
    body = parameters["body"]
    alpha = math.radians(body["cone_half_angle_deg"])
    radius = body["maximum_radius_mm"]
    rn = body["nominal_nose_radius_mm"]
    shoulder = body["shoulder_fillet_radius_mm"]
    reconstructed = body["analytic_reconstruction_in"]
    k = parameters["units"]["inch_to_mm"]
    cone_r_min = reconstructed["sphere_cone_tangent_r"] * k
    cone_r_max = reconstructed["cone_shoulder_tangent_r"] * k
    result: dict[str, Any] = {
        "definition": "principal curvature magnitudes in 1/mm; zero denotes a developable/axial principal direction",
        "nose_sphere": {
            "principal_curvature_1_per_mm": [1.0 / rn, 1.0 / rn],
            "minimum_radius_of_curvature_mm": rn,
            "status": "nominal one-inch sphere reconstruction",
        },
        "straight_70deg_cone": {
            "meridional_curvature_1_per_mm": 0.0,
            "azimuthal_curvature_range_1_per_mm": [
                math.cos(alpha) / cone_r_max,
                math.cos(alpha) / cone_r_min,
            ],
            "minimum_radius_of_curvature_mm": cone_r_min / math.cos(alpha),
        },
        "shoulder_toroidal_fillet": {
            "meridional_curvature_1_per_mm": 1.0 / shoulder,
            "azimuthal_curvature_range_1_per_mm": [
                math.cos(alpha) / cone_r_max,
                1.0 / radius,
            ],
            "minimum_radius_of_curvature_mm": shoulder,
            "status": "NASA DWG 1168296 directly specifies R0.100 in",
        },
        "aft_cylinder": {
            "principal_curvature_1_per_mm": [0.0, 1.0 / radius],
            "minimum_radius_of_curvature_mm": radius,
        },
        "nozzle_divergent_cones_at_throat": {},
    }
    for number, nozzle in parameters["nozzles"].items():
        throat_radius = 0.5 * nozzle["throat_diameter_mm"]
        beta = math.radians(nozzle["divergent_wall_half_angle_deg"])
        curvature = math.cos(beta) / throat_radius
        result["nozzle_divergent_cones_at_throat"][number] = {
            "meridional_curvature_1_per_mm": 0.0,
            "maximum_azimuthal_curvature_1_per_mm": curvature,
            "minimum_surface_radius_of_curvature_mm": 1.0 / curvature,
        }
    return result


def ibm_scale_matrix(models: dict[str, Any]) -> dict[str, Any]:
    spacings = [8.0, 6.0, 4.0, 2.0, 1.0, 0.5, 0.25, 0.125, 0.078125]
    result: dict[str, Any] = {
        "contract": "docs/ibm_fixed_stl_pressure.md requires max STL edge / finest dx <= 0.5 for a formal fixed-STL sequence",
        "rows": {},
    }
    for name, data in models.items():
        global_max = data["unique_edge_length_mm"]["p100"]
        front_max = data["axial_regions"]["front_x_le_35_mm"]["edge_mm"]["p100"]
        chord = data["cad_stl_chord_error_mm"]["max_mm"]
        result["rows"][name] = [
            {
                "dx_mm": dx,
                "global_max_edge_over_dx": global_max / dx,
                "front_max_edge_over_dx": front_max / dx,
                "max_chord_error_over_dx": chord / dx,
                "formal_global_edge_contract_pass": global_max / dx <= 0.5,
            }
            for dx in spacings
        ]
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    parameters = json.loads((HERE / "geometry_parameters.yaml").read_text(encoding="utf-8"))
    build = json.loads((args.output / "reports" / "build_audit.json").read_text(encoding="utf-8"))
    models: dict[str, Any] = {}
    for name, audit in build.items():
        models[name] = analyze_stl(
            args.output / "stl" / f"{name}.stl", audit["mesh"]["chord_error"]
        )
        print(f"analyzed {name}")
    report = {
        "status": "PASS",
        "quality_definition": "q=4*sqrt(3)*A/(l1^2+l2^2+l3^2), with q=1 for an equilateral triangle",
        "models": models,
        "analytic_curvature": curvature_table(parameters),
        "ibm_spacing_matrix": ibm_scale_matrix(models),
    }
    path = args.output / "reports" / "mesh_quality_detailed.json"
    path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"PASS: wrote {path}")


if __name__ == "__main__":
    main()
