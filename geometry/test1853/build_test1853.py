#!/usr/bin/env python3
"""Build auditable NASA Test 1853 external CFD geometry with Gmsh/OCC.

The public report does not contain the original CAD-defined forebody/insert
surfaces or the final assembly drawing.  This builder therefore implements the
explicitly documented reconstruction in ``geometry_parameters.yaml`` and
records every derived placement in the output audit.  Coordinates are mm.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import gmsh
import numpy as np
import trimesh


HERE = Path(__file__).resolve().parent
PARAMETER_FILE = HERE / "geometry_parameters.yaml"
OUTPUT_ROOT = HERE / "outputs"


def report_path(path: Path) -> str:
    try:
        return str(path.relative_to(HERE))
    except ValueError:
        return str(path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_parameters(path: Path = PARAMETER_FILE) -> dict[str, Any]:
    # JSON is a strict subset of YAML 1.2; this avoids an optional PyYAML
    # dependency while retaining the requested .yaml audit artifact.
    with path.open("r", encoding="utf-8") as stream:
        return json.load(stream)


def close(a: float, b: float, tol: float, label: str) -> None:
    if abs(a - b) > tol:
        raise RuntimeError(f"{label}: {a:.12g} != {b:.12g} (tol={tol:g})")


def derived_body(p: dict[str, Any]) -> dict[str, float]:
    b = p["body"]
    k = p["units"]["inch_to_mm"]
    alpha = math.radians(b["cone_half_angle_deg"])
    rn = b["nominal_nose_radius_in"] * k
    radius = b["maximum_radius_in"] * k
    shoulder = b["shoulder_fillet_radius_in"] * k
    xc = rn / math.sin(alpha)
    tip = xc - rn
    sphere_tangent_r = rn * math.cos(alpha)
    sphere_tangent_x = rn * math.cos(alpha) ** 2 / math.sin(alpha)
    cone_tangent_r = radius - shoulder * (1.0 - math.cos(alpha))
    cone_tangent_x = cone_tangent_r / math.tan(alpha)
    fillet_center_r = radius - shoulder
    fillet_center_x = radius / math.tan(alpha) + shoulder * math.tan(alpha / 2.0)
    aft_end = b["aft_end_x_in"] * k
    result = {
        "alpha_rad": alpha,
        "nose_radius": rn,
        "sphere_center_x": xc,
        "nose_tip_x": tip,
        "sphere_cone_tangent_x": sphere_tangent_x,
        "sphere_cone_tangent_r": sphere_tangent_r,
        "cone_shoulder_tangent_x": cone_tangent_x,
        "cone_shoulder_tangent_r": cone_tangent_r,
        "shoulder_radius": shoulder,
        "shoulder_center_x": fillet_center_x,
        "shoulder_center_r": fillet_center_r,
        "body_radius": radius,
        "aft_end_x": aft_end,
    }
    published = b["analytic_reconstruction_in"]
    tol = 2.0e-8 * k
    comparisons = {
        "sphere_center_x": "sphere_center_x",
        "nose_tip_x": "nose_tip_x",
        "sphere_cone_tangent_x": "sphere_cone_tangent_x",
        "sphere_cone_tangent_r": "sphere_cone_tangent_r",
        "cone_shoulder_tangent_x": "cone_shoulder_tangent_x",
        "cone_shoulder_tangent_r": "cone_shoulder_tangent_r",
        "shoulder_center_x": "shoulder_circle_center_x",
        "shoulder_center_r": "shoulder_circle_center_r",
    }
    for computed_key, stored_key in comparisons.items():
        close(result[computed_key], published[stored_key] * k, tol, computed_key)
    close(aft_end, b["aft_end_x_mm"], 1.0e-9, "aft_end_x")
    return result


def nozzle_definition(
    p: dict[str, Any],
    body: dict[str, float],
    number: int,
    theta_deg: float | None,
    geometry_key: str | None = None,
    placement_kind: str | None = None,
    orientation_kind: str = "axial",
) -> dict[str, Any]:
    geometry_key = geometry_key or str(number)
    n = p["nozzles"][geometry_key]
    rt = 0.5 * n["throat_diameter_mm"]
    re = 0.5 * n["virtual_exit_diameter_mm"]
    beta = math.radians(n["divergent_wall_half_angle_deg"])
    placement_kind = placement_kind or (
        "center" if theta_deg is None else "peripheral"
    )
    if placement_kind == "center":
        if theta_deg is not None:
            raise ValueError(f"Center nozzle {number} cannot have an azimuth")
        radial_axis = 0.0
        theta = None
        virtual_x = body["sphere_center_x"] - math.sqrt(
            body["nose_radius"] ** 2 - re**2
        )
        placement = "ideal_sphere_flush_reconstruction"
        y_axis = 0.0
        z_axis = 0.0
        if orientation_kind != "axial":
            raise ValueError(
                f"Center nozzle {number} only supports axial orientation"
            )
    elif placement_kind == "peripheral":
        if theta_deg is None:
            raise ValueError(f"Nozzle {number} requires an azimuth")
        radial_axis = p["nozzle_layout"]["peripheral_axis_radius_mm"]
        theta = float(theta_deg)
        angle = math.radians(theta)
        y_axis = radial_axis * math.sin(angle)
        z_axis = radial_axis * math.cos(angle)
        virtual_x = radial_axis / math.tan(body["alpha_rad"])
        placement = "assumed_axis_cone_intersection_virtual_plane"
    else:
        raise ValueError(
            f"Nozzle {number}: placement_kind must be center or peripheral"
        )

    virtual_center = np.array([virtual_x, y_axis, z_axis], dtype=float)
    if orientation_kind == "axial":
        axis_inward = np.array([1.0, 0.0, 0.0])
        # Preserve the established cutter extent so official fixed-STL builds
        # remain byte-for-byte reproducible.
        cutter_start = np.array(
            [
                body["nose_tip_x"]
                - p["mesh"]["cavity_upstream_extension_mm"],
                y_axis,
                z_axis,
            ],
            dtype=float,
        )
        external_length = virtual_x - cutter_start[0]
    elif orientation_kind == "surface_normal":
        if placement_kind != "peripheral":
            raise ValueError("surface_normal orientation requires a peripheral nozzle")
        theta_rad = math.radians(float(theta))
        radial = np.array(
            [0.0, math.sin(theta_rad), math.cos(theta_rad)], dtype=float
        )
        # The body occupies x-r*cot(alpha) >= 0.  This direction points from
        # the tangent virtual-exit plane into that solid; its negative is the
        # upstream/outward jet direction.
        axis_inward = np.array([math.sin(body["alpha_rad"]), 0.0, 0.0])
        axis_inward -= math.cos(body["alpha_rad"]) * radial
        external_length = float(p["mesh"]["cavity_upstream_extension_mm"])
        cutter_start = virtual_center - external_length * axis_inward
        placement = "local_cone_tangent_exit_plane_surface_normal_axis"
    else:
        raise ValueError(
            f"Nozzle {number}: orientation_kind must be axial or surface_normal"
        )

    divergent_length = (re - rt) / math.tan(beta)
    throat_center = virtual_center + divergent_length * axis_inward
    jet_direction = -axis_inward
    r_up = re + external_length * math.tan(beta)
    return {
        "number": number,
        "geometry_key": geometry_key,
        "placement_kind": placement_kind,
        "orientation_kind": orientation_kind,
        "theta_deg": theta,
        "axis_radius_mm": radial_axis,
        "axis_y_mm": y_axis,
        "axis_z_mm": z_axis,
        "virtual_exit_center_mm": virtual_center.tolist(),
        "throat_center_mm": throat_center.tolist(),
        "cutter_start_center_mm": cutter_start.tolist(),
        "axis_inward_unit": axis_inward.tolist(),
        "jet_direction_unit": jet_direction.tolist(),
        "throat_plane_normal_solid_to_fluid": jet_direction.tolist(),
        "throat_radius_mm": rt,
        "virtual_exit_radius_mm": re,
        "divergent_half_angle_deg": n["divergent_wall_half_angle_deg"],
        "virtual_exit_x_mm": virtual_x,
        "virtual_exit_y_mm": float(virtual_center[1]),
        "virtual_exit_z_mm": float(virtual_center[2]),
        "throat_x_mm": float(throat_center[0]),
        "throat_y_mm": float(throat_center[1]),
        "throat_z_mm": float(throat_center[2]),
        "upstream_extension_x_mm": float(cutter_start[0]),
        "upstream_extension_y_mm": float(cutter_start[1]),
        "upstream_extension_z_mm": float(cutter_start[2]),
        "upstream_extension_radius_mm": r_up,
        "divergent_length_mm": divergent_length,
        "expected_throat_area_mm2": n["throat_area_mm2"],
        "placement_status": placement,
    }


def configurations(p: dict[str, Any]) -> dict[str, dict[str, Any]]:
    a = p["nozzle_layout"]["mapping_variants"]
    quad = p["quad_geometry"]
    quad_n3at120: dict[int, float | None] = {1: None}
    quad_n3at120.update({int(k): float(v) for k, v in a["n3_at_120"].items()})
    quad_n3at240: dict[int, float | None] = {1: None}
    quad_n3at240.update({int(k): float(v) for k, v in a["n3_at_240"].items()})
    ring = p["four_peripheral_ring_geometry"]
    ring_active = {
        number: {
            "theta_deg": float(theta),
            "geometry_key": str(ring["nozzle_geometry"]),
            "placement_kind": "peripheral",
        }
        for number, theta in enumerate(ring["azimuths_deg"], start=1)
    }
    surface_normal = p["single_peripheral_surface_normal_geometry"]
    surface_normal_active = {
        1: {
            "theta_deg": float(surface_normal["azimuth_deg"]),
            "geometry_key": str(surface_normal["nozzle_geometry"]),
            "placement_kind": "peripheral",
            "orientation_kind": str(surface_normal["orientation"]),
        }
    }
    return {
        "test1853_run165_single": {
            "description": "Run 165 hardware: active center Nozzle 1; peripheral plugs flush",
            "active": {1: None},
            "run_ids": [165],
            "mapping": "not_applicable",
        },
        "test1853_run262_263_tri_n3at120": {
            "description": "Run 247/262/263 tri hardware; canonical unresolved mapping N3@120, N4@240",
            "active": {int(k): float(v) for k, v in a["n3_at_120"].items()},
            "run_ids": [247, 262, 263],
            "mapping": "N2@0, N3@120, N4@240 (one of two public-data-consistent mappings)",
        },
        "test1853_run262_263_tri_n3at240": {
            "description": "Run 247/262/263 tri hardware; mirrored unresolved mapping N3@240, N4@120",
            "active": {int(k): float(v) for k, v in a["n3_at_240"].items()},
            "run_ids": [247, 262, 263],
            "mapping": "N2@0, N3@240, N4@120 (one of two public-data-consistent mappings)",
        },
        "test1853_quad_n3at120": {
            "description": "Official Test 1853 quad hardware; center N1 plus peripheral N2/N3/N4; N3@120",
            "active": quad_n3at120,
            "run_ids": quad["run_ids"],
            "mapping": "N1@center, N2@0, N3@120, N4@240 (one of two public-data-consistent mappings)",
        },
        "test1853_quad_n3at240": {
            "description": "Official Test 1853 quad hardware; center N1 plus peripheral N2/N3/N4; N3@240",
            "active": quad_n3at240,
            "run_ids": quad["run_ids"],
            "mapping": "N1@center, N2@0, N3@240, N4@120 (one of two public-data-consistent mappings)",
        },
        "test1853_four_peripheral_ring90": {
            "description": "Parametric four-fold ring: four identical peripheral nozzles and a flush center",
            "active": ring_active,
            "run_ids": [],
            "mapping": "Q1@0, Q2@90, Q3@180, Q4@270 on the 63.5-mm pitch circle",
            "classification": "parametric_extension_not_an_official_test1853_run",
        },
        "test1853_single_peripheral_surface_normal": {
            "description": "Parametric single peripheral nozzle with a forebody-tangent exit plane and local surface-normal axis",
            "active": surface_normal_active,
            "run_ids": [],
            "mapping": "Q1 at theta=0 (+z), R=31.75 mm; jet points 20 deg outward from -x",
            "classification": "parametric_extension_not_an_official_test1853_run",
        },
    }


def make_body_occ(body: dict[str, float]) -> int:
    occ = gmsh.model.occ
    point_values = [
        (body["nose_tip_x"], 0.0),
        (body["sphere_cone_tangent_x"], body["sphere_cone_tangent_r"]),
        (body["cone_shoulder_tangent_x"], body["cone_shoulder_tangent_r"]),
        (body["shoulder_center_x"], body["body_radius"]),
        (body["aft_end_x"], body["body_radius"]),
        (body["aft_end_x"], 0.0),
    ]
    points = [occ.addPoint(x, r, 0.0) for x, r in point_values]
    sphere_center = occ.addPoint(body["sphere_center_x"], 0.0, 0.0)
    shoulder_center = occ.addPoint(
        body["shoulder_center_x"], body["shoulder_center_r"], 0.0
    )
    curves = [
        occ.addCircleArc(points[0], sphere_center, points[1]),
        occ.addLine(points[1], points[2]),
        occ.addCircleArc(points[2], shoulder_center, points[3]),
        occ.addLine(points[3], points[4]),
        occ.addLine(points[4], points[5]),
        occ.addLine(points[5], points[0]),
    ]
    wire = occ.addCurveLoop(curves)
    meridional_face = occ.addPlaneSurface([wire])
    revolved = occ.revolve(
        [(2, meridional_face)], 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0 * math.pi
    )
    volumes = [tag for dim, tag in revolved if dim == 3]
    if len(volumes) != 1:
        raise RuntimeError(f"Expected one revolved body volume, got {volumes}")
    return volumes[0]


def build_occ_solid(
    p: dict[str, Any], name: str, cfg: dict[str, Any]
) -> tuple[int, dict[str, float], list[dict[str, Any]]]:
    gmsh.model.add(name)
    body = derived_body(p)
    body_tag = make_body_occ(body)
    cavities: list[tuple[int, int]] = []
    nozzle_data: list[dict[str, Any]] = []
    for number, specification in cfg["active"].items():
        if isinstance(specification, dict):
            theta = specification.get("theta_deg")
            geometry_key = str(specification.get("geometry_key", number))
            placement_kind = str(specification.get("placement_kind", "peripheral"))
            orientation_kind = str(specification.get("orientation_kind", "axial"))
        else:
            theta = specification
            geometry_key = str(number)
            placement_kind = "center" if theta is None else "peripheral"
            orientation_kind = "axial"
        data = nozzle_definition(
            p,
            body,
            number,
            theta,
            geometry_key=geometry_key,
            placement_kind=placement_kind,
            orientation_kind=orientation_kind,
        )
        nozzle_data.append(data)
        cone = gmsh.model.occ.addCone(
            float(data["upstream_extension_x_mm"]),
            float(data["upstream_extension_y_mm"]),
            float(data["upstream_extension_z_mm"]),
            float(data["throat_x_mm"] - data["upstream_extension_x_mm"]),
            float(data["throat_y_mm"] - data["upstream_extension_y_mm"]),
            float(data["throat_z_mm"] - data["upstream_extension_z_mm"]),
            float(data["upstream_extension_radius_mm"]),
            float(data["throat_radius_mm"]),
        )
        cavities.append((3, cone))
    gmsh.model.occ.synchronize()
    cut, _ = gmsh.model.occ.cut(
        [(3, body_tag)], cavities, removeObject=True, removeTool=True
    )
    gmsh.model.occ.synchronize()
    volume_tags = [tag for dim, tag in cut if dim == 3]
    all_volumes = gmsh.model.getEntities(3)
    if len(volume_tags) != 1 or len(all_volumes) != 1:
        raise RuntimeError(
            f"Boolean result must be exactly one solid: cut={cut}, all={all_volumes}"
        )
    volume = volume_tags[0]
    physical = gmsh.model.addPhysicalGroup(3, [volume])
    gmsh.model.setPhysicalName(3, physical, name)
    return volume, body, nozzle_data


def boundary_orientations(volume: int) -> dict[int, int]:
    result: dict[int, int] = {}
    for dim, signed_tag in gmsh.model.getBoundary(
        [(3, volume)], combined=False, oriented=True, recursive=False
    ):
        if dim != 2:
            continue
        result[abs(signed_tag)] = 1 if signed_tag > 0 else -1
    return result


def identify_throat_surfaces(
    volume: int, nozzle_data: list[dict[str, Any]], tolerance_mm: float
) -> dict[int, int]:
    surfaces = boundary_orientations(volume)
    matches: dict[int, int] = {}
    for data in nozzle_data:
        number = int(data["number"])
        target = np.array(
            data["throat_center_mm"],
            dtype=float,
        )
        candidates: list[int] = []
        for surface in surfaces:
            center = np.array(gmsh.model.occ.getCenterOfMass(2, surface), dtype=float)
            area = float(gmsh.model.occ.getMass(2, surface))
            expected_area = float(data["expected_throat_area_mm2"])
            area_error = abs(area - expected_area) / expected_area
            if (
                np.linalg.norm(center - target) <= 5.0 * tolerance_mm
                and area_error <= 1.0e-4
            ):
                candidates.append(surface)
        if len(candidates) != 1:
            raise RuntimeError(
                f"Nozzle {number}: expected one throat cap at {target}, got {candidates}"
            )
        matches[number] = candidates[0]
    if len(set(matches.values())) != len(matches):
        raise RuntimeError(f"Throat surfaces are not unique: {matches}")
    return matches


def export_cad(output: Path, name: str) -> dict[str, int]:
    cad_dir = output / "cad"
    cad_dir.mkdir(parents=True, exist_ok=True)
    paths = {"brep": cad_dir / f"{name}.brep", "step": cad_dir / f"{name}.step"}
    gmsh.write(str(paths["brep"]))
    gmsh.write(str(paths["step"]))
    return {key: path.stat().st_size for key, path in paths.items()}


def set_mesh_options(p: dict[str, Any], quick: bool) -> None:
    mesh = p["mesh"]
    front = float(mesh["front_target_size_mm"])
    aft = float(mesh["aft_target_size_mm"])
    if quick:
        front *= 2.5
        aft *= 2.0
    x0 = float(mesh["transition_start_x_mm"])
    x1 = float(mesh["transition_end_x_mm"])

    def size_callback(dim: int, tag: int, x: float, y: float, z: float, lc: float) -> float:
        del dim, tag, y, z, lc
        if x <= x0:
            return front
        if x >= x1:
            return aft
        fraction = (x - x0) / (x1 - x0)
        return front + fraction * (aft - front)

    gmsh.model.mesh.setSizeCallback(size_callback)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeMin", min(front, aft))
    gmsh.option.setNumber("Mesh.MeshSizeMax", max(front, aft))
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.ElementOrder", int(mesh["element_order"]))
    gmsh.option.setNumber("Mesh.Binary", 1)
    # Gmsh's parallel 2-D mesher can produce equally valid but byte-different
    # triangulations between runs.  A locked facet count and SHA-256 are part
    # of the fixed-STL IBM contract, so formal outputs are generated with one
    # thread for reproducibility.
    gmsh.option.setNumber("General.NumThreads", int(mesh["reproducible_num_threads"]))


def surface_triangles(
    surface: int, orientation: int, lookup: np.ndarray
) -> np.ndarray:
    # Gmsh meshes OCC boundary entities in their topological orientation as
    # embedded in the volume.  The signed value returned by getBoundary is
    # useful topology metadata, but applying it a second time reverses only a
    # subset of faces and breaks winding consistency.
    del orientation
    element_types, _, element_nodes = gmsh.model.mesh.getElements(2, surface)
    blocks: list[np.ndarray] = []
    for element_type, nodes in zip(element_types, element_nodes):
        name, dim, order, count, _, primary = gmsh.model.mesh.getElementProperties(
            element_type
        )
        if dim != 2 or primary != 3:
            raise RuntimeError(
                f"Unsupported surface element {name}: dim={dim}, order={order}, "
                f"nodes={count}, primary={primary}"
            )
        connectivity = np.asarray(nodes, dtype=np.int64).reshape((-1, count))[:, :3]
        faces = lookup[connectivity]
        if np.any(faces < 0):
            raise RuntimeError(f"Surface {surface} references an unknown mesh node")
        blocks.append(faces)
    return np.vstack(blocks) if blocks else np.empty((0, 3), dtype=np.int64)


def collect_mesh(
    volume: int,
) -> tuple[np.ndarray, dict[int, np.ndarray], dict[int, int]]:
    node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
    node_tags = np.asarray(node_tags, dtype=np.int64)
    vertices = np.asarray(coordinates, dtype=float).reshape((-1, 3))
    if len(node_tags) == 0:
        raise RuntimeError("Surface mesh has no nodes")
    lookup = np.full(int(node_tags.max()) + 1, -1, dtype=np.int64)
    lookup[node_tags] = np.arange(len(node_tags), dtype=np.int64)
    orientations = boundary_orientations(volume)
    by_surface = {
        surface: surface_triangles(surface, sign, lookup)
        for surface, sign in orientations.items()
    }
    return vertices, by_surface, orientations


def triangle_area(vertices: np.ndarray, faces: np.ndarray) -> float:
    tri = vertices[faces]
    return float(
        0.5
        * np.linalg.norm(
            np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1
        ).sum()
    )


def write_stl(path: Path, vertices: np.ndarray, faces: np.ndarray) -> dict[str, Any]:
    path.parent.mkdir(parents=True, exist_ok=True)
    mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)
    payload = mesh.export(file_type="stl")
    if not isinstance(payload, bytes):
        payload = payload.encode("utf-8")
    path.write_bytes(payload)
    return {
        "path": report_path(path),
        "bytes": path.stat().st_size,
        "triangles": int(len(faces)),
        "area_mm2": triangle_area(vertices, faces),
    }


def chord_error_by_surface(
    vertices: np.ndarray, by_surface: dict[int, np.ndarray]
) -> dict[str, float]:
    maximum = 0.0
    sum_sq = 0.0
    count = 0
    for surface, faces in by_surface.items():
        if len(faces) == 0:
            continue
        tri = vertices[faces]
        samples = np.concatenate(
            [
                0.5 * (tri[:, 0] + tri[:, 1]),
                0.5 * (tri[:, 1] + tri[:, 2]),
                0.5 * (tri[:, 2] + tri[:, 0]),
                tri.mean(axis=1),
            ],
            axis=0,
        )
        for start in range(0, len(samples), 50000):
            chunk = samples[start : start + 50000]
            closest, _ = gmsh.model.getClosestPoint(
                2, surface, chunk.reshape(-1).tolist()
            )
            closest_array = np.asarray(closest, dtype=float).reshape((-1, 3))
            distances = np.linalg.norm(chunk - closest_array, axis=1)
            maximum = max(maximum, float(distances.max(initial=0.0)))
            sum_sq += float(np.dot(distances, distances))
            count += len(distances)
    return {
        "sample_count": count,
        "max_mm": maximum,
        "rms_mm": math.sqrt(sum_sq / count) if count else 0.0,
        "sample_scheme": "all triangle edge midpoints and centroids projected to source OCC surface",
    }


def export_surface_meshes(
    p: dict[str, Any],
    output: Path,
    name: str,
    volume: int,
    throat_surfaces: dict[int, int],
    nozzle_data: list[dict[str, Any]],
    quick: bool,
) -> tuple[dict[str, Any], np.ndarray, dict[int, np.ndarray]]:
    vertices, by_surface, _ = collect_mesh(volume)
    all_faces = np.vstack(list(by_surface.values()))
    report: dict[str, Any] = {}
    report["closed_solid"] = write_stl(
        output / "stl" / f"{name}.stl", vertices, all_faces
    )
    cap_tags = set(throat_surfaces.values())
    wall_faces = np.vstack(
        [faces for surface, faces in by_surface.items() if surface not in cap_tags]
    )
    patch_dir = output / "patches" / name
    report["body_wall_open"] = write_stl(
        patch_dir / "body_wall_excluding_throat_caps.stl", vertices, wall_faces
    )
    report["throat_caps"] = {}
    for number, surface in sorted(throat_surfaces.items()):
        faces = by_surface[surface]
        item = write_stl(
            patch_dir / f"nozzle_{number}_throat_cap.stl", vertices, faces
        )
        item["occ_surface_tag"] = surface
        item["cad_area_mm2"] = float(gmsh.model.occ.getMass(2, surface))
        nozzle = next(item for item in nozzle_data if int(item["number"]) == number)
        expected = float(nozzle["expected_throat_area_mm2"])
        item["expected_as_built_area_mm2"] = expected
        item["cad_area_relative_error"] = abs(item["cad_area_mm2"] - expected) / expected
        item["mesh_area_relative_error"] = abs(item["area_mm2"] - expected) / expected
        if item["cad_area_relative_error"] > p["mesh"]["cad_area_relative_tolerance"]:
            raise RuntimeError(f"Nozzle {number} OCC throat area mismatch: {item}")
        mesh_area_tolerance = float(p["mesh"]["mesh_area_relative_tolerance"])
        if quick:
            mesh_area_tolerance = max(mesh_area_tolerance, 0.025)
        if item["mesh_area_relative_error"] > mesh_area_tolerance:
            raise RuntimeError(f"Nozzle {number} mesh throat area mismatch: {item}")
        report["throat_caps"][str(number)] = item
    report["chord_error"] = chord_error_by_surface(vertices, by_surface)
    return report, vertices, by_surface


def make_reference_disk(
    path: Path, x: float, y0: float, z0: float, radius: float, segments: int
) -> None:
    angles = np.linspace(0.0, 2.0 * math.pi, segments, endpoint=False)
    ring = np.column_stack(
        [
            np.full(segments, x),
            y0 + radius * np.cos(angles),
            z0 + radius * np.sin(angles),
        ]
    )
    vertices = np.vstack([[x, y0, z0], ring])
    # Viewed along +x the ring is counter-clockwise; this order gives -x,
    # the upstream/exhaust normal of the reference plane.
    faces = np.array(
        [[0, 1 + ((i + 1) % segments), 1 + i] for i in range(segments)],
        dtype=np.int64,
    )
    write_stl(path, vertices, faces)


def oriented_disk_frame(normal: list[float]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    unit_normal = np.asarray(normal, dtype=float)
    unit_normal /= np.linalg.norm(unit_normal)
    reference = (
        np.array([0.0, 0.0, 1.0])
        if abs(unit_normal[2]) < 0.9
        else np.array([0.0, 1.0, 0.0])
    )
    tangent_1 = np.cross(reference, unit_normal)
    tangent_1 /= np.linalg.norm(tangent_1)
    tangent_2 = np.cross(unit_normal, tangent_1)
    return unit_normal, tangent_1, tangent_2


def make_oriented_reference_disk(
    path: Path,
    center: list[float],
    normal: list[float],
    radius: float,
    segments: int,
) -> None:
    unit_normal, tangent_1, tangent_2 = oriented_disk_frame(normal)
    del unit_normal
    center_array = np.asarray(center, dtype=float)
    angles = np.linspace(0.0, 2.0 * math.pi, segments, endpoint=False)
    ring = center_array + radius * (
        np.cos(angles)[:, None] * tangent_1
        + np.sin(angles)[:, None] * tangent_2
    )
    vertices = np.vstack([center_array, ring])
    faces = np.array(
        [[0, 1 + i, 1 + ((i + 1) % segments)] for i in range(segments)],
        dtype=np.int64,
    )
    write_stl(path, vertices, faces)


def peripheral_lip_rho(
    rho_virtual: float, axis_radius: float, alpha: float, beta: float, psi: float
) -> float:
    q = 1.0 / math.tan(alpha)
    m = 1.0 / math.tan(beta)

    def residual(rho: float) -> float:
        global_r = math.sqrt(
            axis_radius**2
            + rho**2
            + 2.0 * axis_radius * rho * math.cos(psi)
        )
        return m * (rho_virtual - rho) - q * (global_r - axis_radius)

    lo, hi = 0.0, 2.0 * rho_virtual
    while residual(hi) > 0.0:
        hi *= 2.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if residual(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def export_references(
    p: dict[str, Any], output: Path, name: str, nozzle_data: list[dict[str, Any]]
) -> dict[str, Any]:
    ref_dir = output / "reference" / name
    ref_dir.mkdir(parents=True, exist_ok=True)
    segments = int(p["mesh"]["reference_disk_segments"])
    samples = int(p["mesh"]["lip_samples"])
    body = derived_body(p)
    rows: list[dict[str, Any]] = []
    disk_files: list[str] = []
    for data in nozzle_data:
        number = int(data["number"])
        disk = ref_dir / f"nozzle_{number}_virtual_exit_REFERENCE_ONLY.stl"
        if data["orientation_kind"] == "axial":
            make_reference_disk(
                disk,
                float(data["virtual_exit_x_mm"]),
                float(data["axis_y_mm"]),
                float(data["axis_z_mm"]),
                float(data["virtual_exit_radius_mm"]),
                segments,
            )
        else:
            make_oriented_reference_disk(
                disk,
                list(data["virtual_exit_center_mm"]),
                list(data["jet_direction_unit"]),
                float(data["virtual_exit_radius_mm"]),
                segments,
            )
        disk_files.append(report_path(disk))
        if data["orientation_kind"] != "axial":
            _, tangent_1, tangent_2 = oriented_disk_frame(
                list(data["jet_direction_unit"])
            )
            center = np.asarray(data["virtual_exit_center_mm"], dtype=float)
            radius = float(data["virtual_exit_radius_mm"])
            for index in range(samples):
                psi = 2.0 * math.pi * index / samples
                point = center + radius * (
                    math.cos(psi) * tangent_1 + math.sin(psi) * tangent_2
                )
                rows.append(
                    {
                        "nozzle": number,
                        "sample": index,
                        "psi_deg": math.degrees(psi),
                        "x_mm": float(point[0]),
                        "y_mm": float(point[1]),
                        "z_mm": float(point[2]),
                        "local_rho_mm": radius,
                        "status": "tangent_virtual_exit_reference; physical lip is the Boolean passage/OML intersection",
                    }
                )
            continue
        theta = math.radians(float(data["theta_deg"] or 0.0))
        radial_yz = np.array([math.sin(theta), math.cos(theta)])
        tangent_yz = np.array([math.cos(theta), -math.sin(theta)])
        for index in range(samples):
            psi = 2.0 * math.pi * index / samples
            if data["placement_kind"] == "center":
                rho = float(data["virtual_exit_radius_mm"])
                x = float(data["virtual_exit_x_mm"])
            else:
                rho = peripheral_lip_rho(
                    float(data["virtual_exit_radius_mm"]),
                    float(data["axis_radius_mm"]),
                    body["alpha_rad"],
                    math.radians(float(data["divergent_half_angle_deg"])),
                    psi,
                )
                x = float(data["virtual_exit_x_mm"]) + (
                    float(data["virtual_exit_radius_mm"]) - rho
                ) / math.tan(math.radians(float(data["divergent_half_angle_deg"])))
            offset = rho * (
                math.cos(psi) * radial_yz + math.sin(psi) * tangent_yz
            )
            rows.append(
                {
                    "nozzle": number,
                    "sample": index,
                    "psi_deg": math.degrees(psi),
                    "x_mm": x,
                    "y_mm": float(data["axis_y_mm"]) + offset[0],
                    "z_mm": float(data["axis_z_mm"]) + offset[1],
                    "local_rho_mm": rho,
                    "status": data["placement_status"],
                }
            )
    csv_path = ref_dir / "physical_lip_reconstruction.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return {
        "virtual_exit_disks": disk_files,
        "lip_csv": report_path(csv_path),
        "lip_rows": len(rows),
        "warning": "virtual-exit disks are reference planes, not peripheral physical openings; surface-normal physical lips remain Boolean passage/OML intersections",
    }


def validate_closed_stl(path: Path) -> dict[str, Any]:
    mesh = trimesh.load_mesh(path, process=True)
    if not isinstance(mesh, trimesh.Trimesh):
        raise RuntimeError(f"{path}: expected one Trimesh, got {type(mesh)}")
    checks = {
        "vertices": int(len(mesh.vertices)),
        "triangles": int(len(mesh.faces)),
        "watertight": bool(mesh.is_watertight),
        "winding_consistent": bool(mesh.is_winding_consistent),
        "is_volume": bool(mesh.is_volume),
        "body_count": int(mesh.body_count),
        "volume_mm3": float(mesh.volume),
        "bounds_mm": mesh.bounds.tolist(),
    }
    if not (
        checks["watertight"]
        and checks["winding_consistent"]
        and checks["is_volume"]
        and checks["body_count"] == 1
        and checks["volume_mm3"] > 0.0
    ):
        raise RuntimeError(f"Closed-STL validation failed for {path}: {checks}")
    return checks


def validate_reimport(path: Path, expected_bounds: list[list[float]]) -> dict[str, Any]:
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("Geometry.OCCBoundsUseStl", 1)
        gmsh.model.add("reimport_check")
        gmsh.model.occ.importShapes(str(path), highestDimOnly=False)
        gmsh.model.occ.synchronize()
        volumes = gmsh.model.getEntities(3)
        if len(volumes) != 1:
            raise RuntimeError(f"{path}: re-import produced {len(volumes)} solids")
        bounds = np.asarray(gmsh.model.getBoundingBox(3, volumes[0][1])).reshape((2, 3))
        expected = np.asarray(expected_bounds)
        maximum_difference = float(np.max(np.abs(bounds - expected)))
        if maximum_difference > 2.0e-4:
            raise RuntimeError(
                f"{path}: re-import bounds changed by {maximum_difference:g} mm"
            )
        return {
            "solid_count": 1,
            "bounds_mm": bounds.tolist(),
            "max_bounds_difference_mm": maximum_difference,
        }
    finally:
        gmsh.finalize()


def write_manifest(output: Path, audits: dict[str, Any]) -> None:
    reports = output / "reports"
    reports.mkdir(parents=True, exist_ok=True)
    manifest = reports / "build_audit.json"
    manifest.write_text(
        json.dumps(audits, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    lines = [
        "# Test 1853 geometry build audit",
        "",
        "All CAD coordinates are millimetres. `PASS` means that the generated artifact met the checks below; it does not turn a documented public-data reconstruction into the unavailable NASA original CAD.",
        "",
        "| Model | OCC solids | STL triangles | Watertight | Volume | Max CAD–STL chord sample |",
        "|---|---:|---:|---|---|---:|",
    ]
    for name, audit in audits.items():
        validation = audit["stl_validation"]
        lines.append(
            f"| `{name}` | {audit['occ']['solid_count']} | "
            f"{validation['triangles']} | {validation['watertight']} | "
            f"{validation['is_volume']} | {audit['mesh']['chord_error']['max_mm']:.6g} mm |"
        )
    lines += [
        "",
        "## Interpretation boundary",
        "",
        "- Run 165 is the center single-nozzle configuration. Runs 247/262/263 use the same tri-nozzle hardware with run-specific flow and roll conditions.",
        "- Runs 307/315 use the quad hardware: center Nozzle 1 plus peripheral Nozzles 2/3/4. The phi=0 CAD/STL masters remain run-neutral; the preparation tool applies the selected run roll.",
        "- `test1853_four_peripheral_ring90` is a parametric four-fold extension with no center nozzle and no official Test 1853 run ID; it must not be presented as Run 307/315 hardware.",
        "- `test1853_single_peripheral_surface_normal` is a parametric one-nozzle extension. Its exit plane is tangent to the local 70-degree cone and its axis is surface-normal; it is not official Test 1853 hardware.",
        "- The body OML is the documented nominal tangent reconstruction. The public drawings explicitly require missing CAD surfaces.",
        "- Peripheral virtual-exit axial placement is an explicit CFD reconstruction assumption; the physical lip in these solids is the Boolean cone/OML intersection.",
        "- Both Nozzle-3/Nozzle-4 azimuth assignments are exported because the public report does not identify their one-to-one mapping.",
        "- The throat caps are CFD boundary patches. Complete convergent passages, internal plenum, installed sting services, and pressure-port holes are not guessed.",
        "",
        "Detailed numeric results are in `build_audit.json`.",
    ]
    (reports / "geometry_audit.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_run_transforms(p: dict[str, Any], output: Path) -> None:
    reference = output / "reference"
    reference.mkdir(parents=True, exist_ok=True)
    result: dict[str, Any] = {
        "convention": "homogeneous rotation about +x by reported model roll phi; CAD master is phi=0",
        "units": "mm and degrees",
        "runs": {},
    }
    for run, data in p["runs"].items():
        phi = math.radians(float(data["roll_phi_deg"]))
        c, s = math.cos(phi), math.sin(phi)
        result["runs"][run] = {
            "configuration": data["configuration"],
            "roll_phi_deg": data["roll_phi_deg"],
            "model_to_tunnel_homogeneous": [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, c, -s, 0.0],
                [0.0, s, c, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        }
    (reference / "run_roll_transforms.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def build_one(
    p: dict[str, Any], output: Path, name: str, cfg: dict[str, Any], quick: bool
) -> dict[str, Any]:
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 1)
        gmsh.option.setString("Geometry.OCCTargetUnit", "MM")
        # OpenCASCADE's default box is the box of the untrimmed support
        # surface; for the shoulder torus that is intentionally loose.
        gmsh.option.setNumber("Geometry.OCCBoundsUseStl", 1)
        volume, body, nozzle_data = build_occ_solid(p, name, cfg)
        tolerance = float(p["mesh"]["key_dimension_tolerance_mm"])
        throats = identify_throat_surfaces(volume, nozzle_data, tolerance)
        bbox = np.asarray(gmsh.model.getBoundingBox(3, volume), dtype=float).reshape((2, 3))
        expected_front_x = body["nose_tip_x"]
        for nozzle in nozzle_data:
            if nozzle["placement_kind"] == "center":
                # The center nozzle replaces the spherical nose cap.  The
                # most-upstream remaining solid is its flush circular lip.
                expected_front_x = float(nozzle["virtual_exit_x_mm"])
        expected_bounds = np.array(
            [
                [expected_front_x, -body["body_radius"], -body["body_radius"]],
                [body["aft_end_x"], body["body_radius"], body["body_radius"]],
            ]
        )
        max_bounds_error = float(np.max(np.abs(bbox - expected_bounds)))
        if max_bounds_error > tolerance:
            raise RuntimeError(
                f"{name}: OCC bounds differ from expected by {max_bounds_error:g} mm"
            )
        pre_export_volume = float(gmsh.model.occ.getMass(3, volume))
        cad_sizes = export_cad(output, name)
        # Re-import the BREP master before tessellation.  Besides proving that
        # the master is reusable, this normalizes OpenCASCADE face
        # orientations after BooleanDifference; meshing the transient Boolean
        # result directly can carry mixed per-face orientations in Gmsh.
        gmsh.clear()
        gmsh.model.add(f"{name}_mesh_from_brep")
        gmsh.model.occ.importShapes(str(output / "cad" / f"{name}.brep"))
        gmsh.model.occ.synchronize()
        imported_volumes = gmsh.model.getEntities(3)
        if len(imported_volumes) != 1:
            raise RuntimeError(
                f"{name}: BREP master re-import gave {len(imported_volumes)} solids"
            )
        volume = imported_volumes[0][1]
        throats = identify_throat_surfaces(volume, nozzle_data, tolerance)
        imported_volume = float(gmsh.model.occ.getMass(3, volume))
        close(
            imported_volume,
            pre_export_volume,
            max(1.0e-5, abs(pre_export_volume) * 1.0e-10),
            f"{name} BREP re-import volume",
        )
        bbox = np.asarray(gmsh.model.getBoundingBox(3, volume), dtype=float).reshape((2, 3))
        set_mesh_options(p, quick)
        gmsh.model.mesh.generate(2)
        mesh_report, _, _ = export_surface_meshes(
            p, output, name, volume, throats, nozzle_data, quick
        )
        target_chord = float(p["mesh"]["target_max_cad_chord_error_mm"])
        if not quick and mesh_report["chord_error"]["max_mm"] > target_chord:
            raise RuntimeError(
                f"{name}: sampled chord error {mesh_report['chord_error']['max_mm']:.6g} "
                f"mm exceeds {target_chord:g} mm"
            )
        reference_report = export_references(p, output, name, nozzle_data)
        occ_report = {
            "solid_count": len(gmsh.model.getEntities(3)),
            "surface_count": len(boundary_orientations(volume)),
            "volume_mm3": float(gmsh.model.occ.getMass(3, volume)),
            "bounds_mm": bbox.tolist(),
            "expected_bounds_mm": expected_bounds.tolist(),
            "max_bounds_error_mm": max_bounds_error,
            "throat_surface_tags": {str(k): v for k, v in throats.items()},
        }
    finally:
        gmsh.finalize()

    stl_path = output / "stl" / f"{name}.stl"
    stl_validation = validate_closed_stl(stl_path)
    step_path = output / "cad" / f"{name}.step"
    brep_path = output / "cad" / f"{name}.brep"
    step_reimport = validate_reimport(step_path, occ_report["bounds_mm"])
    brep_reimport = validate_reimport(brep_path, occ_report["bounds_mm"])
    return {
        "status": "PASS",
        "description": cfg["description"],
        "run_ids": cfg["run_ids"],
        "mapping": cfg["mapping"],
        "quick_mesh": quick,
        "nozzles": nozzle_data,
        "occ": occ_report,
        "cad_file_bytes": cad_sizes,
        "mesh": mesh_report,
        "reference_outputs": reference_report,
        "stl_validation": stl_validation,
        "step_reimport": step_reimport,
        "brep_reimport": brep_reimport,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=OUTPUT_ROOT, help="output directory"
    )
    parser.add_argument(
        "--config",
        action="append",
        help="model stem to build; repeat as needed (default: all six)",
    )
    parser.add_argument(
        "--quick", action="store_true", help="coarser development mesh; relax chord gate"
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    p = load_parameters()
    all_configs = configurations(p)
    selected = args.config or list(all_configs)
    unknown = sorted(set(selected) - set(all_configs))
    if unknown:
        raise SystemExit(f"Unknown configurations: {unknown}; choose from {list(all_configs)}")
    args.output.mkdir(parents=True, exist_ok=True)
    audits: dict[str, Any] = {}
    for name in selected:
        print(f"\n=== Building {name} ===", flush=True)
        audits[name] = build_one(p, args.output, name, all_configs[name], args.quick)
        audits[name]["parameter_file_sha256"] = sha256(PARAMETER_FILE)
        audits[name]["builder_sha256"] = sha256(Path(__file__).resolve())
    write_manifest(args.output, audits)
    write_run_transforms(p, args.output)
    print(f"\nPASS: generated {len(audits)} Test 1853 model(s) in {args.output}")


if __name__ == "__main__":
    main()
