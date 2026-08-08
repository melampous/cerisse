#!/usr/bin/env python3
"""Audit fixed IBM curve/STL geometry before a production grid sequence.

The audit is intentionally independent of the flow solver.  It checks the
geometric contract needed by a static ghost-cell IBM: finite coordinates,
non-degenerate elements, closed/manifold topology, consistent orientation,
surface-normal closure, and element size relative to requested fluid spacings.

It does not prove absence of triangle/triangle self-intersections.  That test
requires a robust geometry kernel and remains a separate CAD/CGAL check.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import struct
from collections import defaultdict
from pathlib import Path

import numpy as np


class UnionFind:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))
        self.rank = [0] * size

    def find(self, item: int) -> int:
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, first: int, second: int) -> None:
        a, b = self.find(first), self.find(second)
        if a == b:
            return
        if self.rank[a] < self.rank[b]:
            a, b = b, a
        self.parent[b] = a
        if self.rank[a] == self.rank[b]:
            self.rank[a] += 1


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_stl(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Return triangle vertices and stored facet normals."""
    data = path.read_bytes()
    if len(data) >= 84:
        count = struct.unpack_from("<I", data, 80)[0]
        if len(data) == 84 + 50 * count:
            raw = np.empty((count, 12), dtype=float)
            for index in range(count):
                raw[index] = struct.unpack_from(
                    "<12fH", data, 84 + 50 * index
                )[:12]
            return raw[:, 3:].reshape(count, 3, 3), raw[:, :3]

    text = data.decode("ascii", errors="strict")
    normals: list[list[float]] = []
    vertices: list[list[float]] = []
    for lineno, line in enumerate(text.splitlines(), start=1):
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "facet" and len(fields) >= 5 and fields[1] == "normal":
            normals.append([float(value) for value in fields[2:5]])
        elif fields[0] == "vertex" and len(fields) >= 4:
            vertices.append([float(value) for value in fields[1:4]])
    if not vertices or len(vertices) % 3:
        raise ValueError(f"{path}: malformed ASCII STL vertex count")
    triangles = np.asarray(vertices, dtype=float).reshape(-1, 3, 3)
    stored = np.asarray(normals, dtype=float)
    if stored.shape != (triangles.shape[0], 3):
        raise ValueError(f"{path}: facet-normal/triangle count mismatch")
    return triangles, stored


def read_curve(path: Path) -> np.ndarray:
    points: list[list[float]] = []
    for lineno, line in enumerate(path.read_text().splitlines(), start=1):
        fields = line.split()
        if not fields or fields[0].startswith("#"):
            continue
        if len(fields) < 2:
            raise ValueError(f"{path}:{lineno}: expected x y")
        points.append([float(fields[0]), float(fields[1])])
    array = np.asarray(points, dtype=float)
    if array.shape[0] < 3:
        raise ValueError(f"{path}: a closed IBM curve needs at least 3 points")
    scale = max(float(np.ptp(array, axis=0).max()), 1.0)
    if np.linalg.norm(array[0] - array[-1]) <= 1.0e-13 * scale:
        array = array[:-1]
    return array


def quantized_vertex_ids(vertices: np.ndarray, tolerance: float) -> tuple[np.ndarray, int]:
    # Match src/ibm/ibm_backend_bvh.h: round absolute coordinates at 1e-7
    # unless the caller deliberately supplies another tolerance.
    scaled = vertices / tolerance
    # std::round, used by the solver, rounds halfway cases away from zero;
    # np.rint uses ties-to-even and can therefore produce a different weld.
    keys = (
        np.sign(scaled) * np.floor(np.abs(scaled) + 0.5)
    ).astype(np.int64)
    mapping: dict[tuple[int, ...], int] = {}
    ids = np.empty(keys.shape[0], dtype=np.int64)
    for index, key_array in enumerate(keys):
        key = tuple(int(value) for value in key_array)
        if key not in mapping:
            mapping[key] = len(mapping)
        ids[index] = mapping[key]
    return ids, len(mapping)


def stl_audit(path: Path, args: argparse.Namespace) -> dict[str, object]:
    triangles, stored_normals = read_stl(path)
    finite = bool(np.isfinite(triangles).all() and np.isfinite(stored_normals).all())
    if not finite:
        raise ValueError(f"{path}: nonfinite STL coordinates or normals")

    bounds_lo = np.min(triangles.reshape(-1, 3), axis=0)
    bounds_hi = np.max(triangles.reshape(-1, 3), axis=0)
    extent = bounds_hi - bounds_lo
    scale = max(float(np.max(extent)), np.finfo(float).tiny)
    weld_tolerance = args.weld_tolerance

    edge_vectors = np.stack(
        (
            triangles[:, 1] - triangles[:, 0],
            triangles[:, 2] - triangles[:, 1],
            triangles[:, 0] - triangles[:, 2],
        ),
        axis=1,
    )
    edge_lengths = np.linalg.norm(edge_vectors, axis=2)
    area_vectors = 0.5 * np.cross(
        triangles[:, 1] - triangles[:, 0],
        triangles[:, 2] - triangles[:, 0],
    )
    areas = np.linalg.norm(area_vectors, axis=1)
    area_tolerance = args.degenerate_tolerance * scale * scale
    degenerate = areas <= area_tolerance
    total_area = float(np.sum(areas))
    closure = np.sum(area_vectors, axis=0)
    closure_relative = float(np.linalg.norm(closure) / max(total_area, np.finfo(float).tiny))

    flat_vertices = triangles.reshape(-1, 3)
    flat_ids, unique_vertices = quantized_vertex_ids(flat_vertices, weld_tolerance)
    face_ids = flat_ids.reshape(-1, 3)
    weld_collapsed = np.array(
        [len(set(int(value) for value in ids)) < 3 for ids in face_ids],
        dtype=bool,
    )
    edge_faces: dict[tuple[int, int], list[tuple[int, int]]] = defaultdict(list)
    canonical_faces: dict[tuple[int, int, int], int] = defaultdict(int)
    for face, ids in enumerate(face_ids):
        canonical_faces[tuple(sorted(int(value) for value in ids))] += 1
        for start, end in ((ids[0], ids[1]), (ids[1], ids[2]), (ids[2], ids[0])):
            a, b = int(start), int(end)
            edge_faces[(min(a, b), max(a, b))].append((face, 1 if a < b else -1))

    boundary_edges = sum(len(owners) == 1 for owners in edge_faces.values())
    nonmanifold_edges = sum(len(owners) > 2 for owners in edge_faces.values())
    inconsistent_edges = sum(
        len(owners) == 2 and owners[0][1] == owners[1][1]
        for owners in edge_faces.values()
    )
    duplicate_faces = sum(count - 1 for count in canonical_faces.values() if count > 1)

    components = UnionFind(triangles.shape[0])
    for owners in edge_faces.values():
        first = owners[0][0]
        for face, _ in owners[1:]:
            components.union(first, face)
    component_count = len({components.find(face) for face in range(triangles.shape[0])})

    geometric_normals = np.zeros_like(area_vectors)
    valid_area = ~degenerate
    geometric_normals[valid_area] = (
        area_vectors[valid_area] / areas[valid_area, None]
    )
    stored_magnitude = np.linalg.norm(stored_normals, axis=1)
    valid_stored = stored_magnitude > np.finfo(float).tiny
    stored_unit = np.zeros_like(stored_normals)
    stored_unit[valid_stored] = stored_normals[valid_stored] / stored_magnitude[valid_stored, None]
    normal_dot = np.einsum("ij,ij->i", stored_unit, geometric_normals)
    normal_mismatch = valid_area & valid_stored & (normal_dot < 1.0 - 1.0e-6)
    reversed_normals = valid_area & valid_stored & (normal_dot < 0.0)

    signed_volume = float(
        np.sum(
            np.einsum(
                "ij,ij->i",
                triangles[:, 0],
                np.cross(triangles[:, 1], triangles[:, 2]),
            )
        )
        / 6.0
    )

    failures: list[str] = []
    if np.any(degenerate):
        failures.append("degenerate facets")
    if boundary_edges:
        failures.append("open boundary edges")
    if nonmanifold_edges:
        failures.append("non-manifold edges")
    if inconsistent_edges:
        failures.append("inconsistent neighboring face orientation")
    if duplicate_faces:
        failures.append("duplicate facets")
    if np.any(weld_collapsed):
        failures.append("facets collapse under the solver STL weld tolerance")
    if closure_relative > args.normal_closure_tolerance:
        failures.append("surface-normal closure tolerance exceeded")

    ratios = []
    for spacing in args.dx:
        ratio = float(np.max(edge_lengths) / spacing)
        ratios.append({"dx": spacing, "max_edge_over_dx": ratio})
        if args.max_edge_over_dx is not None and ratio > args.max_edge_over_dx:
            failures.append(f"max edge/dx exceeds limit at dx={spacing:g}")

    warnings: list[str] = []
    if signed_volume <= 0.0:
        warnings.append(
            "input faces are inward-oriented; the closed-mesh BVH reader will "
            "reverse them before building IBM normals"
        )
        if args.require_outward_input and not args.allow_inward:
            failures.append("input STL is not outward-oriented")
    if np.any(~valid_stored) or np.any(normal_mismatch):
        warnings.append(
            "stored facet normals are invalid/inconsistent; the BVH reader "
            "currently recomputes normals from oriented vertices"
        )

    return {
        "kind": "stl",
        "sha256": file_sha256(path),
        "triangles": int(triangles.shape[0]),
        "unique_vertices": unique_vertices,
        "bounds_lo": bounds_lo.tolist(),
        "bounds_hi": bounds_hi.tolist(),
        "weld_tolerance": weld_tolerance,
        "edge_length": {
            "min": float(np.min(edge_lengths)),
            "mean": float(np.mean(edge_lengths)),
            "max": float(np.max(edge_lengths)),
        },
        "area": total_area,
        "signed_volume": signed_volume,
        "normal_closure": closure.tolist(),
        "normal_closure_relative": closure_relative,
        "degenerate_facets": int(np.count_nonzero(degenerate)),
        "weld_collapsed_facets": int(np.count_nonzero(weld_collapsed)),
        "boundary_edges": int(boundary_edges),
        "nonmanifold_edges": int(nonmanifold_edges),
        "inconsistent_edges": int(inconsistent_edges),
        "duplicate_facets": int(duplicate_faces),
        "connected_components": component_count,
        "zero_stored_normals": int(np.count_nonzero(~valid_stored)),
        "mismatched_stored_normals": int(np.count_nonzero(normal_mismatch)),
        "reversed_stored_normals": int(np.count_nonzero(reversed_normals)),
        "fluid_spacing_ratios": ratios,
        "self_intersection_checked": False,
        "warnings": warnings,
        "failures": failures,
        "passed": not failures,
    }


def orientation(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    ab = b - a
    ac = c - a
    return float(ab[0] * ac[1] - ab[1] * ac[0])


def proper_segment_intersection(
    a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray, tolerance: float
) -> bool:
    o1, o2 = orientation(a, b, c), orientation(a, b, d)
    o3, o4 = orientation(c, d, a), orientation(c, d, b)
    return o1 * o2 < -tolerance * tolerance and o3 * o4 < -tolerance * tolerance


def curve_audit(path: Path, args: argparse.Namespace) -> dict[str, object]:
    points = read_curve(path)
    if not np.isfinite(points).all():
        raise ValueError(f"{path}: nonfinite curve coordinates")
    following = np.roll(points, -1, axis=0)
    segments = following - points
    lengths = np.linalg.norm(segments, axis=1)
    bounds_lo, bounds_hi = np.min(points, axis=0), np.max(points, axis=0)
    scale = max(float(np.max(bounds_hi - bounds_lo)), np.finfo(float).tiny)
    length_tolerance = args.degenerate_tolerance * scale
    degenerate = lengths <= length_tolerance
    signed_area = 0.5 * float(
        np.sum(points[:, 0] * following[:, 1] - following[:, 0] * points[:, 1])
    )
    perimeter = float(np.sum(lengths))
    normal_vectors = np.column_stack((segments[:, 1], -segments[:, 0]))
    closure = np.sum(normal_vectors, axis=0)
    closure_relative = float(np.linalg.norm(closure) / max(perimeter, np.finfo(float).tiny))

    intersections = 0
    if args.check_self_intersections:
        count = points.shape[0]
        tolerance = args.degenerate_tolerance * scale
        for first in range(count):
            for second in range(first + 1, count):
                if second in (first, (first + 1) % count) or first == (second + 1) % count:
                    continue
                if proper_segment_intersection(
                    points[first], following[first],
                    points[second], following[second], tolerance
                ):
                    intersections += 1

    failures: list[str] = []
    if np.any(degenerate):
        failures.append("degenerate curve segments")
    if closure_relative > args.normal_closure_tolerance:
        failures.append("curve-normal closure tolerance exceeded")
    if intersections:
        failures.append("self-intersecting curve")

    ratios = []
    for spacing in args.dx:
        ratio = float(np.max(lengths) / spacing)
        ratios.append({"dx": spacing, "max_edge_over_dx": ratio})
        if args.max_edge_over_dx is not None and ratio > args.max_edge_over_dx:
            failures.append(f"max segment/dx exceeds limit at dx={spacing:g}")

    warnings: list[str] = []
    if signed_area <= 0.0:
        warnings.append(
            "input curve is clockwise; the BVH reader will reverse it to CCW "
            "before building outward normals"
        )
        if args.require_outward_input and not args.allow_inward:
            failures.append("input curve is not counter-clockwise")

    return {
        "kind": "curve",
        "sha256": file_sha256(path),
        "points": int(points.shape[0]),
        "bounds_lo": bounds_lo.tolist(),
        "bounds_hi": bounds_hi.tolist(),
        "segment_length": {
            "min": float(np.min(lengths)),
            "mean": float(np.mean(lengths)),
            "max": float(np.max(lengths)),
        },
        "perimeter": perimeter,
        "signed_area": signed_area,
        "normal_closure": closure.tolist(),
        "normal_closure_relative": closure_relative,
        "degenerate_segments": int(np.count_nonzero(degenerate)),
        "self_intersection_checked": bool(args.check_self_intersections),
        "proper_self_intersections": intersections,
        "fluid_spacing_ratios": ratios,
        "warnings": warnings,
        "failures": failures,
        "passed": not failures,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("geometry", type=Path)
    parser.add_argument("--kind", choices=("auto", "stl", "curve"), default="auto")
    parser.add_argument("--dx", type=float, nargs="*", default=[],
                        help="Finest or sequence of physical fluid-cell spacings.")
    parser.add_argument("--max-edge-over-dx", type=float,
                        help="Fail if the maximum geometry edge/dx exceeds this value.")
    parser.add_argument("--weld-tolerance", type=float, default=1.0e-7,
                        help="STL vertex welding tolerance (default: 1e-7, "
                             "matching the BVH reader).")
    parser.add_argument("--degenerate-tolerance", type=float, default=1.0e-13)
    parser.add_argument("--normal-closure-tolerance", type=float, default=1.0e-8)
    parser.add_argument("--allow-inward", action="store_true")
    parser.add_argument("--require-outward-input", action="store_true",
                        help="Fail instead of relying on the BVH reader's "
                             "automatic orientation correction.")
    parser.add_argument("--check-self-intersections", action="store_true",
                        help="O(N^2) proper-intersection check for 2-D curves.")
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()

    if not args.geometry.is_file():
        parser.error(f"missing geometry: {args.geometry}")
    if any(value <= 0.0 or not math.isfinite(value) for value in args.dx):
        parser.error("--dx values must be finite and positive")
    if args.max_edge_over_dx is not None and args.max_edge_over_dx <= 0.0:
        parser.error("--max-edge-over-dx must be positive")
    if args.weld_tolerance <= 0.0 or not math.isfinite(args.weld_tolerance):
        parser.error("--weld-tolerance must be finite and positive")

    kind = args.kind
    if kind == "auto":
        kind = "stl" if args.geometry.suffix.lower() == ".stl" else "curve"
    result = stl_audit(args.geometry, args) if kind == "stl" else curve_audit(args.geometry, args)
    report = {"file": str(args.geometry.resolve()), **result}

    rendered = json.dumps(report, indent=2)
    print(rendered)
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(rendered + "\n")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
