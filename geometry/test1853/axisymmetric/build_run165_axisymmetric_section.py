#!/usr/bin/env python3
"""Build and audit the Run-165 single-nozzle axisymmetric meridional section.

The generated Cerisse preview geometry is one CCW polygon in ``(r, x)`` order,
SI metres, with implicit last-to-first closure.  It is derived from the same
public-data reconstruction equations as the audited 3-D OCC solid.  The script
also intersects the prepared SI STL with ``y=0`` as an independent tessellation
check; the STL slice is never used to define the 2-D curve.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import struct
from dataclasses import asdict, dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-test1853-axisym")
import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["svg.hashsalt"] = "nasa-test1853-run165-axisymmetric-v1"
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np


HERE = Path(__file__).resolve().parent
TEST1853 = HERE.parent
PARAMETER_FILE = TEST1853 / "geometry_parameters.yaml"
BUILD_AUDIT = TEST1853 / "outputs/reports/build_audit.json"
PREPARED_STL = (
    TEST1853
    / "outputs/ibm_prepared/test1853_run165_single_master_SI"
    / "test1853_run165_single_master_SI.stl"
)


SOURCE_CLASS = {
    "axis_closure": "RZ_NUMERICAL_CLOSURE",
    "throat_jet_patch": "NASA_AS_BUILT_PLUS_DERIVED_X",
    "nozzle_divergent_wall": "NASA_AS_BUILT_PLUS_PLACEMENT_ASSUMPTION",
    "sphere_oml": "NOMINAL_PUBLIC_DATA_RECONSTRUCTION",
    "cone_oml": "NASA_NOMINAL_PLUS_DERIVED_TANGENCY",
    "shoulder_oml": "NASA_RADIUS_PLUS_DERIVED_TANGENCY",
    "cylinder_oml": "NASA_DIMENSIONED",
    "aft_cap": "ARTIFICIAL_IBM_CLOSURE",
}


@dataclass(frozen=True)
class Geometry:
    cone_half_angle_deg: float
    nozzle_divergent_half_angle_deg: float
    nose_radius_mm: float
    sphere_center_x_mm: float
    nominal_nose_tip_x_mm: float
    sphere_cone_x_mm: float
    sphere_cone_r_mm: float
    cone_shoulder_x_mm: float
    cone_shoulder_r_mm: float
    shoulder_radius_mm: float
    shoulder_center_x_mm: float
    shoulder_center_r_mm: float
    body_radius_mm: float
    aft_x_mm: float
    virtual_exit_x_mm: float
    virtual_exit_r_mm: float
    throat_x_mm: float
    throat_r_mm: float


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def derive_geometry(parameters: dict) -> Geometry:
    body = parameters["body"]
    nozzle = parameters["nozzles"]["1"]
    inch = float(parameters["units"]["inch_to_mm"])
    cone_half_angle_deg = float(body["cone_half_angle_deg"])
    alpha = math.radians(cone_half_angle_deg)
    nozzle_divergent_half_angle_deg = float(
        nozzle["divergent_wall_half_angle_deg"]
    )
    beta = math.radians(nozzle_divergent_half_angle_deg)
    rn = float(body["nominal_nose_radius_in"]) * inch
    radius = float(body["maximum_radius_in"]) * inch
    shoulder = float(body["shoulder_fillet_radius_in"]) * inch
    sphere_center_x = rn / math.sin(alpha)
    nose_tip_x = sphere_center_x - rn
    sphere_cone_r = rn * math.cos(alpha)
    sphere_cone_x = rn * math.cos(alpha) ** 2 / math.sin(alpha)
    cone_shoulder_r = radius - shoulder * (1.0 - math.cos(alpha))
    cone_shoulder_x = cone_shoulder_r / math.tan(alpha)
    shoulder_center_r = radius - shoulder
    shoulder_center_x = (
        radius / math.tan(alpha) + shoulder * math.tan(alpha / 2.0)
    )
    virtual_exit_r = 0.5 * float(nozzle["virtual_exit_diameter_mm"])
    throat_r = 0.5 * float(nozzle["throat_diameter_mm"])
    virtual_exit_x = sphere_center_x - math.sqrt(rn**2 - virtual_exit_r**2)
    throat_x = virtual_exit_x + (virtual_exit_r - throat_r) / math.tan(beta)
    return Geometry(
        cone_half_angle_deg=cone_half_angle_deg,
        nozzle_divergent_half_angle_deg=nozzle_divergent_half_angle_deg,
        nose_radius_mm=rn,
        sphere_center_x_mm=sphere_center_x,
        nominal_nose_tip_x_mm=nose_tip_x,
        sphere_cone_x_mm=sphere_cone_x,
        sphere_cone_r_mm=sphere_cone_r,
        cone_shoulder_x_mm=cone_shoulder_x,
        cone_shoulder_r_mm=cone_shoulder_r,
        shoulder_radius_mm=shoulder,
        shoulder_center_x_mm=shoulder_center_x,
        shoulder_center_r_mm=shoulder_center_r,
        body_radius_mm=radius,
        aft_x_mm=float(body["aft_end_x_mm"]),
        virtual_exit_x_mm=virtual_exit_x,
        virtual_exit_r_mm=virtual_exit_r,
        throat_x_mm=throat_x,
        throat_r_mm=throat_r,
    )


def validate_against_build_audit(g: Geometry, audit: dict) -> float:
    run = audit["test1853_run165_single"]
    nozzle = run["nozzles"][0]
    checks = {
        "virtual_exit_x_mm": float(nozzle["virtual_exit_x_mm"]),
        "virtual_exit_r_mm": float(nozzle["virtual_exit_radius_mm"]),
        "throat_x_mm": float(nozzle["throat_x_mm"]),
        "throat_r_mm": float(nozzle["throat_radius_mm"]),
        "aft_x_mm": float(run["occ"]["bounds_mm"][1][0]),
        "body_radius_mm": float(run["occ"]["bounds_mm"][1][1]),
    }
    maximum = 0.0
    for name, expected in checks.items():
        error = abs(float(getattr(g, name)) - expected)
        maximum = max(maximum, error)
        if error > 1.0e-9:
            raise RuntimeError(f"{name}: derived={getattr(g, name)} audit={expected}")
    return maximum


def build_segments(g: Geometry) -> dict[str, np.ndarray]:
    """Return full segment point arrays in (r_mm, x_mm), CCW traversal order."""
    alpha = math.radians(g.cone_half_angle_deg)
    beta = math.radians(g.nozzle_divergent_half_angle_deg)

    x_sphere = np.linspace(g.virtual_exit_x_mm, g.sphere_cone_x_mm, 150)
    r_sphere = np.sqrt(
        np.maximum(
            0.0,
            g.nose_radius_mm**2 - (x_sphere - g.sphere_center_x_mm) ** 2,
        )
    )
    x_shoulder = np.linspace(g.cone_shoulder_x_mm, g.shoulder_center_x_mm, 150)
    r_shoulder = g.shoulder_center_r_mm + np.sqrt(
        np.maximum(
            0.0,
            g.shoulder_radius_mm**2
            - (x_shoulder - g.shoulder_center_x_mm) ** 2,
        )
    )

    # Each array includes both endpoints.  Points are emitted below as
    # segment[:-1], so the next segment owns the shared vertex and the final
    # aft-cap edge is represented by the reader's implicit last-to-first edge.
    return {
        "axis_closure": np.array(
            [[0.0, g.aft_x_mm], [0.0, g.throat_x_mm]], dtype=float
        ),
        "throat_jet_patch": np.array(
            [[0.0, g.throat_x_mm], [g.throat_r_mm, g.throat_x_mm]], dtype=float
        ),
        "nozzle_divergent_wall": np.array(
            [
                [g.throat_r_mm, g.throat_x_mm],
                [g.virtual_exit_r_mm, g.virtual_exit_x_mm],
            ],
            dtype=float,
        ),
        "sphere_oml": np.column_stack((r_sphere, x_sphere)),
        "cone_oml": np.array(
            [
                [g.sphere_cone_r_mm, g.sphere_cone_x_mm],
                [g.cone_shoulder_r_mm, g.cone_shoulder_x_mm],
            ],
            dtype=float,
        ),
        "shoulder_oml": np.column_stack((r_shoulder, x_shoulder)),
        "cylinder_oml": np.array(
            [
                [g.body_radius_mm, g.shoulder_center_x_mm],
                [g.body_radius_mm, g.aft_x_mm],
            ],
            dtype=float,
        ),
        "aft_cap": np.array(
            [[g.body_radius_mm, g.aft_x_mm], [0.0, g.aft_x_mm]], dtype=float
        ),
    }


def assemble_profile(
    segments: dict[str, np.ndarray],
) -> tuple[np.ndarray, list[str]]:
    points: list[np.ndarray] = []
    outgoing: list[str] = []
    for name, values in segments.items():
        for value in values[:-1]:
            points.append(value)
            outgoing.append(name)
    profile = np.asarray(points, dtype=float)
    if len(profile) != len(outgoing):
        raise AssertionError("profile label mismatch")
    if np.any(np.linalg.norm(np.roll(profile, -1, axis=0) - profile, axis=1) <= 1e-14):
        raise RuntimeError("profile contains a zero-length edge")
    return profile, outgoing


def polygon_area_centroid_r(profile_rx_mm: np.ndarray) -> tuple[float, float]:
    r = profile_rx_mm[:, 0]
    x = profile_rx_mm[:, 1]
    r_next = np.roll(r, -1)
    x_next = np.roll(x, -1)
    cross = r * x_next - r_next * x
    area = 0.5 * float(np.sum(cross))
    centroid_r = float(np.sum((r + r_next) * cross) / (6.0 * area))
    return area, centroid_r


def segments_intersect(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> bool:
    def orient(p: np.ndarray, q: np.ndarray, r: np.ndarray) -> float:
        return float((q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]))

    o1, o2 = orient(a, b, c), orient(a, b, d)
    o3, o4 = orient(c, d, a), orient(c, d, b)
    return o1 * o2 < 0.0 and o3 * o4 < 0.0


def polygon_is_simple(profile: np.ndarray) -> bool:
    n = len(profile)
    for i in range(n):
        a, b = profile[i], profile[(i + 1) % n]
        for j in range(i + 2, n):
            if i == 0 and j == n - 1:
                continue
            c, d = profile[j], profile[(j + 1) % n]
            if segments_intersect(a, b, c, d):
                return False
    return True


def circle_sqrt_antiderivative(t: float, radius: float) -> float:
    ratio = max(-1.0, min(1.0, t / radius))
    root = math.sqrt(max(0.0, radius**2 - t**2))
    return 0.5 * (t * root + radius**2 * math.asin(ratio))


def exact_volume_mm3(g: Geometry) -> float:
    """Analytic revolution volume of the current reconstructed solid."""
    x0, x1 = g.virtual_exit_x_mm, g.sphere_cone_x_mm
    t0, t1 = x0 - g.sphere_center_x_mm, x1 - g.sphere_center_x_mm
    sphere = g.nose_radius_mm**2 * (x1 - x0) - (t1**3 - t0**3) / 3.0

    tan_alpha = math.tan(math.radians(g.cone_half_angle_deg))
    cone = tan_alpha**2 * (
        g.cone_shoulder_x_mm**3 - g.sphere_cone_x_mm**3
    ) / 3.0

    t0 = g.cone_shoulder_x_mm - g.shoulder_center_x_mm
    t1 = 0.0
    rc, rs = g.shoulder_center_r_mm, g.shoulder_radius_mm
    shoulder = (
        (rc**2 + rs**2) * (t1 - t0)
        - (t1**3 - t0**3) / 3.0
        + 2.0
        * rc
        * (
            circle_sqrt_antiderivative(t1, rs)
            - circle_sqrt_antiderivative(t0, rs)
        )
    )
    cylinder = g.body_radius_mm**2 * (g.aft_x_mm - g.shoulder_center_x_mm)
    beta = math.radians(g.nozzle_divergent_half_angle_deg)
    cavity = (
        g.virtual_exit_r_mm**3 - g.throat_r_mm**3
    ) / (3.0 * math.tan(beta))
    return math.pi * (sphere + cone + shoulder + cylinder - cavity)


def read_binary_stl_vertices(path: Path) -> tuple[np.ndarray, int]:
    size = path.stat().st_size
    with path.open("rb") as stream:
        header = stream.read(84)
    if len(header) != 84:
        raise RuntimeError(f"short STL header: {path}")
    count = struct.unpack("<I", header[80:84])[0]
    if size != 84 + 50 * count:
        raise RuntimeError(f"expected binary STL size {84 + 50 * count}, got {size}")
    dtype = np.dtype(
        [
            ("normal", "<f4", (3,)),
            ("vertices", "<f4", (3, 3)),
            ("attribute", "<u2"),
        ],
        align=False,
    )
    records = np.memmap(path, dtype=dtype, mode="r", offset=84, shape=(count,))
    return np.asarray(records["vertices"], dtype=np.float64), count


def intersect_stl_y0(
    vertices_m: np.ndarray, plane_tol_m: float = 1.0e-12, weld_tol_m: float = 1.0e-10
) -> tuple[np.ndarray, dict]:
    """Return unique signed (x,z) segments for the y=0 meridional plane."""
    raw: list[tuple[np.ndarray, np.ndarray]] = []
    coplanar_triangles = 0
    for tri in vertices_m:
        distances = tri[:, 1]
        on = np.abs(distances) <= plane_tol_m
        n_on = int(np.count_nonzero(on))
        if n_on == 3:
            coplanar_triangles += 1
            continue
        if n_on == 2:
            ids = np.flatnonzero(on)
            raw.append((tri[ids[0], [0, 2]], tri[ids[1], [0, 2]]))
            continue
        if n_on == 1:
            i_on = int(np.flatnonzero(on)[0])
            others = [i for i in range(3) if i != i_on]
            d0, d1 = distances[others[0]], distances[others[1]]
            if d0 * d1 < 0.0:
                t = -d0 / (d1 - d0)
                cross = tri[others[0]] + t * (tri[others[1]] - tri[others[0]])
                raw.append((tri[i_on, [0, 2]], cross[[0, 2]]))
            continue
        positive = distances > plane_tol_m
        negative = distances < -plane_tol_m
        if not (np.any(positive) and np.any(negative)):
            continue
        points: list[np.ndarray] = []
        for i, j in ((0, 1), (1, 2), (2, 0)):
            di, dj = distances[i], distances[j]
            if di * dj < 0.0:
                t = -di / (dj - di)
                cross = tri[i] + t * (tri[j] - tri[i])
                points.append(cross[[0, 2]])
        if len(points) == 2:
            raw.append((points[0], points[1]))

    vertex_accum: dict[tuple[int, int], list[np.ndarray]] = {}
    edge_keys: set[tuple[tuple[int, int], tuple[int, int]]] = set()
    for a, b in raw:
        if np.linalg.norm(a - b) <= weld_tol_m:
            continue
        ka = tuple(np.rint(a / weld_tol_m).astype(np.int64))
        kb = tuple(np.rint(b / weld_tol_m).astype(np.int64))
        if ka == kb:
            continue
        vertex_accum.setdefault(ka, []).append(a)
        vertex_accum.setdefault(kb, []).append(b)
        edge_keys.add(tuple(sorted((ka, kb))))

    representative = {
        key: np.mean(np.asarray(values), axis=0) for key, values in vertex_accum.items()
    }
    edges = np.asarray(
        [[representative[a], representative[b]] for a, b in sorted(edge_keys)],
        dtype=float,
    )
    degree = {key: 0 for key in representative}
    adjacency = {key: set() for key in representative}
    for a, b in edge_keys:
        degree[a] += 1
        degree[b] += 1
        adjacency[a].add(b)
        adjacency[b].add(a)
    components = 0
    unseen = set(representative)
    while unseen:
        components += 1
        stack = [unseen.pop()]
        while stack:
            current = stack.pop()
            for neighbor in adjacency[current]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    stack.append(neighbor)
    qa = {
        "plane": "y=0",
        "plane_tolerance_m": plane_tol_m,
        "weld_tolerance_m": weld_tol_m,
        "raw_segments": len(raw),
        "unique_segments": len(edge_keys),
        "unique_vertices": len(representative),
        "connected_components": components,
        "minimum_degree": min(degree.values()),
        "maximum_degree": max(degree.values()),
        "non_degree_two_vertices": sum(value != 2 for value in degree.values()),
        "coplanar_triangles": coplanar_triangles,
    }
    if components != 1 or qa["non_degree_two_vertices"] != 0:
        raise RuntimeError(f"STL section is not one closed degree-2 loop: {qa}")
    return edges, qa


def positive_half_segments_mm(edges_xz_m: np.ndarray) -> np.ndarray:
    clipped: list[np.ndarray] = []
    for edge in edges_xz_m:
        a, b = edge.copy()
        za, zb = a[1], b[1]
        if za < 0.0 and zb < 0.0:
            continue
        if za < 0.0 <= zb:
            t = -za / (zb - za)
            a = a + t * (b - a)
            a[1] = 0.0
        elif zb < 0.0 <= za:
            t = -zb / (za - zb)
            b = b + t * (a - b)
            b[1] = 0.0
        if np.linalg.norm(a - b) > 1.0e-12:
            clipped.append(np.asarray([a, b]) * 1000.0)
    return np.asarray(clipped, dtype=float)


def distances_to_segment_sets(
    query_xr_mm: np.ndarray, segment_sets_rx_mm: list[np.ndarray]
) -> np.ndarray:
    result = np.full(len(query_xr_mm), np.inf)
    for values_rx in segment_sets_rx_mm:
        values = values_rx[:, [1, 0]]  # x,r
        a, b = values[:-1], values[1:]
        vector = b - a
        denom = np.sum(vector * vector, axis=1)
        for start in range(0, len(query_xr_mm), 256):
            query = query_xr_mm[start : start + 256]
            delta = query[:, None, :] - a[None, :, :]
            t = np.sum(delta * vector[None, :, :], axis=2) / denom[None, :]
            t = np.clip(t, 0.0, 1.0)
            closest = a[None, :, :] + t[:, :, None] * vector[None, :, :]
            distance = np.linalg.norm(query[:, None, :] - closest, axis=2)
            result[start : start + len(query)] = np.minimum(
                result[start : start + len(query)], np.min(distance, axis=1)
            )
    return result


def write_profile_files(
    output: Path, profile: np.ndarray, outgoing: list[str]
) -> tuple[Path, Path]:
    dat_path = output / "run165_axisymmetric_section_preview.dat"
    csv_path = output / "run165_axisymmetric_section_profile.csv"
    np.savetxt(dat_path, profile * 1.0e-3, fmt="%.15e")
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "point_id",
                "outgoing_segment",
                "source_class",
                "r_m",
                "x_m",
                "r_mm",
                "x_mm",
            ]
        )
        for index, (point, name) in enumerate(zip(profile, outgoing)):
            writer.writerow(
                [
                    index,
                    name,
                    SOURCE_CLASS[name],
                    f"{point[0] * 1.0e-3:.15e}",
                    f"{point[1] * 1.0e-3:.15e}",
                    f"{point[0]:.12f}",
                    f"{point[1]:.12f}",
                ]
            )
    return dat_path, csv_path


def render_section(
    output: Path,
    g: Geometry,
    profile: np.ndarray,
    segments: dict[str, np.ndarray],
    stl_half_segments_mm: np.ndarray,
    stl_max_distance_mm: float,
) -> tuple[Path, Path]:
    colors = {
        "axis_closure": "#616161",
        "throat_jet_patch": "#f28e2b",
        "nozzle_divergent_wall": "#d62728",
        "sphere_oml": "#1769aa",
        "cone_oml": "#1769aa",
        "shoulder_oml": "#2ca02c",
        "cylinder_oml": "#1769aa",
        "aft_cap": "#9467bd",
    }
    fig = plt.figure(figsize=(16, 8.5), constrained_layout=True)
    grid = fig.add_gridspec(1, 2, width_ratios=(2.05, 1.0))
    full = fig.add_subplot(grid[0, 0])
    zoom = fig.add_subplot(grid[0, 1])
    x_profile, r_profile = profile[:, 1], profile[:, 0]
    full.fill(x_profile, r_profile, color="#dfe5eb", zorder=0)
    full.fill(x_profile, -r_profile, color="#dfe5eb", zorder=0)
    zoom.fill(x_profile, r_profile, color="#dfe5eb", zorder=0)

    # Actual triangulated-STL intersection, drawn behind the analytic profile.
    for edge in stl_half_segments_mm:
        full.plot(edge[:, 0], edge[:, 1], color="0.55", lw=0.45, alpha=0.8)
        full.plot(edge[:, 0], -edge[:, 1], color="0.55", lw=0.45, alpha=0.8)
        zoom.plot(edge[:, 0], edge[:, 1], color="0.55", lw=0.5, alpha=0.85)

    for name, values in segments.items():
        x, r = values[:, 1], values[:, 0]
        style = "--" if name == "axis_closure" else "-"
        width = 3.0 if name == "throat_jet_patch" else 2.0
        full.plot(x, r, style, color=colors[name], lw=width, zorder=3)
        full.plot(x, -r, style, color=colors[name], lw=width, zorder=3)
        zoom.plot(x, r, style, color=colors[name], lw=width, zorder=3)

    for axis in (full, zoom):
        axis.axvline(0.0, color="0.35", lw=0.9, ls=":")
        axis.grid(alpha=0.18)
        axis.set_xlabel("x from theoretical sharp point [mm]")
        axis.set_aspect("equal", adjustable="box")
    full.axhline(0.0, color="0.35", lw=0.8)
    full.set_ylabel("signed radius [mm]")
    zoom.set_ylabel("r [mm] (Cerisse RZ half-plane)")
    full.set_xlim(-10.0, g.aft_x_mm + 7.0)
    full.set_ylim(-70.0, 70.0)
    zoom.set_xlim(-4.0, 31.0)
    zoom.set_ylim(-1.0, 67.0)
    full.set_title("Run 165 single-nozzle meridional section")
    zoom.set_title("Axisymmetric IBM profile — forebody/nozzle zoom")

    zoom.annotate(
        f"reconstructed N1 lip\nx={g.virtual_exit_x_mm:.6f}\nRe={g.virtual_exit_r_mm:.5f}",
        (g.virtual_exit_x_mm, g.virtual_exit_r_mm),
        xytext=(6.6, 5.8),
        textcoords="data",
        fontsize=8,
        arrowprops={"arrowstyle": "->", "lw": 0.8},
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 1.5},
    )
    zoom.annotate(
        f"N1 throat face (future jet BC)\nx={g.throat_x_mm:.6f}\nRt={g.throat_r_mm:.5f}",
        (g.throat_x_mm, 0.55 * g.throat_r_mm),
        xytext=(17.0, 1.4),
        textcoords="data",
        fontsize=8,
        arrowprops={"arrowstyle": "->", "lw": 0.8},
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 1.5},
    )
    zoom.annotate(
        "sphere / 70-deg cone",
        (g.sphere_cone_x_mm, g.sphere_cone_r_mm),
        xytext=(8.7, 13.2),
        textcoords="data",
        fontsize=8,
        arrowprops={"arrowstyle": "->", "lw": 0.8},
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 1.5},
    )
    zoom.annotate(
        "R2.54 shoulder / cylinder",
        (g.shoulder_center_x_mm, g.body_radius_mm),
        xytext=(-122, -42),
        textcoords="offset points",
        fontsize=8,
        arrowprops={"arrowstyle": "->", "lw": 0.8},
    )
    full.text(0.0, -68.0, "TSP x=0", ha="center", va="bottom", fontsize=8)
    full.annotate(
        "artificial flat aft closure",
        (g.aft_x_mm, 0.7 * g.body_radius_mm),
        xytext=(-132, -5),
        textcoords="offset points",
        fontsize=8,
        arrowprops={"arrowstyle": "->", "lw": 0.8},
    )

    legend = [
        Patch(facecolor="#dfe5eb", edgecolor="none", label="IBM solid"),
        Line2D([0], [0], color="#1769aa", lw=2, label="outer mold line"),
        Line2D(
            [0],
            [0],
            color="#2ca02c",
            lw=2,
            label="R2.54 shoulder reconstruction",
        ),
        Line2D([0], [0], color="#d62728", lw=2, label="N1 divergent wall"),
        Line2D(
            [0], [0], color="#f28e2b", lw=3, label="N1 throat face (future BC)"
        ),
        Line2D([0], [0], color="#9467bd", lw=2, label="artificial aft cap"),
        Line2D([0], [0], color="#616161", lw=2, ls="--", label="RZ axis closure"),
        Line2D(
            [0],
            [0],
            color="0.55",
            lw=1,
            label="prepared reconstructed SI STL, y=0 slice",
        ),
    ]
    full.legend(handles=legend, loc="lower center", ncol=2, fontsize=8)
    fig.suptitle(
        "NASA Test 1853 Run 165 — AXISYMMETRIC GEOMETRY PREVIEW "
        "(flow BCs not yet configured/tested)\n"
        "public-data reconstruction, angle of attack = 0° only; "
        "max sampled STL-slice endpoint → analytic-profile distance "
        f"{stl_max_distance_mm * 1000.0:.2f} µm (one-sided)",
        fontsize=13,
    )
    png = output / "run165_axisymmetric_meridional_section.png"
    svg = output / "run165_axisymmetric_meridional_section.svg"
    fig.savefig(png, dpi=240)
    # Suppress Matplotlib's run-time dc:date field so the audited SVG hash is
    # reproducible when the same source geometry is regenerated.
    fig.savefig(svg, metadata={"Date": None})
    plt.close(fig)
    return png, svg


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=HERE / "outputs")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    parameters = load_json(PARAMETER_FILE)
    build_audit = load_json(BUILD_AUDIT)
    g = derive_geometry(parameters)
    audit_coordinate_error = validate_against_build_audit(g, build_audit)
    segments = build_segments(g)
    profile, outgoing = assemble_profile(segments)
    area_mm2, centroid_r_mm = polygon_area_centroid_r(profile)
    if area_mm2 <= 0.0:
        raise RuntimeError(f"profile must be CCW, signed area={area_mm2}")
    if not polygon_is_simple(profile):
        raise RuntimeError("profile is self-intersecting")
    if float(np.min(profile[:, 0])) < 0.0:
        raise RuntimeError("RZ profile contains negative radius")

    exact_volume = exact_volume_mm3(g)
    occ_volume = float(build_audit["test1853_run165_single"]["occ"]["volume_mm3"])
    polyline_volume = 2.0 * math.pi * area_mm2 * centroid_r_mm

    stl_vertices, facets = read_binary_stl_vertices(PREPARED_STL)
    stl_edges, stl_section_qa = intersect_stl_y0(stl_vertices)
    stl_half = positive_half_segments_mm(stl_edges)
    stl_query = np.unique(np.round(stl_half.reshape((-1, 2)), 10), axis=0)
    physical_segments = [
        values
        for name, values in segments.items()
        if name != "axis_closure"
    ]
    distances = distances_to_segment_sets(stl_query, physical_segments)
    dat_path, csv_path = write_profile_files(output, profile, outgoing)
    png_path, svg_path = render_section(
        output,
        g,
        profile,
        segments,
        stl_half,
        float(np.max(distances)),
    )

    metadata = {
        "status": "PREVIEW_GEOMETRY_QA_PASS",
        "purpose": (
            "Run-165 single-nozzle, zero-angle-of-attack axisymmetric "
            "meridional section"
        ),
        "solver_status": (
            "The polygon is reader-compatible by construction but has not yet been "
            "advanced as a dedicated Cerisse RZ case."
        ),
        "coordinate_convention": {
            "dat_columns": ["r_m", "x_from_TSP_m"],
            "axis": "r=0",
            "flow": "freestream +x; retropropulsive jet -x",
            "winding": "CCW",
            "closure": "implicit last vertex to first vertex",
        },
        "scope_limitations": [
            "Valid only for the Run 165 single center nozzle at zero angle of attack.",
            "The nominal R25.4 nose and N1 virtual-exit axial placement are public-data reconstruction assumptions, not unpublished as-built CAD.",
            "The aft cap and r=0 edge are numerical closures; exclude them from aerodynamic wall loads.",
            "No complete convergent passage, plenum, sting, pressure-port holes, or nonzero-angle physics is added.",
            "RZ total forces require 2*pi*r weighting in post-processing; the current built-in 2-D surface measure is line length.",
        ],
        "source_sha256": {
            "geometry_parameters.yaml": sha256(PARAMETER_FILE),
            "build_audit.json": sha256(BUILD_AUDIT),
            "prepared_run165_SI.stl": sha256(PREPARED_STL),
            "generator": sha256(Path(__file__).resolve()),
        },
        "geometry_parameters": asdict(g),
        "profile": {
            "input_vertices": len(profile),
            "signed_area_mm2": area_mm2,
            "centroid_r_mm": centroid_r_mm,
            "minimum_r_mm": float(np.min(profile[:, 0])),
            "maximum_r_mm": float(np.max(profile[:, 0])),
            "minimum_x_mm": float(np.min(profile[:, 1])),
            "maximum_x_mm": float(np.max(profile[:, 1])),
            "simple_polygon": True,
            "source_classes": SOURCE_CLASS,
            "sphere_samples": len(segments["sphere_oml"]),
            "shoulder_samples": len(segments["shoulder_oml"]),
        },
        "volume_cross_check": {
            "analytic_axisymmetric_mm3": exact_volume,
            "audited_OCC_mm3": occ_volume,
            "analytic_minus_OCC_mm3": exact_volume - occ_volume,
            "analytic_OCC_relative_difference": abs(exact_volume - occ_volume)
            / occ_volume,
            "preview_polyline_mm3": polyline_volume,
            "polyline_OCC_relative_difference": abs(polyline_volume - occ_volume)
            / occ_volume,
        },
        "audit_coordinate_max_error_mm": audit_coordinate_error,
        "stl_section_cross_check": {
            "source_facets": facets,
            "distance_metric_definition": (
                "Maximum Euclidean distance from each unique endpoint of the "
                "positive-half y=0 prepared reconstructed SI-STL section to the "
                "analytic profile polyline. This is a one-sided section metric, "
                "not a global STL-surface/CAD Hausdorff error."
            ),
            **stl_section_qa,
            "positive_half_segments": len(stl_half),
            "unique_positive_half_endpoints": len(stl_query),
            "STL_endpoint_to_analytic_profile_max_mm": float(np.max(distances)),
            "STL_endpoint_to_analytic_profile_rms_mm": float(
                np.sqrt(np.mean(distances**2))
            ),
        },
        "outputs_sha256": {
            dat_path.name: sha256(dat_path),
            csv_path.name: sha256(csv_path),
            png_path.name: sha256(png_path),
            svg_path.name: sha256(svg_path),
        },
    }
    metadata_path = output / "run165_axisymmetric_section_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {dat_path}")
    print(f"wrote {csv_path}")
    print(f"wrote {png_path}")
    print(f"wrote {svg_path}")
    print(f"wrote {metadata_path}")
    print(
        "QA "
        f"vertices={len(profile)} area={area_mm2:.9f}mm2 "
        f"volume_rel={abs(exact_volume-occ_volume)/occ_volume:.3e} "
        f"stl_max={np.max(distances)*1000.0:.3f}um "
        f"section_edges={stl_section_qa['unique_segments']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
