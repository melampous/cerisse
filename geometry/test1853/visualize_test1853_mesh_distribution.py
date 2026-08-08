#!/usr/bin/env python3
"""Visualize Test 1853 surface tessellation and a proposed IBM AMR zoning.

The two quantities are intentionally kept separate:

* ``h_STL`` is a triangle-edge metric on the fixed geometry surface.
* ``dx`` is the Cartesian AMReX cell spacing in the fluid volume.

The AMR figure is a design proposal, not a plot of a computed AMR hierarchy.
An actual hierarchy must be rendered from a Test 1853 plotfile after regridding.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-test1853-mesh-distribution")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np
import pyvista as pv


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "outputs"
MODEL = "test1853_run262_263_tri_n3at120"

ZONE_NAMES = {
    0: "active nozzles",
    1: "forebody",
    2: "R2.54 shoulder",
    3: "forward cylinder",
    4: "35-60 mm size transition",
    5: "aft cylinder",
    6: "artificial aft cap",
}


def load_parameters() -> dict[str, Any]:
    return json.loads((HERE / "geometry_parameters.yaml").read_text(encoding="utf-8"))


def derived_body(parameters: dict[str, Any]) -> dict[str, float]:
    source = parameters["body"]["analytic_reconstruction_in"]
    inch = float(parameters["units"]["inch_to_mm"])
    return {
        "sphere_center_x": float(source["sphere_center_x"]) * inch,
        "nose_tip_x": float(source["nose_tip_x"]) * inch,
        "sphere_cone_tangent_x": float(source["sphere_cone_tangent_x"]) * inch,
        "sphere_cone_tangent_r": float(source["sphere_cone_tangent_r"]) * inch,
        "cone_shoulder_tangent_x": float(source["cone_shoulder_tangent_x"]) * inch,
        "cone_shoulder_tangent_r": float(source["cone_shoulder_tangent_r"]) * inch,
        "shoulder_center_x": float(source["shoulder_circle_center_x"]) * inch,
        "shoulder_center_r": float(source["shoulder_circle_center_r"]) * inch,
        "nose_radius": float(parameters["body"]["nominal_nose_radius_mm"]),
        "shoulder_radius": float(parameters["body"]["shoulder_fillet_radius_mm"]),
        "body_radius": float(parameters["body"]["maximum_radius_mm"]),
        "aft_end_x": float(parameters["body"]["aft_end_x_mm"]),
        "alpha": math.radians(float(parameters["body"]["cone_half_angle_deg"])),
    }


def body_radius_at_x(x: np.ndarray, body: dict[str, float]) -> np.ndarray:
    result = np.full_like(x, np.nan, dtype=float)
    sphere = (x >= body["nose_tip_x"]) & (x <= body["sphere_cone_tangent_x"])
    result[sphere] = np.sqrt(
        np.maximum(
            0.0,
            body["nose_radius"] ** 2
            - (x[sphere] - body["sphere_center_x"]) ** 2,
        )
    )
    cone = (x > body["sphere_cone_tangent_x"]) & (
        x <= body["cone_shoulder_tangent_x"]
    )
    result[cone] = x[cone] * math.tan(body["alpha"])
    shoulder = (x > body["cone_shoulder_tangent_x"]) & (
        x <= body["shoulder_center_x"]
    )
    result[shoulder] = body["shoulder_center_r"] + np.sqrt(
        np.maximum(
            0.0,
            body["shoulder_radius"] ** 2
            - (x[shoulder] - body["shoulder_center_x"]) ** 2,
        )
    )
    cylinder = (x > body["shoulder_center_x"]) & (x <= body["aft_end_x"])
    result[cylinder] = body["body_radius"]
    return result


def throat_data(output: Path) -> list[dict[str, float]]:
    metadata_path = (
        output
        / "ibm_prepared"
        / f"{MODEL}_master_SI"
        / f"{MODEL}_master_SI_metadata.json"
    )
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    result = []
    for item in metadata["active_nozzle_throats"]:
        center = [1000.0 * float(value) for value in item["throat_center_world_m"]]
        result.append(
            {
                "number": int(item["number"]),
                "x": center[0],
                "y": center[1],
                "z": center[2],
                "radius": 1000.0 * float(item["throat_radius_m"]),
            }
        )
    return result


def face_arrays(mesh: pv.PolyData) -> dict[str, np.ndarray]:
    faces = np.asarray(mesh.faces).reshape(-1, 4)[:, 1:]
    vertices = np.asarray(mesh.points, dtype=float)
    triangles = vertices[faces]
    edge = np.stack(
        [
            np.linalg.norm(triangles[:, 1] - triangles[:, 0], axis=1),
            np.linalg.norm(triangles[:, 2] - triangles[:, 1], axis=1),
            np.linalg.norm(triangles[:, 0] - triangles[:, 2], axis=1),
        ],
        axis=1,
    )
    cross = np.cross(
        triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]
    )
    twice_area = np.linalg.norm(cross, axis=1)
    normal = cross / np.maximum(twice_area[:, None], np.finfo(float).tiny)
    area = 0.5 * twice_area
    return {
        "edge": edge,
        "hmax": np.max(edge, axis=1),
        "hmean": np.mean(edge, axis=1),
        "area": area,
        "heq": np.sqrt(4.0 * area / math.sqrt(3.0)),
        "centroid": np.mean(triangles, axis=1),
        "normal": normal,
    }


def classify_zones(
    arrays: dict[str, np.ndarray], throats: list[dict[str, float]], body: dict[str, float]
) -> np.ndarray:
    centroid = arrays["centroid"]
    normal = arrays["normal"]
    nearest_axis = np.full(len(centroid), np.inf)
    radial_normal = np.zeros(len(centroid))
    for throat in throats:
        dy = centroid[:, 1] - throat["y"]
        dz = centroid[:, 2] - throat["z"]
        distance = np.sqrt(dy * dy + dz * dz)
        dot = (normal[:, 1] * dy + normal[:, 2] * dz) / np.maximum(
            distance, np.finfo(float).tiny
        )
        update = distance < nearest_axis
        nearest_axis[update] = distance[update]
        radial_normal[update] = dot[update]

    nozzle = (
        (centroid[:, 0] < 24.5)
        & (nearest_axis < 7.5)
        & (radial_normal < -0.2)
    )
    for throat in throats:
        transverse = np.sqrt(
            (centroid[:, 1] - throat["y"]) ** 2
            + (centroid[:, 2] - throat["z"]) ** 2
        )
        nozzle |= (
            (np.abs(centroid[:, 0] - throat["x"]) < 0.05)
            & (transverse < throat["radius"] + 0.1)
        )

    zone = np.full(len(centroid), -1, dtype=np.int32)
    zone[nozzle] = 0
    remaining = ~nozzle
    zone[remaining & (centroid[:, 0] < body["cone_shoulder_tangent_x"])] = 1
    zone[
        remaining
        & (centroid[:, 0] >= body["cone_shoulder_tangent_x"])
        & (centroid[:, 0] < body["shoulder_center_x"])
    ] = 2
    zone[
        remaining
        & (centroid[:, 0] >= body["shoulder_center_x"])
        & (centroid[:, 0] < 35.0)
    ] = 3
    flat_cap = (
        np.abs(centroid[:, 0] - body["aft_end_x"]) <= 1.0e-3
    ) & (normal[:, 0] >= 0.999)
    zone[(centroid[:, 0] >= 35.0) & (centroid[:, 0] < 60.0)] = 4
    zone[(centroid[:, 0] >= 60.0) & ~flat_cap] = 5
    zone[flat_cap] = 6
    if np.any(zone < 0):
        raise RuntimeError(f"zone classification left {np.count_nonzero(zone < 0)} faces")
    return zone


def attach_metrics(
    mesh: pv.PolyData, arrays: dict[str, np.ndarray], zone: np.ndarray
) -> None:
    mesh.cell_data["hmax_mm"] = arrays["hmax"].astype(np.float32)
    mesh.cell_data["hmean_mm"] = arrays["hmean"].astype(np.float32)
    mesh.cell_data["equilateral_area_spacing_mm"] = arrays["heq"].astype(np.float32)
    mesh.cell_data["area_mm2"] = arrays["area"].astype(np.float32)
    mesh.cell_data["centroid_x_mm"] = arrays["centroid"][:, 0].astype(np.float32)
    mesh.cell_data["visualization_zone_id"] = zone


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def distribution(values: np.ndarray) -> dict[str, float]:
    percentiles = np.percentile(values, [0, 50, 90, 95, 99, 100])
    return {
        key: float(value)
        for key, value in zip(
            ["minimum", "p50", "p90", "p95", "p99", "maximum"], percentiles
        )
    }


def build_report(
    mesh: pv.PolyData,
    arrays: dict[str, np.ndarray],
    zone: np.ndarray,
    output: Path,
) -> dict[str, Any]:
    regions: dict[str, Any] = {}
    for zone_id, name in ZONE_NAMES.items():
        selected = zone == zone_id
        regions[name] = {
            "zone_id": zone_id,
            "faces": int(np.count_nonzero(selected)),
            "fraction_of_faces": float(np.mean(selected)),
            "surface_area_mm2": float(np.sum(arrays["area"][selected])),
            "face_max_edge_mm": distribution(arrays["hmax"][selected]),
            "all_face_edges_mm": distribution(arrays["edge"][selected].reshape(-1)),
        }

    hmax = arrays["hmax"]
    classes = [0.0, 0.5, 0.75, 1.0, 1.5, 2.25, 3.0, np.inf]
    class_report = []
    for lower, upper in zip(classes[:-1], classes[1:]):
        selected = (hmax >= lower) & (hmax < upper)
        class_report.append(
            {
                "lower_mm_inclusive": lower,
                "upper_mm_exclusive": None if np.isinf(upper) else upper,
                "faces": int(np.count_nonzero(selected)),
                "fraction": float(np.mean(selected)),
            }
        )

    finest = 0.20833333333333334
    proposed_levels = {
        "L3_nozzle_forebody_shoulder_x_le_40_mm": {
            "level": 3,
            "dx_mm": finest,
            "relative_volume_cell_density": 1.0,
        },
        "L2_forward_cylinder_40_to_70_mm": {
            "level": 2,
            "dx_mm": 2.0 * finest,
            "relative_volume_cell_density": 1.0 / 8.0,
        },
        "L1_mid_cylinder_70_to_130_mm": {
            "level": 1,
            "dx_mm": 4.0 * finest,
            "relative_volume_cell_density": 1.0 / 64.0,
        },
        "L0_aft_cylinder_and_cap_x_gt_130_mm": {
            "level": 0,
            "dx_mm": 8.0 * finest,
            "relative_volume_cell_density": 1.0 / 512.0,
        },
    }
    return {
        "status": "PASS",
        "model": MODEL,
        "triangles": int(mesh.n_cells),
        "vertices": int(mesh.n_points),
        "provenance": {
            "source_stl": str((output / "stl" / f"{MODEL}.stl").relative_to(HERE)),
            "source_stl_sha256": sha256(output / "stl" / f"{MODEL}.stl"),
            "generator": str(Path(__file__).resolve().relative_to(HERE)),
            "generator_sha256": sha256(Path(__file__).resolve()),
            "interactive_vtp_sha256": sha256(
                output / "visualization" / f"{MODEL}_surface_mesh_metrics.vtp"
            ),
            "numpy_version": np.__version__,
            "pyvista_version": pv.__version__,
            "matplotlib_version": matplotlib.__version__,
        },
        "surface_metric": {
            "name": "h_STL",
            "definition": "maximum of the three edge lengths of each STL triangle",
            "units": "mm",
            "regions": regions,
            "hmax_classes": class_report,
        },
        "proposed_ibm_amr": {
            "status": "DESIGN PROPOSAL, NOT A COMPUTED AMR HIERARCHY",
            "coordinate": "body-frame x measured from TSP; length in mm",
            "levels": proposed_levels,
            "nozzle_near_field": (
                "L3 within 10 mm and an L2 parent envelope within 12 mm of each "
                "discrete nozzle axis for -20 <= x <= throat+5 mm"
            ),
            "plume_extension": "L2 within 12 mm of each discrete nozzle axis for -80 <= x < -20 mm",
            "meridional_figure_note": (
                "The x-r panel is an azimuthal projection for visualization only; "
                "it must not be implemented as an axisymmetric annulus or torus."
            ),
            "bow_shock": "geometry-anchored L1 startup envelope plus density-gradient tracking capped at L2",
            "aft_cap_option": "keep only the outer 8-10 mm annulus at L1 if rim/wake sensitivity requires it",
            "implementation_blocker": (
                "Current ib.amr_support_buffer=1 seeds every IBM interface on every "
                "level, so user_tagging alone cannot impose the proposed axial level cap."
            ),
            "required_change": (
                "Use one body-frame local_max_level policy for all user/flow tags and the "
                "baseline support seed; allow support-closure promotion wherever another "
                "tag creates a child grid near IBM, retain runtime halo dilation and strict "
                "support audits, and evaluate gradient tags only on valid-fluid stencils."
            ),
        },
        "files": {
            "interactive_surface_vtp": str(
                (output / "visualization" / f"{MODEL}_surface_mesh_metrics.vtp").relative_to(
                    HERE
                )
            ),
            "surface_overview_png": str(
                (output / "figures" / "tri_surface_mesh_density_overview.png").relative_to(
                    HERE
                )
            ),
            "surface_axial_png": str(
                (output / "figures" / "tri_surface_mesh_density_axial.png").relative_to(
                    HERE
                )
            ),
            "surface_forebody_nozzles_png": str(
                (
                    output
                    / "figures"
                    / "tri_surface_mesh_density_forebody_nozzles.png"
                ).relative_to(HERE)
            ),
            "surface_local_details_png": str(
                (
                    output
                    / "figures"
                    / "tri_surface_mesh_density_local_details.png"
                ).relative_to(HERE)
            ),
            "single_tri_surface_comparison_png": str(
                (
                    output
                    / "figures"
                    / "single_tri_surface_mesh_density_comparison.png"
                ).relative_to(HERE)
            ),
            "ibm_dx_proposal_png": str(
                (output / "figures" / "proposed_ibm_dx_distribution.png").relative_to(
                    HERE
                )
            ),
        },
    }


def add_surface(
    plotter: pv.Plotter,
    mesh: pv.DataSet,
    *,
    scalar_bar: bool,
    edges: bool,
    edge_opacity: float = 0.18,
) -> None:
    kwargs: dict[str, Any] = {
        "scalars": "hmax_mm",
        "preference": "cell",
        "cmap": "viridis_r",
        "clim": (0.35, 3.0),
        "log_scale": True,
        "smooth_shading": False,
        "interpolate_before_map": False,
        "show_edges": edges,
        "edge_color": "#202020",
        "line_width": 0.3,
        "edge_opacity": edge_opacity,
        "show_scalar_bar": scalar_bar,
        "lighting": True,
        "ambient": 0.35,
        "diffuse": 0.65,
    }
    if scalar_bar:
        kwargs["scalar_bar_args"] = {
            "title": "hmax [mm]",
            "n_labels": 6,
            "vertical": True,
            "fmt": "%.2g",
        }
    plotter.add_mesh(mesh, **kwargs)


def set_camera(
    plotter: pv.Plotter,
    position: tuple[float, float, float],
    focal: tuple[float, float, float],
    up: tuple[float, float, float] = (0.0, 0.0, 1.0),
    zoom: float = 1.0,
) -> None:
    plotter.camera_position = [position, focal, up]
    plotter.enable_parallel_projection()
    if zoom != 1.0:
        plotter.camera.zoom(zoom)


def render_surface_overview(mesh: pv.PolyData, figures: Path) -> None:
    fore = mesh.extract_cells(np.asarray(mesh.cell_data["centroid_x_mm"]) <= 40.0)
    shoulder = mesh.extract_cells(
        (np.asarray(mesh.cell_data["centroid_x_mm"]) >= 15.0)
        & (np.asarray(mesh.cell_data["centroid_x_mm"]) <= 42.0)
    )
    aft = mesh.extract_cells(np.asarray(mesh.cell_data["centroid_x_mm"]) >= 225.0)

    plotter = pv.Plotter(
        shape=(2, 2), off_screen=True, window_size=(2200, 1450), border=False
    )
    plotter.set_background("white")

    plotter.subplot(0, 0)
    add_surface(plotter, mesh, scalar_bar=True, edges=False)
    plotter.add_text("A. Full body — actual STL h_STL", font_size=13, color="black")
    set_camera(plotter, (135.0, -520.0, 0.0), (135.0, 0.0, 0.0))

    plotter.subplot(0, 1)
    add_surface(plotter, fore, scalar_bar=False, edges=True)
    plotter.add_text(
        "B. Nozzles + forebody — true triangle edges", font_size=13, color="black"
    )
    set_camera(plotter, (-190.0, -190.0, 145.0), (12.0, 0.0, 0.0), zoom=1.15)

    plotter.subplot(1, 0)
    add_surface(plotter, shoulder, scalar_bar=False, edges=True)
    plotter.add_text(
        "C. Cone / R2.54 shoulder / cylinder transition",
        font_size=13,
        color="black",
    )
    set_camera(plotter, (27.0, -340.0, 0.0), (27.0, 0.0, 0.0), zoom=1.18)

    plotter.subplot(1, 1)
    add_surface(plotter, aft, scalar_bar=False, edges=True)
    plotter.add_text(
        "D. Aft cylinder + artificial flat cap", font_size=13, color="black"
    )
    set_camera(plotter, (520.0, -210.0, 135.0), (252.0, 0.0, 0.0), zoom=1.05)

    plotter.screenshot(figures / "tri_surface_mesh_density_overview.png")
    plotter.close()

    closeup = pv.Plotter(
        shape=(1, 2), off_screen=True, window_size=(2200, 1050), border=False
    )
    closeup.set_background("white")
    closeup.subplot(0, 0)
    add_surface(closeup, fore, scalar_bar=True, edges=True)
    closeup.add_text("Upstream isometric", font_size=14, color="black")
    set_camera(closeup, (-190.0, -190.0, 145.0), (12.0, 0.0, 0.0), zoom=1.2)
    closeup.subplot(0, 1)
    add_surface(closeup, fore, scalar_bar=False, edges=True)
    closeup.add_text("Looking downstream (+x)", font_size=14, color="black")
    set_camera(closeup, (-270.0, 0.0, 0.0), (13.0, 0.0, 0.0), zoom=1.15)
    closeup.screenshot(figures / "tri_surface_mesh_density_forebody_nozzles.png")
    closeup.close()

    centroid_x = np.asarray(mesh.cell_data["centroid_x_mm"])
    centers = np.asarray(mesh.cell_centers().points)
    nozzle_2 = mesh.extract_cells(
        (centroid_x <= 25.0)
        & (np.sqrt(centers[:, 1] ** 2 + (centers[:, 2] - 31.75) ** 2) <= 10.0)
    )
    shoulder_patch = mesh.extract_cells(
        (centroid_x >= 18.0)
        & (centroid_x <= 35.0)
        & (centers[:, 1] < 0.0)
        & (np.abs(centers[:, 2]) <= 22.0)
    )
    cylinder_patch = mesh.extract_cells(
        (centroid_x >= 100.0)
        & (centroid_x <= 135.0)
        & (centers[:, 1] < 0.0)
        & (np.abs(centers[:, 2]) <= 18.0)
    )
    cap_patch = mesh.extract_cells(
        (centroid_x >= 267.0)
        & (np.sqrt(centers[:, 1] ** 2 + centers[:, 2] ** 2) <= 22.0)
    )

    detail = pv.Plotter(
        shape=(2, 2), off_screen=True, window_size=(2000, 1500), border=False
    )
    detail.set_background("white")
    detail.subplot(0, 0)
    add_surface(detail, nozzle_2, scalar_bar=True, edges=True, edge_opacity=0.55)
    detail.add_text("A. N2 nozzle local triangles", font_size=13, color="black")
    set_camera(detail, (-70.0, 0.0, 31.75), (17.0, 0.0, 31.75), zoom=1.45)
    detail.subplot(0, 1)
    add_surface(detail, shoulder_patch, scalar_bar=False, edges=True, edge_opacity=0.55)
    detail.add_text("B. R2.54 shoulder patch", font_size=13, color="black")
    set_camera(detail, (27.0, -230.0, 0.0), (27.0, -61.0, 0.0), zoom=1.55)
    detail.subplot(1, 0)
    add_surface(detail, cylinder_patch, scalar_bar=False, edges=True, edge_opacity=0.65)
    detail.add_text("C. Aft-cylinder wall patch", font_size=13, color="black")
    set_camera(detail, (117.5, -260.0, 0.0), (117.5, -63.5, 0.0), zoom=1.25)
    detail.subplot(1, 1)
    add_surface(detail, cap_patch, scalar_bar=False, edges=True, edge_opacity=0.65)
    detail.add_text("D. Flat-cap centre patch", font_size=13, color="black")
    set_camera(detail, (390.0, 0.0, 0.0), (268.0, 0.0, 0.0), zoom=1.3)
    detail.screenshot(figures / "tri_surface_mesh_density_local_details.png")
    detail.close()


def render_configuration_comparison(
    single: pv.PolyData, tri: pv.PolyData, figures: Path
) -> None:
    single_front = single.extract_cells(
        np.asarray(single.cell_data["centroid_x_mm"]) <= 40.0
    )
    tri_front = tri.extract_cells(np.asarray(tri.cell_data["centroid_x_mm"]) <= 40.0)
    plotter = pv.Plotter(
        shape=(1, 2), off_screen=True, window_size=(2200, 1050), border=False
    )
    plotter.set_background("white")
    plotter.subplot(0, 0)
    add_surface(plotter, single_front, scalar_bar=True, edges=True)
    plotter.add_text("Run 165 — centre N1 active", font_size=14, color="black")
    set_camera(plotter, (-270.0, 0.0, 0.0), (13.0, 0.0, 0.0), zoom=1.15)
    plotter.subplot(0, 1)
    add_surface(plotter, tri_front, scalar_bar=False, edges=True)
    plotter.add_text("Runs 262/263 — N2/N3/N4 active", font_size=14, color="black")
    set_camera(plotter, (-270.0, 0.0, 0.0), (13.0, 0.0, 0.0), zoom=1.15)
    plotter.screenshot(figures / "single_tri_surface_mesh_density_comparison.png")
    plotter.close()


def render_axial_profile(
    arrays: dict[str, np.ndarray], body: dict[str, float], figures: Path
) -> None:
    x = arrays["centroid"][:, 0]
    h = arrays["hmax"]
    area = arrays["area"]
    lower = math.floor(float(np.min(x)))
    upper = math.ceil(float(np.max(x)))
    edges = np.arange(lower, upper + 1.0, 1.0)
    centers = 0.5 * (edges[:-1] + edges[1:])
    p50 = np.full_like(centers, np.nan)
    p95 = np.full_like(centers, np.nan)
    p100 = np.full_like(centers, np.nan)
    density = np.full_like(centers, np.nan)
    for index in range(len(centers)):
        selected = (x >= edges[index]) & (x < edges[index + 1])
        if not np.any(selected):
            continue
        p50[index], p95[index], p100[index] = np.percentile(
            h[selected], [50, 95, 100]
        )
        density[index] = np.count_nonzero(selected) / np.sum(area[selected])

    target = np.full_like(centers, 2.0)
    target[centers <= 35.0] = 0.4
    transition = (centers > 35.0) & (centers < 60.0)
    target[transition] = 0.4 + (centers[transition] - 35.0) * (1.6 / 25.0)

    fig, axes = plt.subplots(2, 1, figsize=(14.5, 8.8), sharex=True, constrained_layout=True)
    axes[0].plot(centers, p50, label="face hmax p50", lw=1.6, color="#1b9e77")
    axes[0].plot(centers, p95, label="face hmax p95", lw=1.3, color="#d95f02")
    axes[0].plot(centers, p100, label="face hmax maximum", lw=0.9, color="#7570b3")
    axes[0].plot(
        centers,
        target,
        label="Gmsh requested size lc(x)",
        lw=2.0,
        ls="--",
        color="black",
    )
    axes[0].set_ylabel("surface triangle scale [mm]")
    axes[0].set_ylim(0.25, 3.05)
    axes[0].grid(alpha=0.2)
    axes[0].legend(ncol=2, fontsize=9)
    axes[0].set_title(
        "Actual Test 1853 tri STL surface density along body axis (1 mm bins)"
    )

    axes[1].plot(centers, density, color="#377eb8", lw=1.4)
    axes[1].set_yscale("log")
    axes[1].set_ylabel("triangles / mm$^2$ in axial bin")
    axes[1].set_xlabel("body-frame x from theoretical sharp point [mm]")
    axes[1].grid(alpha=0.2, which="both")

    stations = [
        (body["nose_tip_x"], "nose tip"),
        (body["sphere_cone_tangent_x"], "sphere/cone"),
        (body["cone_shoulder_tangent_x"], "cone/shoulder"),
        (body["shoulder_center_x"], "shoulder/cylinder"),
        (35.0, "fine target ends"),
        (60.0, "aft target begins"),
        (body["aft_end_x"], "flat cap"),
    ]
    x_offsets = [-9, 9, -8, 8, 0, 0, 0]
    for (station, label), x_offset in zip(stations, x_offsets):
        for axis in axes:
            axis.axvline(station, color="0.42", lw=0.65, ls=":" if station not in (35.0, 60.0) else "--")
        axes[0].annotate(
            label,
            xy=(station, 3.02),
            xytext=(x_offset, -3),
            textcoords="offset points",
            rotation=90,
            va="top",
            ha="left",
            fontsize=7,
            color="0.25",
        )

    fig.savefig(figures / "tri_surface_mesh_density_axial.png", dpi=220)
    plt.close(fig)


def target_level_map(
    body: dict[str, float], throats: list[dict[str, float]]
) -> tuple[np.ndarray, np.ndarray, np.ma.MaskedArray, np.ndarray]:
    x = np.linspace(-180.0, 290.0, 1500)
    radius = np.linspace(0.0, 145.0, 520)
    xx, rr = np.meshgrid(x, radius)
    level = np.zeros_like(xx, dtype=np.int8)

    # Geometry-anchored startup envelope for the bow shock and plume.
    startup = (xx >= -175.0) & (xx <= 65.0) & (rr <= 120.0)
    level[startup] = np.maximum(level[startup], 1)

    body_r = body_radius_at_x(x, body)
    body_r_2d = np.broadcast_to(body_r, xx.shape)
    valid = np.isfinite(body_r_2d)
    distance = np.abs(rr - body_r_2d)

    # Surface level caps.  Band widths are schematic; production widths are
    # computed from the IBM interpolation halo and verified by strict audit.
    front_surface = valid & (xx <= 40.0) & (distance <= 8.0)
    level[front_surface] = 3
    forward_cylinder = valid & (xx > 40.0) & (xx <= 70.0) & (distance <= 9.0)
    level[forward_cylinder] = np.maximum(level[forward_cylinder], 2)
    mid_cylinder = valid & (xx > 70.0) & (xx <= 130.0) & (distance <= 12.0)
    level[mid_cylinder] = np.maximum(level[mid_cylinder], 1)

    # Optional L1 annulus around the artificial aft rim; cap interior remains L0.
    rim = (
        (np.abs(xx - body["aft_end_x"]) <= 10.0)
        & (rr >= body["body_radius"] - 10.0)
        & (rr <= body["body_radius"] + 10.0)
    )
    level[rim] = np.maximum(level[rim], 1)

    # This is the tri-peripheral configuration: all three active axes lie on
    # the same r=31.75 mm pitch circle.  Do not draw a fictitious centre jet.
    radial_axes = [float(np.mean([math.hypot(t["y"], t["z"]) for t in throats]))]
    max_throat_x = max(t["x"] for t in throats)
    for axis_radius in radial_axes:
        parent_envelope = (
            (xx >= -80.0)
            & (xx <= max_throat_x + 5.0)
            & (np.abs(rr - axis_radius) <= 12.0)
        )
        level[parent_envelope] = np.maximum(level[parent_envelope], 2)
        nozzle = (
            (xx >= -20.0)
            & (xx <= max_throat_x + 5.0)
            & (np.abs(rr - axis_radius) <= 10.0)
        )
        level[nozzle] = 3

    # Mask solid body, except for a schematic cut through peripheral nozzle N2.
    solid = valid & (rr < body_r_2d)
    n2 = min(throats, key=lambda item: abs(item["y"]))
    cavity = (
        (xx >= 8.0)
        & (xx <= n2["x"] + 0.2)
        & (np.abs(rr - math.hypot(n2["y"], n2["z"])) <= 6.5)
    )
    solid &= ~cavity
    masked = np.ma.masked_where(solid, level)
    return xx, rr, masked, body_r


def render_dx_proposal(
    parameters: dict[str, Any],
    body: dict[str, float],
    throats: list[dict[str, float]],
    figures: Path,
) -> None:
    xx, rr, level, body_r = target_level_map(body, throats)
    colors = ["#5e4fa2", "#3288bd", "#99d594", "#fdae61"]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)
    finest = 0.20833333333333334
    dx_by_level = {0: 8.0 * finest, 1: 4.0 * finest, 2: 2.0 * finest, 3: finest}

    fig = plt.figure(figsize=(17.0, 10.2), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, width_ratios=[2.45, 1.0], height_ratios=[1.15, 1.0])
    side = fig.add_subplot(grid[:, 0])
    front = fig.add_subplot(grid[0, 1])
    axial = fig.add_subplot(grid[1, 1])

    image = side.pcolormesh(xx, rr, level, cmap=cmap, norm=norm, shading="auto")
    x_profile = xx[0]
    side.fill_between(
        x_profile,
        0.0,
        np.nan_to_num(body_r, nan=0.0),
        where=np.isfinite(body_r),
        color="#d9d9d9",
        edgecolor="black",
        linewidth=0.8,
        label="solid model",
        zorder=4,
    )
    # Re-show the N2 cavity schematically after filling the body.
    n2 = min(throats, key=lambda item: abs(item["y"]))
    n2_r = math.hypot(n2["y"], n2["z"])
    side.fill_betweenx(
        [n2_r - 6.5, n2_r + 6.5],
        8.0,
        n2["x"],
        color=colors[3],
        alpha=0.95,
        zorder=5,
    )
    side.plot([8.0, n2["x"]], [n2_r, n2_r], color="black", lw=0.8, ls="--", zorder=6)
    side.text(-77.0, n2_r + 15.0, "peripheral plume L2", fontsize=8)
    side.text(
        -55.0,
        n2_r + 13.0,
        "nozzle near-field L3",
        fontsize=8,
        zorder=7,
    )
    side.text(
        -8.0,
        77.0,
        "r-projection only — NOT a 3-D annular/torus tag;\nproduction uses distance to each discrete nozzle axis",
        fontsize=8,
        color="#7f2704",
    )
    side.text(-150.0, 112.0, "bow-shock startup envelope: at least L1\nflow sensor may promote to L2", fontsize=8)

    for station, label in [(40.0, "40"), (70.0, "70"), (130.0, "130"), (body["aft_end_x"], "268.14")]:
        side.axvline(station, color="black", lw=0.7, ls=":")
        side.text(station + 2.0, 137.0, f"x={label} mm", rotation=90, va="top", fontsize=7)

    side.set_xlim(-180.0, 290.0)
    side.set_ylim(0.0, 145.0)
    side.set_xlabel("body-frame x from TSP [mm]; upstream is negative")
    side.set_ylabel("radial distance r [mm]")
    side.set_title("A. Proposed IBM Cartesian AMR zoning — meridional design mask")
    side.grid(alpha=0.12)
    colorbar = fig.colorbar(image, ax=side, orientation="horizontal", pad=0.06, shrink=0.78)
    colorbar.set_ticks([0, 1, 2, 3])
    colorbar.set_ticklabels(
        [
            f"L0  dx={dx_by_level[0]:.4f} mm",
            f"L1  dx={dx_by_level[1]:.4f} mm",
            f"L2  dx={dx_by_level[2]:.4f} mm",
            f"L3  dx={dx_by_level[3]:.4f} mm",
        ]
    )

    body_radius = float(parameters["body"]["maximum_radius_mm"])
    front.add_patch(plt.Circle((0.0, 0.0), body_radius, fc="#eeeeee", ec="black", lw=1.0))
    for throat in throats:
        outer = plt.Circle(
            (throat["y"], throat["z"]), 12.0, fc=colors[2], ec="none", alpha=0.7
        )
        inner = plt.Circle(
            (throat["y"], throat["z"]), 10.0, fc=colors[3], ec="black", lw=0.7, alpha=0.95
        )
        throat_circle = plt.Circle(
            (throat["y"], throat["z"]),
            throat["radius"],
            fc="white",
            ec="black",
            lw=0.8,
        )
        front.add_patch(outer)
        front.add_patch(inner)
        front.add_patch(throat_circle)
        front.text(throat["y"], throat["z"], f"N{throat['number']}", ha="center", va="center", fontsize=8)
    front.set_aspect("equal")
    front.set_xlim(-72.0, 72.0)
    front.set_ylim(-72.0, 72.0)
    front.set_xlabel("y [mm]")
    front.set_ylabel("z [mm]")
    front.set_title("B. Tri-nozzle transverse zones\n10 mm L3 core + 12 mm L2 envelope")
    front.grid(alpha=0.15)

    xline = np.linspace(body["nose_tip_x"], body["aft_end_x"], 800)
    target_level = np.select(
        [xline <= 40.0, xline <= 70.0, xline <= 130.0], [3, 2, 1], default=0
    )
    axial.step(xline, np.full_like(xline, 3), where="mid", color="#d73027", ls="--", lw=1.7, label="current default support seed")
    axial.step(xline, target_level, where="mid", color="#1a9850", lw=2.2, label="proposed baseline surface target")
    axial.fill_between(xline, 0, target_level, step="mid", color="#1a9850", alpha=0.12)
    axial.set_yticks([0, 1, 2, 3])
    axial.set_yticklabels(
        [
            f"L0\n{dx_by_level[0]:.3f} mm",
            f"L1\n{dx_by_level[1]:.3f} mm",
            f"L2\n{dx_by_level[2]:.3f} mm",
            f"L3\n{dx_by_level[3]:.3f} mm",
        ]
    )
    axial.set_xlim(body["nose_tip_x"], body["aft_end_x"])
    axial.set_ylim(-0.15, 3.25)
    axial.set_xlabel("surface station x from TSP [mm]")
    axial.set_ylabel("target AMR level")
    axial.set_title("C. Why a support-buffer level cap is required")
    axial.grid(alpha=0.2)
    axial.legend(loc="lower left", fontsize=8)
    axial.text(
        138.0,
        2.55,
        "Default ib.amr_support_buffer=1\nseeds every interface at every level.\nuser_tagging alone cannot make\nthe aft body coarse.",
        fontsize=8,
        color="#a50026",
        va="top",
    )
    axial.text(
        138.0,
        0.55,
        "Support-closure promotion still overrides\nthe baseline target where another tag\ncreates a child grid near IBM.",
        fontsize=7,
        color="#006837",
        va="top",
    )
    axial.text(
        138.0,
        1.30,
        "At equal physical volume:\nL0 has 1/512 of L3 cell density.",
        fontsize=8,
        color="#006837",
        va="top",
    )

    fig.suptitle(
        "Test 1853 IBM dx proposal — analytic target only; actual AMReX boxes must be checked from a regridded plotfile",
        fontsize=13,
    )
    fig.savefig(figures / "proposed_ibm_dx_distribution.png", dpi=220)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output = args.output.resolve()
    figures = output / "figures"
    visualization = output / "visualization"
    reports = output / "reports"
    figures.mkdir(parents=True, exist_ok=True)
    visualization.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)

    parameters = load_parameters()
    body = derived_body(parameters)
    throats = throat_data(output)
    mesh_path = output / "stl" / f"{MODEL}.stl"
    mesh = pv.read(mesh_path)
    if not isinstance(mesh, pv.PolyData):
        mesh = mesh.extract_surface()
    mesh = mesh.triangulate()
    arrays = face_arrays(mesh)
    zone = classify_zones(arrays, throats, body)
    attach_metrics(mesh, arrays, zone)

    vtp = visualization / f"{MODEL}_surface_mesh_metrics.vtp"
    mesh.save(vtp, binary=True)
    render_surface_overview(mesh, figures)

    single = pv.read(output / "stl" / "test1853_run165_single.stl")
    if not isinstance(single, pv.PolyData):
        single = single.extract_surface()
    single = single.triangulate()
    single_arrays = face_arrays(single)
    single.cell_data["hmax_mm"] = single_arrays["hmax"].astype(np.float32)
    single.cell_data["centroid_x_mm"] = single_arrays["centroid"][:, 0].astype(
        np.float32
    )
    render_configuration_comparison(single, mesh, figures)

    render_axial_profile(arrays, body, figures)
    render_dx_proposal(parameters, body, throats, figures)

    report = build_report(mesh, arrays, zone, output)
    report_path = reports / "mesh_distribution_visualization.json"
    report_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"PASS: wrote {vtp}")
    print(f"PASS: wrote visualization figures to {figures}")
    print(f"PASS: wrote {report_path}")


if __name__ == "__main__":
    main()
